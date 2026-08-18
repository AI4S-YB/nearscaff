"""Selective gap filling of nearscaff scaffolds with long reads.

Only gaps flanked by high-confidence contigs (tiers ``protein``/``asm5``
by default, from ``nearscaff_tiered.paf``) are eligible.  Short-read
approaches were benchmarked and removed: blind k-mer walks (abyss-sealer)
and short-read clip extension both closed 0/1001 eligible gaps on real
data — contig ends are tandem repeats.  Long reads fix this:

Recipe (validated on a HiFi benchmark):
  1. build a mini reference of the two ~1.5 kb flanks of every eligible gap
  2. map long reads to it with minimap2 -> span candidates (reads hitting
     BOTH flanks) + tails (clips extending into the gap from one flank)
  3. span fills: per-gap consensus (median-length candidate)
  4. tail overlap closure: ava-pb all-vs-all on tails, dovetail L x R
     overlaps; the overlap segment must map to a UNIQUE locus in the
     scaffolds (repeat-overlap pseudo-closures are the dominant failure
     mode and are filtered out)
"""
import gzip
import logging
import os
import subprocess

from nearscaff.agp import AGPSeqLine, AGPGapLine, AGPReader, AGPWriter

logger = logging.getLogger("nearscaff")

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def read_fasta(path: str) -> dict:
    """Read a FASTA file into {id: seq}; id = first token of the header."""
    seqs: dict[str, list] = {}
    name = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_tiers(tiered_paf_path: str) -> dict:
    """Read {contig: tier} from nearscaff_tiered.paf (nc:Z:<tier> tag)."""
    tiers: dict[str, str] = {}
    with open(tiered_paf_path) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            contig = fields[0]
            for tag in fields[12:]:
                if tag.startswith("nc:Z:"):
                    tiers[contig] = tag[5:]
                    break
    return tiers


def classify_gaps(agp_lines: list, tiers: dict,
                  allowed: set = frozenset({"protein", "asm5"})) -> set:
    """Return indices of AGP gap lines eligible for filling.

    A gap is eligible iff the nearest sequence components on both sides
    (within the same scaffold object) both have a tier in *allowed*.
    """
    eligible = set()
    for i, line in enumerate(agp_lines):
        if not isinstance(line, AGPGapLine):
            continue
        left = right = None
        for j in range(i - 1, -1, -1):
            if agp_lines[j].object_name != line.object_name:
                break
            if isinstance(agp_lines[j], AGPSeqLine):
                left = agp_lines[j].component_id
                break
        for j in range(i + 1, len(agp_lines)):
            if agp_lines[j].object_name != line.object_name:
                break
            if isinstance(agp_lines[j], AGPSeqLine):
                right = agp_lines[j].component_id
                break
        if (tiers.get(left) in allowed) and (tiers.get(right) in allowed):
            eligible.add(i)
    return eligible


def _component_seq(line: AGPSeqLine, contig_seqs: dict) -> str:
    """Component sequence in scaffold orientation (honours beg/end + strand)."""
    seq = contig_seqs[line.component_id][line.component_beg - 1:line.component_end]
    if line.orientation == "-":
        seq = _revcomp(seq)
    return seq


def write_outputs(agp_lines: list, eligible: set, fills: dict,
                  contig_seqs: dict,
                  out_agp: str, out_fasta: str,
                  pe_resize: dict | None = None,
                  trims: dict | None = None,
                  fill_prefix: str = "gapfill",
                  extensions: dict | None = None,
                  placements: dict | None = None) -> dict:
    """Write the gap-filled AGP + scaffold FASTA; return the report dict.

    *fills* is {agp_idx: fill_or_None}; gaps with a non-None fill are
    closed (empty fill = collapsed, the gap line is dropped).
    *pe_resize* is {agp_idx: estimated_size}: unfilled gaps whose size is
    estimated from PE spans are re-written as type-N gaps of the
    estimated length (evidence: paired-ends).
    *trims* is {agp_idx: overlap_bp}: end-join sites — the gap line is
    dropped and the preceding component is shortened by overlap_bp.
    *extensions* is {agp_idx: (left_seq, right_seq)}: one-sided fills —
    the gap stays open but its known content is written as W components
    adjacent to the flank(s) (left_seq hugs the left flank, right_seq
    hugs the right flank), with the residual gap (original length,
    placeholder) in between.
    *placements* is {agp_idx: [block, ...]}: ref-guided placements —
    the gap line is replaced by an ordered list of blocks, each either
    a sequence string (written as a W component, e.g. transcript exons)
    or a (None, est_len) tuple (estimated-N gap row of length est_len:
    intron/intergenic spacer, linkage "no").  Priority when several
    mechanisms target the same gap: fills > placements > extensions >
    pe_resize.
    All other lines pass through unchanged (coordinates and part numbers
    are re-numbered consistently).
    """
    # Map: agp gap index -> fill sequence (None stays a gap)
    fill_for_gap: dict[int, str] = {idx: seq for idx, seq in fills.items()
                                    if seq is not None}
    pe_resize = pe_resize or {}
    trims = trims or {}
    extensions = extensions or {}
    placements = placements or {}
    # preceding component index -> bp to trim off its 3' end
    trim_before: dict[int, int] = {}
    for gidx, ov in trims.items():
        j = gidx - 1
        while j >= 0 and not isinstance(agp_lines[j], AGPSeqLine):
            j -= 1
        if j >= 0:
            trim_before[j] = trim_before.get(j, 0) + ov

    fill_seqs: dict[str, str] = {}
    new_lines: list = []
    pos = 0
    part = 0
    cur = None
    n_closed = 0
    bases_filled = 0
    n_resized = 0
    n_endjoin = 0
    n_extended = 0
    n_ext_applied = 0
    n_placed = 0

    for idx, line in enumerate(agp_lines):
        if line.object_name != cur:
            cur = line.object_name
            pos = 1
            part = 0
        part += 1
        if isinstance(line, AGPGapLine) and idx in trims:
            part -= 1            # end-join: drop the gap line
            n_endjoin += 1
            continue
        if isinstance(line, AGPGapLine) and idx in fill_for_gap:
            seq = fill_for_gap[idx]
            fill_id = f"{cur}_{fill_prefix}{part}"
            fill_seqs[fill_id] = seq
            if len(seq) == 0:
                # Collapsed gap: emit nothing, just drop the gap line.
                part -= 1
                n_closed += 1
                continue
            new_lines.append(AGPSeqLine(
                cur, pos, pos + len(seq) - 1, part,
                "W", fill_id, 1, len(seq), "+"))
            pos += len(seq)
            n_closed += 1
            bases_filled += len(seq)
        elif isinstance(line, AGPGapLine):
            if idx in placements:
                blocks = placements[idx]
                n_placed += 1
                for bi, item in enumerate(blocks):
                    if bi > 0:
                        part += 1
                    if isinstance(item, tuple):
                        est = item[1]
                        new_lines.append(AGPGapLine(
                            cur, pos, pos + est - 1, part,
                            "N", est, line.gap_type, "no", "na"))
                        pos += est
                    else:
                        bseq = item
                        fill_id = f"{cur}_{fill_prefix}{part}"
                        fill_seqs[fill_id] = bseq
                        new_lines.append(AGPSeqLine(
                            cur, pos, pos + len(bseq) - 1, part,
                            "W", fill_id, 1, len(bseq), "+"))
                        pos += len(bseq)
                        bases_filled += len(bseq)
            elif idx in extensions:
                n_ext_applied += 1
                lseq, rseq = extensions[idx]
                if lseq:
                    part_l = part
                    fill_id = f"{cur}_{fill_prefix}{part_l}"
                    fill_seqs[fill_id] = lseq
                    new_lines.append(AGPSeqLine(
                        cur, pos, pos + len(lseq) - 1, part_l,
                        "W", fill_id, 1, len(lseq), "+"))
                    pos += len(lseq)
                    bases_filled += len(lseq)
                    n_extended += 1
                    part += 1
                new_lines.append(AGPGapLine(
                    cur, pos, pos + line.gap_length - 1, part,
                    line.component_type, line.gap_length, line.gap_type,
                    line.linkage, line.linkage_evidence))
                pos += line.gap_length
                if rseq:
                    part += 1
                    fill_id = f"{cur}_{fill_prefix}{part}"
                    fill_seqs[fill_id] = rseq
                    new_lines.append(AGPSeqLine(
                        cur, pos, pos + len(rseq) - 1, part,
                        "W", fill_id, 1, len(rseq), "+"))
                    pos += len(rseq)
                    bases_filled += len(rseq)
                    n_extended += 1
            elif idx in pe_resize:
                est = pe_resize[idx]
                new_lines.append(AGPGapLine(
                    cur, pos, pos + est - 1, part,
                    "N", est, line.gap_type, line.linkage, "paired-ends"))
                pos += est
                n_resized += 1
            else:
                new_lines.append(AGPGapLine(
                    cur, pos, pos + line.gap_length - 1, part,
                    line.component_type, line.gap_length, line.gap_type,
                    line.linkage, line.linkage_evidence))
                pos += line.gap_length
        else:
            trim = trim_before.get(idx, 0)
            span = line.object_end - line.object_beg + 1 - trim
            # the trimmed end is the scaffold-side 3' end of the component;
            # for a "-" component that is the contig's 5' side (beg)
            if line.orientation == "-":
                cbeg, cend = line.component_beg + trim, line.component_end
            else:
                cbeg, cend = line.component_beg, line.component_end - trim
            new_lines.append(AGPSeqLine(
                cur, pos, pos + span - 1, part,
                line.component_type, line.component_id,
                cbeg, cend, line.orientation))
            pos += span

    with open(out_agp, "w") as f:
        f.write(AGPWriter().format(new_lines))

    with open(out_fasta, "w") as f:
        cur = None
        buf: list[str] = []
        for line in new_lines:
            if line.object_name != cur:
                if cur is not None:
                    _write_fasta_record(f, cur, buf)
                cur = line.object_name
                buf = []
            if isinstance(line, AGPSeqLine):
                if line.component_id in fill_seqs:
                    buf.append(fill_seqs[line.component_id])
                else:
                    buf.append(_component_seq(line, contig_seqs))
            else:
                buf.append("N" * line.gap_length)
        if cur is not None:
            _write_fasta_record(f, cur, buf)

    report = {
        "gaps_total": sum(1 for l in agp_lines if isinstance(l, AGPGapLine)),
        "gaps_eligible": len(eligible),
        "gaps_closed": n_closed,
        "gaps_resized": n_resized,
        "gaps_endjoined": n_endjoin,
        "gaps_extended": n_ext_applied,
        "gaps_placed": n_placed,
        "bases_filled": bases_filled,
    }
    return report


def _write_fasta_record(f, name: str, chunks: list) -> None:
    seq = "".join(chunks)
    f.write(f">{name}\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i + 60] + "\n")


# ---------------------------------------------------------------------------
# Long-read recruitment + closure
# ---------------------------------------------------------------------------

LR_FLANK = 1500
LR_MIN_CLIP = 300
LR_EDGE = 10
LR_MAX_TAILS = 30
LR_MIN_OVLP = 500
LR_MIN_OVLP_IDENT = 0.85


def scaffold_sequences(agp_lines: list, contig_seqs: dict) -> dict:
    """Reconstruct {scaffold: sequence} from AGP + contigs (gaps as Ns,
    orientation-aware, uppercased)."""
    scafs: dict[str, list] = {}
    for line in agp_lines:
        if isinstance(line, AGPSeqLine):
            scafs.setdefault(line.object_name, []).append(
                _component_seq(line, contig_seqs))
        elif isinstance(line, AGPGapLine):
            scafs.setdefault(line.object_name, []).append("N" * line.gap_length)
    return {k: "".join(v).upper() for k, v in scafs.items()}


def build_flank_fasta(agp_lines: list, eligible: set, scaf_seqs: dict,
                      flank: int = LR_FLANK, max_n_frac: float = 0.2) -> tuple:
    """Per eligible gap emit two flank records from the scaffold sequence.

    Returns (records, flank_len) where records = [(name, seq)] with names
    "gid|L" / "gid|R" and flank_len = {name: length}.  Flanks containing
    too many Ns (inside other unfilled gaps) are dropped.
    """
    records: list[tuple[str, str]] = []
    flank_len: dict[str, int] = {}
    for idx, line in enumerate(agp_lines):
        if not isinstance(line, AGPGapLine) or idx not in eligible:
            continue
        gid = f"{line.object_name}:{line.object_beg}-{line.object_end}"
        s = scaf_seqs[line.object_name]
        gb, ge = line.object_beg - 1, line.object_end
        for side, seq in (("L", s[max(0, gb - flank):gb]),
                          ("R", s[ge:ge + flank])):
            if not seq:
                continue
            n_frac = seq.count("N") / len(seq)
            if n_frac > max_n_frac:
                continue
            name = f"{gid}|{side}"
            records.append((name, seq))
            flank_len[name] = len(seq)
    return records, flank_len


def _write_records_fa(records: list, path: str) -> None:
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def run_minimap2(target_fa: str, reads: list, out_paf: str, preset: str,
                 threads: int = 4, extra: list | None = None,
                 minimap2: str = "minimap2") -> str:
    cmd = [minimap2, "-x", preset, "-t", str(threads)]
    if extra:
        cmd += list(extra)
    cmd += [target_fa] + list(reads)
    logger.info("Running: %s", " ".join(cmd))
    try:
        with open(out_paf, "w") as out:
            result = subprocess.run(cmd, stdout=out,
                                    stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        raise RuntimeError("minimap2 not found on PATH")
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed: {result.stderr.strip()[-1000:]}")
    return out_paf


def _fetch_read_seqs(fastqs, wanted: set) -> dict:
    """Fetch {name: seq} for *wanted* read names from FASTQ/FASTA(.gz)
    files (format auto-detected per file from the first record)."""
    seqs: dict[str, str] = {}
    if isinstance(fastqs, str):
        fastqs = [fastqs]
    for fq in fastqs:
        if not wanted - set(seqs):
            break
        opener = gzip.open if fq.endswith(".gz") else open
        with opener(fq, "rt") as f:
            h = f.readline()
            if h.startswith(">"):                # FASTA (multi-line)
                name = h[1:].split()[0]
                chunks: list[str] = []
                for line in f:
                    if line.startswith(">"):
                        if name in wanted:
                            seqs[name] = "".join(chunks)
                        name = line[1:].split()[0]
                        chunks = []
                    else:
                        chunks.append(line.strip())
                if name in wanted:
                    seqs[name] = "".join(chunks)
            else:                                # FASTQ (4-line)
                while h:
                    if not h.strip():
                        h = f.readline()
                        continue
                    name = h[1:].split()[0]
                    seq = f.readline().strip()
                    f.readline()
                    f.readline()
                    if name in wanted:
                        seqs[name] = seq
                    h = f.readline()
    return seqs


def parse_span_fills(paf_path: str, reads_fq, min_reads: int = 1,
                     max_shortfall: int = 300,
                     max_depth_factor: float | None = 3.0) -> dict:
    """Extract fill sequences from reads hitting BOTH flanks of a gap.

    A spanning read has two alignments (to gid|L and gid|R).  Each
    alignment must reach within *max_shortfall* bp of the gap-side edge
    of its flank; the read coordinate of the gap boundary is obtained by
    extrapolating the alignment to the edge (the shortfall is flank
    sequence already present in the scaffold and is NOT included in the
    fill).  The segment between the two boundaries (genomic forward
    orientation) is a fill candidate; per gap the median-length candidate
    wins.  Returns {gid: fill_seq}.

    Depth gate: the per-gap candidate count is the spanning depth.  At
    unique loci it approximates the read coverage; repeat gaps attract
    reads from the whole family and show inflated depth.  When
    *max_depth_factor* is set, gaps whose depth exceeds
    factor x (median depth over gaps) are left open (over-collapsed);
    repeat gaps within a reasonable depth are closed with the
    median-length (family consensus) fill.
    """
    flank_len: dict[str, int] = {}
    hits: dict[tuple, dict] = {}
    with open(paf_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[2] == "*" or p[5] == "*":
                continue
            gid, side = p[5].rsplit("|", 1)
            flank_len[p[5]] = int(p[6])
            ts, te = int(p[7]), int(p[8])
            if side == "L":
                if te < flank_len[p[5]] - max_shortfall:
                    continue    # alignment stops far short of the gap edge
            else:
                if ts > max_shortfall:
                    continue
            key = (gid, p[0])
            h = hits.setdefault(key, {})
            if side not in h:
                # store strand, qs, qe, ts, te
                h[side] = (p[4], int(p[2]), int(p[3]), ts, te)
    seqs = _fetch_read_seqs(reads_fq, {qn for _g, qn in hits})

    cands: dict[str, list] = {}
    for (gid, qn), h in hits.items():
        if set(h) != {"L", "R"}:
            continue
        s = seqs.get(qn)
        if s is None:
            continue
        strand_l, qs_l, qe_l, _ts_l, te_l = h["L"]
        strand_r, qs_r, qe_r, ts_r, _te_r = h["R"]
        if strand_l != strand_r:
            continue
        len_l = flank_len[f"{gid}|L"]
        if strand_l == "+":
            # gap boundaries in read coords (extrapolated to flank edges)
            b0 = qe_l + (len_l - te_l)   # gap start
            b1 = qs_r - ts_r             # gap end
            if b1 <= b0 or b0 < 0 or b1 > len(s):
                continue
            seg = s[b0:b1]
        else:
            b0 = qe_r + ts_r             # gap end (read coords)
            b1 = qs_l - (len_l - te_l)   # gap start
            if b1 <= b0 or b0 < 0 or b1 > len(s):
                continue
            seg = _revcomp(s[b0:b1])
        if "N" not in seg.upper():
            cands.setdefault(gid, []).append(seg)

    # depth gate
    skip: set = set()
    if max_depth_factor and cands:
        depths = sorted(len(v) for v in cands.values())
        baseline = depths[(len(depths) - 1) // 2] or 1
        cutoff = max_depth_factor * max(baseline, min_reads)
        for gid, segs in cands.items():
            if len(segs) > cutoff:
                skip.add(gid)
        if skip:
            logger.info("Span depth gate: baseline %d reads, rejecting %d "
                        "over-collapsed gaps (>%.1fx)", baseline, len(skip),
                        max_depth_factor)

    fills: dict[str, str] = {}
    for gid, segs in cands.items():
        if gid in skip or len(segs) < min_reads:
            continue
        segs.sort(key=len)
        fills[gid] = segs[len(segs) // 2]
    return fills


def parse_tails(paf_path: str, flank_len: dict, reads_fq,
                min_clip: int = LR_MIN_CLIP, edge: int = LR_EDGE,
                max_tails: int = LR_MAX_TAILS) -> dict:
    """Collect clip tails extending into gaps from flank-edge alignments.

    Returns {gid: {"L": [tail], "R": [tail]}} in genomic fwd orientation.
    """
    events: list = []
    with open(paf_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[2] == "*" or p[5] == "*":
                continue
            qn, ql, qs, qe, strand = p[0], int(p[1]), int(p[2]), int(p[3]), p[4]
            tn, ts, te = p[5], int(p[7]), int(p[8])
            gid, side = tn.rsplit("|", 1)
            L = flank_len[tn]
            if side == "L" and abs(te - L) <= edge:
                if strand == "+" and ql - qe >= min_clip:
                    events.append((gid, "L", qn, qe, "+"))
                elif strand == "-" and qs >= min_clip:
                    events.append((gid, "L", qn, qs, "-"))
            elif side == "R" and ts <= edge:
                if strand == "+" and qs >= min_clip:
                    events.append((gid, "R", qn, qs, "+"))
                elif strand == "-" and ql - qe >= min_clip:
                    events.append((gid, "R", qn, qe, "-"))

    seqs = _fetch_read_seqs(reads_fq, {qn for _g, _s, qn, _p, _st in events})

    tails: dict[str, dict] = {}
    for gid, side, qn, pos, strand in events:
        s = seqs.get(qn)
        if s is None:
            continue
        if side == "L":
            tail = s[pos:] if strand == "+" else _revcomp(s[:pos])
        else:
            # R tails are stored deep-in-gap -> flank; the left clip of a
            # "+" read is already in genomic forward orientation
            tail = s[:pos] if strand == "+" else _revcomp(s[pos:])
        if len(tail) < min_clip:
            continue
        t = tails.setdefault(gid, {"L": [], "R": []})
        if len(t[side]) < max_tails:
            t[side].append(tail)
    return tails


def overlap_closure_fills(tails: dict, scaf_seqs: dict, work_prefix: str,
                          threads: int = 4, min_ovlp: int = LR_MIN_OVLP,
                          min_ident: float = LR_MIN_OVLP_IDENT,
                          unique_check: bool = True,
                          uniq_preset: str = "map-hifi",
                          minimap2: str = "minimap2") -> dict:
    """Close gaps by L x R tail overlap.

    Only dovetail (L-suffix x R-prefix) overlaps are used:
    fill = L_tail[:qe] + R_tail[te:].  With *unique_check* (default, used
    for long reads) the overlap segment must additionally map to a unique
    locus in the scaffolds, filtering repeat-driven pseudo-closures; the
    uniqueness mapping uses *uniq_preset* ("splice" for cDNA tails, whose
    sequence is not contiguous in the genome).
    Short-read mode disables it (closure rate over correctness).
    Returns {gid: fill_seq}.
    """
    records = []
    for gid in sorted(tails):
        for side in ("L", "R"):
            for j, t in enumerate(tails[gid][side]):
                records.append((f"{gid}|{side}{j}", t))
    if not records:
        return {}
    tails_fa = work_prefix + ".tails.fa"
    _write_records_fa(records, tails_fa)
    tail_seqs = dict(records)

    ava_paf = work_prefix + ".tails.ava.paf"
    cmd = [minimap2, "-x", "ava-pb", "-t", str(threads), tails_fa, tails_fa]
    logger.info("Running: %s", " ".join(cmd))
    with open(ava_paf, "w") as out:
        result = subprocess.run(cmd, stdout=out,
                                stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        logger.warning("ava-pb failed (%s); no overlap closures",
                       result.stderr.strip()[-500:])
        return {}

    # best dovetail overlap per gap
    best: dict[str, tuple] = {}
    with open(ava_paf) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[2] == "*" or p[5] == "*":
                continue
            q, t = p[0], p[5]
            gid_q, side_q = q.rsplit("|", 1)
            gid_t, side_t = t.rsplit("|", 1)
            if gid_q != gid_t or side_q == side_t:
                continue
            if not (side_q.startswith("L") and side_t.startswith("R")):
                continue
            qe, ts_, te = int(p[3]), int(p[7]), int(p[8])
            ovlen = int(p[10])
            ident = int(p[9]) / max(ovlen, 1)
            # dovetail: overlap reaches L tail's 3' end and R tail's 5' end
            if abs(qe - len(tail_seqs[q])) > LR_EDGE or ts_ > LR_EDGE:
                continue
            if ovlen < min_ovlp or ident < min_ident:
                continue
            if gid_q not in best or ovlen > best[gid_q][0]:
                best[gid_q] = (ovlen, q, t, qe, ts_, te)

    if not best:
        return {}
    if not unique_check:
        return {gid: tail_seqs[q][:qe] + tail_seqs[t][te:]
                for gid, (ovlen, q, t, qe, ts_, te) in best.items()}

    # uniqueness of the overlap segment in the scaffolds
    ov_fa = work_prefix + ".overlap.fa"
    with open(ov_fa, "w") as out:
        for gid, (ovlen, q, t, qe, ts_, te) in best.items():
            out.write(f">{gid}\n{tail_seqs[q][qe - ovlen:qe]}\n")
    scaf_fa = work_prefix + ".scaffolds.fa"
    _write_records_fa(sorted(scaf_seqs.items()), scaf_fa)
    uniq_paf = work_prefix + ".overlap.uniq.paf"
    cmd = [minimap2, "-x", uniq_preset,
           "-N", "100", "-t", str(threads), scaf_fa, ov_fa]
    try:
        with open(uniq_paf, "w") as out:
            result = subprocess.run(cmd, stdout=out,
                                    stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            # fall back to a widely-supported preset if the local
            # minimap2 lacks the requested one
            cmd[2] = "map-ont" if uniq_preset != "map-ont" else "splice"
            with open(uniq_paf, "w") as out:
                subprocess.run(cmd, stdout=out,
                               stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        raise RuntimeError("minimap2 not found on PATH")

    # uniqueness by LOCUS (same-locus sub-alignments/secondaries don't
    # count against; repeat overlaps land on many distinct loci)
    loci: dict[str, set] = {}
    with open(uniq_paf) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[5] == "*":
                continue
            gid = p[0]
            loci.setdefault(gid, set()).add((p[5], int(p[7]) // 500))

    fills: dict[str, str] = {}
    for gid, (ovlen, q, t, qe, ts_, te) in best.items():
        n_loci = len(loci.get(gid, set()))
        if n_loci != 1:
            logger.debug("Overlap closure %s rejected: maps to %d loci "
                         "(repeat overlap)", gid, n_loci)
            continue
        fills[gid] = tail_seqs[q][:qe] + tail_seqs[t][te:]
    return fills


def parse_pe_spans(paf1: str, paf2: str, insert: int = 500,
                   max_shortfall: int = 300, min_pairs: int = 2) -> dict:
    """Estimate gap sizes from paired-end spans (short-read mode).

    A proper FR pair (R1 '+' on the L flank, R2 '-' on the R flank, both
    primary, mapq >= 20, within *max_shortfall* of the gap edge) implies
    gap size ~ insert - flank_portion_R1 - flank_portion_R2.  Per gap the
    median over >= min_pairs pairs is returned.  NOTE: short-read PE
    spans at repeat-rich junctions have a HIGH false-positive rate —
    use for gap resizing only, never for sequence fills.

    Returns {gid: est_size} (clamped to >= 1).
    """
    def load(paf):
        best: dict[str, dict] = {}   # read -> {gid|side: (strand, ts, te, tlen)}
        with open(paf) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 12 or p[2] == "*" or p[5] == "*":
                    continue
                if "tp:A:P" not in p[12:] or int(p[11]) < 20:
                    continue
                best.setdefault(p[0], {})[p[5]] = (
                    p[4], int(p[7]), int(p[8]), int(p[6]))
        return best

    m1, m2 = load(paf1), load(paf2)
    est: dict[str, list] = {}
    for qn, a1 in m1.items():
        a2 = m2.get(qn)
        if not a2:
            continue
        for t1, (s1, ts1, te1, L1) in a1.items():
            gid, side1 = t1.rsplit("|", 1)
            for t2, (s2, ts2, te2, L2) in a2.items():
                g2, side2 = t2.rsplit("|", 1)
                if g2 != gid or side2 == side1:
                    continue
                if side1 == "L":
                    if not (s1 == "+" and s2 == "-"):
                        continue
                    short1, short2 = L1 - te1, ts2
                    portion1 = (te1 - ts1) + short1
                    portion2 = (te2 - ts2) + short2
                else:
                    if not (s1 == "-" and s2 == "+"):
                        continue
                    short1, short2 = ts1, L2 - te2
                    portion1 = (te1 - ts1) + short1
                    portion2 = (te2 - ts2) + short2
                if short1 > max_shortfall or short2 > max_shortfall:
                    continue
                gap_est = insert - portion1 - portion2
                if gap_est > insert:   # inconsistent pair, skip
                    continue
                est.setdefault(gid, []).append(max(gap_est, 1))
    out: dict[str, int] = {}
    for gid, sizes in est.items():
        if len(sizes) >= min_pairs:
            sizes.sort()
            out[gid] = sizes[len(sizes) // 2]
    return out


def parse_endjoins(paf_path: str, min_reads: int = 2,
                   max_shortfall: int = 300, abut_window: int = 30,
                   min_trim: int = 10) -> dict:
    """Detect end-join sites: gaps whose flanking contigs actually overlap.

    For reads hitting BOTH flanks (same logic as span fills), the
    extrapolated gap size can be <= 0: the two contigs overlap.  Returns
    {gid: overlap_bp} where overlap_bp <= abut_window means the contigs
    abut (join with an empty fill) and larger values mean a real overlap
    (trim overlap_bp from the left contig end before joining, after the
    end-sequence identity check in run_gapfill).
    """
    flank_len: dict[str, int] = {}
    hits: dict[tuple, dict] = {}
    with open(paf_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[2] == "*" or p[5] == "*":
                continue
            gid, side = p[5].rsplit("|", 1)
            flank_len[p[5]] = int(p[6])
            ts, te = int(p[7]), int(p[8])
            if side == "L":
                if te < flank_len[p[5]] - max_shortfall:
                    continue
            else:
                if ts > max_shortfall:
                    continue
            key = (gid, p[0])
            h = hits.setdefault(key, {})
            if side not in h:
                h[side] = (p[4], int(p[2]), int(p[3]), ts, te)

    ests: dict[str, list] = {}
    for (gid, qn), h in hits.items():
        if set(h) != {"L", "R"}:
            continue
        strand_l, qs_l, qe_l, _ts_l, te_l = h["L"]
        strand_r, qs_r, qe_r, ts_r, _te_r = h["R"]
        if strand_l != strand_r:
            continue
        len_l = flank_len[f"{gid}|L"]
        if strand_l == "+":
            est = (qs_r - ts_r) - (qe_l + (len_l - te_l))
        else:
            est = (qs_l - (len_l - te_l)) - (qe_r + ts_r)
        ests.setdefault(gid, []).append(est)

    out: dict[str, int] = {}
    for gid, vals in ests.items():
        if len(vals) < min_reads:
            continue
        vals.sort()
        med = vals[len(vals) // 2]
        if med < -abut_window and -med >= min_trim:
            out[gid] = -med                   # overlap: trim this many bp
        elif med <= abut_window:
            out[gid] = 0                      # abut: empty fill
    return out


def _end_identity(scaf_seq: str, gb: int, ge: int, overlap: int,
                  k: int = 15) -> float:
    """k-mer containment of the right contig start in the left contig end
    (fraction of right-end k-mers present in the left end)."""
    left = scaf_seq[max(0, gb - overlap):gb]
    right = scaf_seq[ge:ge + overlap]
    if not left or not right:
        return 0.0
    kmers = {left[i:i + k] for i in range(len(left) - k + 1)}
    n = 0
    hit = 0
    for i in range(len(right) - k + 1):
        n += 1
        if right[i:i + k] in kmers:
            hit += 1
    return hit / max(n, 1)


def run_gapfill(agp_path: str, query_fasta: str, tiered_paf: str,
                output_dir: str, reads: list,
                method: str = "lr",
                lr_preset: str = "map-hifi",
                threads: int = 4,
                fill_tiers: set = frozenset({"protein", "asm5"}),
                min_span_reads: int = 1,
                overlap_closure: bool = True,
                sr_insert: int = 500,
                max_depth_factor: float | None = 3.0,
                minimap2: str = "minimap2",
                recruit_extra: list | None = None) -> dict:
    """Gap filling via flank mini-reference recruitment.

    method "lr": long reads (HiFi/ONT) — span fills + tail-overlap fills.
    method "sr": short reads — single-read span fills + PE-span gap
                 resizing (PE spans only ESTIMATE gap sizes; false-positive
                 rate at repeat-rich junctions is high, use with caution).
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(agp_path) as f:
        agp_lines = AGPReader().parse(f.read())
    tiers = read_tiers(tiered_paf)
    eligible = classify_gaps(agp_lines, tiers, allowed=fill_tiers)
    n_gaps = sum(1 for l in agp_lines if isinstance(l, AGPGapLine))
    logger.info("%d gaps total, %d eligible for filling (flanks in %s)",
                n_gaps, len(eligible), "/".join(sorted(fill_tiers)))

    contig_seqs = read_fasta(query_fasta)
    scaf_seqs = scaffold_sequences(agp_lines, contig_seqs)
    records, flank_len = build_flank_fasta(agp_lines, eligible, scaf_seqs)
    flanks_fa = os.path.join(output_dir, "gapfill.flanks.fa")
    _write_records_fa(records, flanks_fa)
    logger.info("Flank mini-reference: %d records -> %s", len(records), flanks_fa)

    if method == "sr":
        if len(reads) < 2:
            raise ValueError("short-read mode requires -1 and -2 (paired-end)")
        paf1 = os.path.join(output_dir, "gapfill.flanks.r1.paf")
        paf2 = os.path.join(output_dir, "gapfill.flanks.r2.paf")
        extra = recruit_extra or ["-N", "20"]
        run_minimap2(flanks_fa, [reads[0]], paf1, "sr", threads=threads,
                     extra=extra, minimap2=minimap2)
        run_minimap2(flanks_fa, [reads[1]], paf2, "sr", threads=threads,
                     extra=extra, minimap2=minimap2)
        span = parse_span_fills(paf1, [reads[0]], min_reads=min_span_reads,
                                max_depth_factor=max_depth_factor)
        span2 = parse_span_fills(paf2, [reads[1]], min_reads=min_span_reads,
                                 max_depth_factor=max_depth_factor)
        for gid, seq in span2.items():
            span.setdefault(gid, seq)
        logger.info("Single-read span closures: %d gaps", len(span))
        pe_est = parse_pe_spans(paf1, paf2, insert=sr_insert)
        logger.info("PE-span size estimates (>=2 pairs): %d gaps "
                    "(CAUTION: high FP rate at repeat junctions)", len(pe_est))
        # SR tail-overlap closure: merge mate tails; dovetail merge gives
        # full gap sequence for gaps up to ~2x read length.  Uniqueness
        # check disabled — closure rate over correctness (documented).
        tails: dict = {}
        for pafi, fqi in ((paf1, reads[0]), (paf2, reads[1])):
            for gid, sides in parse_tails(pafi, flank_len, [fqi],
                                          min_clip=30, max_tails=50).items():
                cur = tails.setdefault(gid, {"L": [], "R": []})
                for side in ("L", "R"):
                    cur[side].extend(sides[side])
        prefix = os.path.join(output_dir, "gapfill")
        ovl = overlap_closure_fills(tails, scaf_seqs, prefix,
                                    threads=threads, min_ovlp=25,
                                    min_ident=0.85, unique_check=False,
                                    minimap2=minimap2)
        logger.info("Tail-overlap closures (no uniqueness filter, "
                    "CAUTION: high FP rate): %d gaps", len(ovl))
    else:
        paf = os.path.join(output_dir, "gapfill.flanks.paf")
        run_minimap2(flanks_fa, reads, paf, lr_preset, threads=threads,
                     extra=recruit_extra or ["-N", "20"], minimap2=minimap2)
        span = parse_span_fills(paf, reads, min_reads=min_span_reads,
                                max_depth_factor=max_depth_factor)
        logger.info("Span closures: %d gaps", len(span))
        pe_est = {}
        ovl = {}
        if overlap_closure:
            tails = parse_tails(paf, flank_len, reads)
            prefix = os.path.join(output_dir, "gapfill")
            ovl = overlap_closure_fills(tails, scaf_seqs, prefix,
                                        threads=threads, minimap2=minimap2)
            logger.info("Tail-overlap closures (unique overlap): %d gaps",
                        len(ovl))

    # ---- end-join detection (contigs that actually overlap) ----
    endjoin: dict = {}
    for pafi in ([paf1, paf2] if method == "sr" else [paf]):
        for gid, ov in parse_endjoins(pafi).items():
            endjoin[gid] = max(endjoin.get(gid, 0), ov)
    trims: dict = {}          # gid -> overlap bp (left contig end trimmed)
    n_abut = 0
    for gid, ov in endjoin.items():
        scaf, rng = gid.rsplit(":", 1)
        beg_s, end_s = rng.split("-")
        s = scaf_seqs[scaf]
        gb, ge = int(beg_s) - 1, int(end_s)
        if ov == 0:
            n_abut += 1
        elif _end_identity(s, gb, ge, ov) >= 0.8:
            trims[gid] = ov
    logger.info("End-join sites: %d abut, %d overlap-trim "
                "(identity >= 80%%)", n_abut, len(trims))

    fills = dict(ovl)
    fills.update(span)   # span fills win on conflict
    for gid in [g for g, ov in endjoin.items() if ov == 0]:
        fills.setdefault(gid, "")   # abut: collapse the gap
    # gid -> agp index
    gid_to_idx = {f"{l.object_name}:{l.object_beg}-{l.object_end}": i
                  for i, l in enumerate(agp_lines)
                  if isinstance(l, AGPGapLine)}
    fills_idx = {gid_to_idx[g]: s for g, s in fills.items() if g in gid_to_idx}
    trims_idx = {gid_to_idx[g]: v for g, v in trims.items()
                 if g in gid_to_idx and g not in fills}

    out_agp = os.path.join(output_dir, "nearscaff.gapfill.agp")
    out_fasta = os.path.join(output_dir, "nearscaff.gapfill.scaffolds.fa")
    report = write_outputs(agp_lines, eligible, fills_idx, contig_seqs,
                           out_agp, out_fasta,
                           pe_resize={gid_to_idx[g]: v for g, v in pe_est.items()
                                      if g in gid_to_idx},
                           trims=trims_idx)
    report["span_closed"] = len(span)
    report["overlap_closed"] = len(ovl)
    report["pe_resized"] = len(pe_est)
    report["endjoin_abut"] = n_abut
    report["endjoin_trim"] = len(trims_idx)
    logger.info("Gap filling complete: %d/%d eligible gaps closed "
                "(%d bp filled; %d gaps resized by PE spans; "
                "%d end-joins)", report["gaps_closed"],
                report["gaps_eligible"], report["bases_filled"],
                report.get("gaps_resized", 0), report.get("gaps_endjoined", 0))
    logger.info("Gap-filled AGP: %s", out_agp)
    logger.info("Gap-filled FASTA: %s", out_fasta)
    return report
