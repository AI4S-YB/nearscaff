"""Transcriptome-guided filling of genic scaffold gaps (rna-fill).

Scaffolding raises BUSCO but genic regions remain full of Ns, so
downstream annotation barely improves.  This module maps transcript
sequences (Iso-Seq / ONT cDNA / assembled transcripts; short-read
RNA-seq must be assembled first) back onto the scaffolds and fills the
gaps that transcripts can cross.

Only gaps flanked by high-confidence contigs (tiers ``protein``/``asm5``
by default) are eligible — same rule as ``gapfill``.  The transcript
evidence itself restricts filling to genic regions.

Recipe (reuses the gapfill flank mini-reference machinery):
  1. build a mini reference of the two ~2 kb flanks of every eligible gap
  2. map transcripts to it with minimap2 (splice preset)
  3. span fills: a transcript hitting BOTH flanks contributes its middle
     cDNA segment.  A near-zero middle means the two flanks are adjacent
     exons — the gap is (mostly) intronic.
  4. tail-overlap closure: dovetail L x R clips, with the uniqueness
     check mapped back with a splice preset (cDNA tails are not
     contiguous in the genome)

Fill modes:
  - ``cdna`` (default): fill whatever transcripts span; if the gap
    contains introns they are compressed out (the fill is cDNA).  Abut
    closures (near-zero middle) join the flanks directly.
  - ``exon-only``: only fill gaps believed to lie inside a single exon —
    gaps whose visible edges look like canonical splice sites (GT..AG on
    either strand) or whose transcripts show a near-zero middle are
    skipped.  Recommended when ``--internal-n`` is on (component-internal
    N runs usually cover introns; cDNA fills would compress them out and
    distort downstream ab initio annotation).

Read recruitment (``--reads`` + ``--proteins``):
  When raw RNA reads and reference proteins are given instead of
  assembled transcripts, rna-fill first runs a recruitment phase:
  miniprot maps the proteins onto the scaffolds, mRNAs whose translation
  contains X (CDS crosses an N run) become "broken gene" bait loci
  (± --bait-pad), and together with the gap flanks they form a combined
  bait mini-reference.  Reads mapping to the bait (and their mates) are
  extracted for targeted de novo assembly.  Assemble the recruited reads
  with any transcript assembler, then re-run rna-fill with ``-T``.

Fills are cDNA sequence by construction: in ``cdna`` mode the inserted
sequence is real exonic sequence but intron content inside the gap is
lost.  Use the output for annotation completeness, not for analyses
that need true intronic sequence.

Ref-guided whole-gene placement (``--ref``):
  Gaps still unfilled after every flank-based mechanism get one more
  pass.  The gap is bracketed on the reference genome through its
  flanking contigs' tiered-PAF alignments; reference genes lying fully
  inside the bracket are necessarily missing from the query (the region
  is all N), so transcripts assigned to them (by mapping the
  transcripts against the reference) are written into the gap in
  reference order: exon components of real transcript sequence
  separated by intron/intergenic spacers whose lengths come from the
  reference (clamped).  Intron placeholders carry the canonical splice
  consensus at their edges (GT..AG / CT..AC by strand — splicing is
  proven by the transcript's CIGAR, only the intronic bases are
  unknown), so miniprot and ab initio annotators can splice across
  them; intergenic and bracket-edge spacers are plain N gap rows.  Gene coordinates come from ``--ref-gff``
  when given, else from miniprot mapping ``--proteins`` back onto the
  reference.  Unlike cdna fills, placed genes keep their intron
  structure as N placeholders, so ab initio annotation is not distorted
  — but the intronic sequence itself is unknown, only estimated in
  length.
"""
import logging
import os
import re

from nearscaff.agp import AGPSeqLine, AGPGapLine, AGPReader
from nearscaff.gapfill import (read_fasta, read_tiers, classify_gaps,
                               scaffold_sequences, build_flank_fasta,
                               run_minimap2, _fetch_read_seqs, parse_tails,
                               overlap_closure_fills, write_outputs,
                               _write_records_fa, _revcomp)

logger = logging.getLogger("nearscaff")

TX_FLANK = 2000
TX_MAX_SHORTFALL = 300
TX_ABUT_WINDOW = 30
TX_MIN_OVLP = 300
TX_MIN_OVLP_IDENT = 0.85
TX_MIN_EXT = 500        # one-sided extension: min tail into the gap
TX_MAX_EXT = 20000      # ... capped at this length (deep end truncated)

PLACE_MIN_SPACER = 20   # ref-derived intron spacer: min estimated length
PLACE_MAX_SPACER = 50000  # ... capped (ref distances are only estimates)
PLACE_MIN_EDGE = 50     # bracket-edge/intergenic spacers: min length
PLACE_MAX_BRACKET = 1_000_000
PLACE_MAX_GENES = 50
PLACE_MIN_OVLP = 0.5
PLACE_INTACT_COV = 0.8  # proteins with an intact query hit (cov >= this,
                        # no X) are not placed again (anti-duplication)


_N_RUN = re.compile(r"[Nn]+")


def explode_internal_ns(agp_lines: list, contig_seqs: dict,
                        min_len: int = 20, min_edge: int = 50) -> tuple:
    """Split components at internal N runs into sub-rows + gap rows.

    Short-read assemblies carry N runs *inside* their scaffolds (the
    assembler's estimated spacers); when such a run lands in a CDS the
    transcript/annotation is truncated, yet these Ns are invisible to
    the AGP-gap machinery.  Exploding them into real AGP gap rows lets
    the whole span/abut/extension pipeline treat them like ordinary
    inter-component gaps, and the output AGP stays honest: filled bases
    appear as ``gapfill_tx`` components between the split pieces.

    An N run is exploded only when it is >= *min_len* bp and keeps >=
    *min_edge* bp of component sequence on BOTH sides (runs touching a
    component boundary are adjacent to an AGP gap anyway).  Returns
    (new_agp_lines, n_exploded); part numbers are renumbered per object.
    """
    from nearscaff.gapfill import _component_seq

    out: list = []
    n_exploded = 0
    part = 0
    cur = None
    for line in agp_lines:
        if line.object_name != cur:
            cur = line.object_name
            part = 0
        if not isinstance(line, AGPSeqLine) or \
                line.component_id not in contig_seqs:
            part += 1
            line.part_number = part
            out.append(line)
            continue
        seq = _component_seq(line, contig_seqs)
        runs = [(m.start(), m.end()) for m in _N_RUN.finditer(seq)
                if m.end() - m.start() >= min_len
                and m.start() >= min_edge
                and len(seq) - m.end() >= min_edge]
        if not runs:
            part += 1
            line.part_number = part
            out.append(line)
            continue
        # split into alternating pieces / gap rows (scaffold coords run
        # from line.object_beg; component coords honour orientation)
        cb, ce = line.component_beg, line.component_end

        def piece_coords(a: int, b: int) -> tuple:
            """oriented 0-based [a, b) -> AGP 1-based comp (beg, end)."""
            if line.orientation == "-":
                return ce - b + 1, ce - a
            return cb + a, cb + b - 1

        prev = 0
        for r0, r1 in runs:
            beg1, end1 = piece_coords(prev, r0)
            part += 1
            out.append(AGPSeqLine(
                line.object_name,
                line.object_beg + prev, line.object_beg + r0 - 1,
                part, line.component_type, line.component_id,
                beg1, end1, line.orientation))
            nlen = r1 - r0
            part += 1
            out.append(AGPGapLine(
                line.object_name,
                line.object_beg + r0, line.object_beg + r1 - 1,
                part, "N", nlen, "scaffold", "yes", "na"))
            n_exploded += 1
            prev = r1
        beg1, end1 = piece_coords(prev, len(seq))
        part += 1
        out.append(AGPSeqLine(
            line.object_name,
            line.object_beg + prev, line.object_end,
            part, line.component_type, line.component_id,
            beg1, end1, line.orientation))
    return out, n_exploded


def find_broken_loci(gff_path: str, trans_path: str, scaf_lens: dict,
                     pad: int = 2000) -> dict:
    """Loci whose miniprot translation contains X (CDS crosses an N run).

    *gff_path* / *trans_path* come from ``miniprot --gff`` and
    ``--trans`` on the same scaffold set.  Returns {scaffold: [(beg, end),
    ...]} (1-based, padded by *pad*, clamped, merged).
    """
    x_prots: set = set()
    seen: set = set()
    cur = None
    with open(trans_path) as fh:
        for line in fh:
            if line.startswith("##STA"):
                if cur and cur not in seen:
                    seen.add(cur)
                    if "X" in line:
                        x_prots.add(cur)
                cur = None
            elif not line.startswith("#"):
                cur = line.split("\t")[0]
    loci: dict[str, list] = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "mRNA":
                continue
            target = ""
            for f in p[8].split(";"):
                if f.startswith("Target="):
                    target = f[7:].split()[0]
            if target in x_prots and p[0] in scaf_lens:
                b = max(1, int(p[3]) - pad)
                e = min(scaf_lens[p[0]], int(p[4]) + pad)
                loci.setdefault(p[0], []).append((b, e))
    merged: dict[str, list] = {}
    for scaf, ivs in loci.items():
        ivs.sort()
        cur_ivs = []
        for b, e in ivs:
            if cur_ivs and b <= cur_ivs[-1][1]:
                cur_ivs[-1][1] = max(cur_ivs[-1][1], e)
            else:
                cur_ivs.append([b, e])
        merged[scaf] = [(b, e) for b, e in cur_ivs]
    return merged


def extract_read_pairs(paf_path: str, fq_files: list, out_files: list) -> int:
    """Extract reads (and mates) named in *paf_path* from paired FASTQ
    files.  Mates share the same base name (first token of the header).
    Returns the number of pairs written."""
    wanted: set = set()
    with open(paf_path) as f:
        for line in f:
            p = line.split("\t")
            if len(p) < 12 or p[2] == "*":
                continue
            wanted.add(p[0])
    import gzip as _gz
    for fq_in, fq_out in zip(fq_files, out_files):
        opener = _gz.open if fq_in.endswith(".gz") else open
        with opener(fq_in, "rt") as fi, open(fq_out, "w") as fo:
            while True:
                h = fi.readline()
                if not h:
                    break
                if not h.strip():
                    continue
                name = h[1:].split()[0]
                seq = fi.readline()
                plus = fi.readline()
                qual = fi.readline()
                if name in wanted:
                    fo.write(h + seq + plus + qual)
    return len(wanted)


def _ensure_query_miniprot(scaf_seqs: dict, proteins: str,
                           output_dir: str, miniprot: str = "miniprot",
                           threads: int = 4) -> tuple:
    """miniprot proteins -> scaffolds (--gff + --trans), cached in
    *output_dir*; returns (gff_path, trans_path)."""
    import subprocess
    scaf_fa = os.path.join(output_dir, "rnafill.scaffolds.fa")
    if not os.path.exists(scaf_fa):
        _write_records_fa([(n, s) for n, s in scaf_seqs.items()], scaf_fa)

    gff = os.path.join(output_dir, "rnafill.miniprot.gff")
    faa = os.path.join(output_dir, "rnafill.miniprot.trans.faa")
    for extra, out in ((["--gff"], gff), (["--trans"], faa)):
        if not os.path.exists(out):
            cmd = [miniprot, "-I", "-t", str(threads)] + extra + \
                  [scaf_fa, proteins]
            logger.info("Running: %s > %s", " ".join(cmd), out)
            with open(out, "w") as fh:
                r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE,
                                   text=True)
            if r.returncode != 0:
                raise RuntimeError(
                    f"miniprot failed: {r.stderr.strip()[-1000:]}")
    return gff, faa


def find_intact_proteins(trans_path: str,
                         min_cov: float = 0.8) -> set:
    """Proteins with at least one intact hit on the query: aligned
    coverage >= *min_cov* and translation free of X (no N crossing).

    Parses miniprot ``--trans`` output: PAF-like lines carry the protein
    length (field 2) and aligned protein span (fields 3-4); each is
    followed by a ``##STA`` line with the translated sequence.
    """
    intact: set = set()
    cur = None
    cov = 0.0
    with open(trans_path) as fh:
        for line in fh:
            if line.startswith("##STA"):
                if cur and cov >= min_cov and "X" not in line:
                    intact.add(cur)
                cur = None
            elif not line.startswith("#"):
                p = line.split("\t")
                cur = p[0]
                cov = (int(p[3]) - int(p[2])) / max(int(p[1]), 1)
    return intact


def run_recruit(agp_lines: list, contig_seqs: dict, scaf_seqs: dict,
                eligible: set, output_dir: str, reads: list,
                proteins: str, flank: int = TX_FLANK, bait_pad: int = 2000,
                threads: int = 4, minimap2: str = "minimap2",
                miniprot: str = "miniprot") -> dict:
    """Read-recruitment phase: combined bait of gap flanks + broken-gene
    loci, then extract the recruited read pairs for targeted assembly.

    Writes ``rnafill.recruit_1/2.fastq`` (or .fastq + .fastq for SE) plus
    the bait reference and intermediate miniprot outputs into
    *output_dir*.  Returns a stats dict.
    """
    import subprocess

    gff, faa = _ensure_query_miniprot(scaf_seqs, proteins, output_dir,
                                      miniprot=miniprot, threads=threads)

    scaf_lens = {n: len(s) for n, s in scaf_seqs.items()}
    broken = find_broken_loci(gff, faa, scaf_lens, pad=bait_pad)
    n_broken = sum(len(v) for v in broken.values())
    logger.info("Broken-gene bait loci (X-containing mRNAs +/- %d bp): %d",
                bait_pad, n_broken)

    # combined bait: gap flanks + broken loci
    records, _flank_len = build_flank_fasta(agp_lines, eligible, scaf_seqs,
                                            flank=flank)
    for scaf, ivs in sorted(broken.items()):
        s = scaf_seqs[scaf]
        for i, (b, e) in enumerate(ivs):
            records.append((f"{scaf}|broken{i}", s[b - 1:e]))
    bait_fa = os.path.join(output_dir, "rnafill.bait.fa")
    _write_records_fa(records, bait_fa)
    logger.info("Combined bait: %d flank + %d broken records -> %s",
                len(records) - n_broken, n_broken, bait_fa)

    bait_paf = os.path.join(output_dir, "rnafill.bait.paf")
    run_minimap2(bait_fa, list(reads), bait_paf, "splice", threads=threads,
                 extra=["-N", "20"], minimap2=minimap2)

    out1 = os.path.join(output_dir, "rnafill.recruit_1.fastq")
    out_files = [out1]
    if len(reads) > 1:
        out_files.append(os.path.join(output_dir, "rnafill.recruit_2.fastq"))
    n_names = extract_read_pairs(bait_paf, list(reads), out_files)
    logger.info("Recruited %d read names -> %s", n_names,
                " ".join(out_files))
    return {"bait_records": len(records), "broken_loci": n_broken,
            "recruited_names": n_names,
            "recruit_fastqs": out_files}


def _canonical_splice_edges(scaf_seq: str, gb: int, ge: int) -> bool:
    """True if the visible dinucleotides at the gap edges look like
    canonical splice sites — the left flank ends with the donor (GT) or
    the right flank starts with the acceptor (AG), on either strand —
    meaning the contig breaks are inside an intron and the gap content
    is (mostly) intronic.

    *gb*/*ge* are the 0-based half-open gap coordinates in the scaffold.
    Note the break must be exactly at the dinucleotide for the signal to
    be visible; breaks deeper inside an intron show nothing.  This is a
    one-sided heuristic: True means "definitely intronic", False means
    "no evidence either way".
    """
    left = scaf_seq[max(0, gb - 2):gb].upper()
    right = scaf_seq[ge:ge + 2].upper()
    if len(left) < 2 or len(right) < 2:
        return False
    return ((left.endswith("GT") and right.startswith("AG")) or
            (left.endswith("CT") and right.startswith("AC")))


def _gid_coords(gid: str) -> tuple:
    """gid "scaf:beg-end" (1-based) -> (scaf, gb, ge) 0-based half-open."""
    scaf, rng = gid.rsplit(":", 1)
    beg_s, end_s = rng.split("-")
    return scaf, int(beg_s) - 1, int(end_s)


def parse_tx_span_fills(paf_path: str, tx_files, scaf_seqs: dict,
                        fill_mode: str = "cdna",
                        min_tx: int = 1,
                        max_shortfall: int = TX_MAX_SHORTFALL,
                        abut_window: int = TX_ABUT_WINDOW,
                        max_depth_factor: float | None = 3.0,
                        detail: dict | None = None) -> tuple:
    """Extract fill sequences from transcripts hitting BOTH gap flanks.

    Same recruitment logic as gapfill.parse_span_fills (alignments must
    reach within *max_shortfall* of the gap edge; boundaries are
    extrapolated), but the middle segment is cDNA and is classified:

      middle <= abut_window  -> "abut" candidate (empty fill): the two
                                flanks are adjacent exons, the gap is
                                (mostly) intronic
      middle >  abut_window  -> fill candidate (the cDNA segment)

    In ``exon-only`` mode, fill candidates whose gap edges look like
    canonical splice sites (_canonical_splice_edges) are rejected as
    intronic; abut candidates are rejected as well.  Per gap the
    median-length candidate wins; the depth gate rejects over-collapsed
    gaps (multi-copy gene families piling up transcripts).

    When *detail* (a dict) is passed, each accepted closure is recorded
    as detail[gid] = (kind, transcript, mid_len) with kind "span"/"abut".

    Returns (fills, stats): fills = {gid: seq} ("" for abut), stats has
    counts per outcome.
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
            # prefer the primary alignment per side: paralogous flanks
            # recruit secondary hits that would pair L and R from
            # different loci into a cross-locus pseudo-fill
            primary = "tp:A:P" in p[12:]
            if side not in h or (primary and not h[side][5]):
                # store strand, qs, qe, ts, te, is_primary
                h[side] = (p[4], int(p[2]), int(p[3]), ts, te, primary)
    seqs = _fetch_read_seqs(tx_files, {qn for _g, qn in hits})

    cands: dict[str, list] = {}
    for (gid, qn), h in hits.items():
        if set(h) != {"L", "R"}:
            continue
        s = seqs.get(qn)
        if s is None:
            continue
        strand_l, qs_l, qe_l, _ts_l, te_l, _pl = h["L"]
        strand_r, qs_r, qe_r, ts_r, _te_r, _pr = h["R"]
        if strand_l != strand_r:
            continue
        len_l = flank_len[f"{gid}|L"]
        if strand_l == "+":
            # gap boundaries in transcript coords (extrapolated to edges)
            b0 = qe_l + (len_l - te_l)   # gap start
            b1 = qs_r - ts_r             # gap end
        else:
            b0 = qe_r + ts_r             # gap start (transcript coords)
            b1 = qs_l - (len_l - te_l)   # gap end
        mid = b1 - b0
        if mid <= 0:
            # flanks overlap / are adjacent in the transcript: intronic gap
            cands.setdefault(gid, []).append(("", max(mid, 0), qn))
            continue
        if b0 < 0 or b1 > len(s):
            continue
        seg = s[b0:b1]
        if strand_l == "-":
            seg = _revcomp(seg)
        if "N" not in seg.upper():
            cands.setdefault(gid, []).append((seg, mid, qn))

    # depth gate
    skip: set = set()
    if max_depth_factor and cands:
        depths = sorted(len(v) for v in cands.values())
        baseline = depths[(len(depths) - 1) // 2] or 1
        cutoff = max_depth_factor * max(baseline, min_tx)
        for gid, lst in cands.items():
            if len(lst) > cutoff:
                skip.add(gid)
        if skip:
            logger.info("Span depth gate: baseline %d transcripts, "
                        "rejecting %d over-collapsed gaps (>%.1fx)",
                        baseline, len(skip), max_depth_factor)

    fills: dict[str, str] = {}
    stats = {"span_filled": 0, "abut": 0, "intron_skipped": 0,
             "depth_rejected": len(skip), "multi_candidate_gaps": 0}
    for gid, lst in cands.items():
        if gid in skip or len(lst) < min_tx:
            continue
        if len(lst) > 1:
            stats["multi_candidate_gaps"] += 1
        lst.sort(key=lambda x: x[1])
        seg, mid, qn = lst[len(lst) // 2]
        if mid <= abut_window:
            if fill_mode == "exon-only":
                stats["intron_skipped"] += 1
                continue
            fills[gid] = ""
            stats["abut"] += 1
            if detail is not None:
                detail[gid] = ("abut", qn, mid)
            continue
        if fill_mode == "exon-only":
            scaf, gb, ge = _gid_coords(gid)
            if _canonical_splice_edges(scaf_seqs[scaf], gb, ge):
                stats["intron_skipped"] += 1
                continue
        fills[gid] = seg
        stats["span_filled"] += 1
        if detail is not None:
            detail[gid] = ("span", qn, mid)
    return fills, stats


def parse_tx_extensions(paf_path: str, tx_files,
                        min_tail: int = TX_MIN_EXT,
                        max_tail: int = TX_MAX_EXT,
                        min_clip: int = 300, edge: int = 10,
                        max_tails: int = 30) -> dict:
    """One-sided extensions: transcripts that cover a gap edge (primary
    alignment reaching the edge) and extend into the gap, with NO
    requirement of evidence from the other side.  The gap stays open;
    the covered part is written next to its flank (see write_outputs
    *extensions*).

    Tail orientations follow parse_tails: L tails run flank -> deep gap
    (written adjacent to the left flank), R tails run deep gap -> flank
    (written adjacent to the right flank).  Tails longer than *max_tail*
    are truncated at the deep end.  Per side the median-length candidate
    wins (isoforms).

    Returns {gid: (left_seq_or_None, right_seq_or_None, source_tx)}.
    """
    flank_len: dict[str, int] = {}
    hits: dict[tuple, tuple] = {}
    with open(paf_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[2] == "*" or p[5] == "*":
                continue
            if "tp:A:P" not in p[12:]:      # primary only: one-sided
                continue                    # evidence is weak by nature
            gid, side = p[5].rsplit("|", 1)
            flank_len[p[5]] = int(p[6])
            ql, qs, qe = int(p[1]), int(p[2]), int(p[3])
            strand, ts, te = p[4], int(p[7]), int(p[8])
            L = flank_len[p[5]]
            if side == "L" and abs(te - L) <= edge:
                if strand == "+" and ql - qe >= min_clip:
                    hits[(gid, "L", p[0])] = (qe, strand)
                elif strand == "-" and qs >= min_clip:
                    hits[(gid, "L", p[0])] = (qs, strand)
            elif side == "R" and ts <= edge:
                if strand == "+" and qs >= min_clip:
                    hits[(gid, "R", p[0])] = (qs, strand)
                elif strand == "-" and ql - qe >= min_clip:
                    hits[(gid, "R", p[0])] = (qe, strand)

    seqs = _fetch_read_seqs(tx_files, {qn for _g, _s, qn in hits})

    tails: dict[str, dict] = {}
    for (gid, side, qn), (pos, strand) in hits.items():
        s = seqs.get(qn)
        if s is None:
            continue
        if side == "L":
            tail = s[pos:] if strand == "+" else _revcomp(s[:pos])
            tail = tail[:max_tail]          # truncate the deep end
        else:
            tail = s[:pos] if strand == "+" else _revcomp(s[pos:])
            tail = tail[len(tail) - max_tail:]
        if len(tail) < min_tail or "N" in tail.upper():
            continue
        t = tails.setdefault(gid, {"L": [], "R": []})
        if len(t[side]) < max_tails:
            t[side].append((tail, qn))

    out: dict = {}
    for gid, t in tails.items():
        left = right = None
        qn = ""
        if t["L"]:
            t["L"].sort(key=lambda x: len(x[0]))
            left, qn = t["L"][len(t["L"]) // 2]
        if t["R"]:
            t["R"].sort(key=lambda x: len(x[0]))
            right, qnr = t["R"][len(t["R"]) // 2]
            qn = qn or qnr
        if left or right:
            out[gid] = (left, right, qn)
    return out


# ---------------------------------------------------------------------
# ref-guided whole-gene placement
#
# When a gene's entire locus is N in the query there is no flank for a
# transcript to anchor to, yet a close reference knows where the gene
# sits and nearscaff's tiered PAF already maps query contigs onto the
# reference.  The functions below bracket each unfilled gap on the
# reference (via the flanking contigs' alignments), find the reference
# genes lying fully inside the bracket (those genes are necessarily
# missing from the query — the region is all N), assign transcripts to
# them by mapping the transcripts against the reference, and write the
# transcript sequence into the gap as exon components separated by
# estimated-N intron/intergenic spacers.
# ---------------------------------------------------------------------

def read_ref_genes(gff_path: str) -> dict:
    """Ref gene loci from a GFF3: {chr: [(beg0, end, gene_id, strand)]}.

    mRNA/transcript lines are preferred; if none are present the file is
    re-scanned for gene lines.  Coordinates are 0-based half-open.
    """
    def _scan(kinds: set) -> dict:
        genes: dict[str, list] = {}
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] not in kinds:
                    continue
                gid = ""
                for f in p[8].split(";"):
                    if f.startswith("ID="):
                        gid = f[3:]
                        break
                genes.setdefault(p[0], []).append(
                    (int(p[3]) - 1, int(p[4]), gid or f"{p[0]}:{p[3]}-{p[4]}",
                     p[6]))
        return genes

    genes = _scan({"mRNA", "transcript"})
    if not genes:
        genes = _scan({"gene"})
    for v in genes.values():
        v.sort()
    return genes


def derive_ref_genes(ref_fa: str, proteins: str, work_dir: str,
                     miniprot: str = "miniprot", threads: int = 4) -> dict:
    """Ref gene loci by mapping reference proteins back onto the
    reference genome with miniprot (cached as rnafill.refgenes.gff)."""
    import subprocess
    gff = os.path.join(work_dir, "rnafill.refgenes.gff")
    if not os.path.exists(gff):
        cmd = [miniprot, "-I", "-t", str(threads), "--gff", ref_fa, proteins]
        logger.info("Running: %s > %s", " ".join(cmd), gff)
        with open(gff, "w") as fh:
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE,
                               text=True)
        if r.returncode != 0:
            raise RuntimeError(f"miniprot failed: {r.stderr.strip()[-1000:]}")
    genes: dict[str, list] = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "mRNA":
                continue
            target = ""
            for f in p[8].split(";"):
                if f.startswith("Target="):
                    target = f[7:].split()[0]
            genes.setdefault(p[0], []).append(
                (int(p[3]) - 1, int(p[4]), target or f"{p[0]}:{p[3]}-{p[4]}",
                 p[6]))
    for v in genes.values():
        v.sort()
    return genes


def contig_ref_index(tiered_paf: str) -> dict:
    """{contig: [(qs, qe, strand, ref_chr, ts, te, is_primary)]} from
    nearscaff_tiered.paf (all alignments, not just the tier tag)."""
    idx: dict[str, list] = {}
    with open(tiered_paf) as f:
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[5] == "*":
                continue
            idx.setdefault(p[0], []).append(
                (int(p[2]), int(p[3]), p[4], p[5],
                 int(p[7]), int(p[8]), "tp:A:P" in p[12:]))
    return idx


def _map_to_ref(aligns: list, coord: int) -> tuple | None:
    """Map a 1-based contig coordinate to (ref_chr, ref_pos0, strand)
    via the covering alignment (primary preferred, then longest)."""
    c0 = coord - 1
    best = None
    for qs, qe, strand, tname, ts, te, primary in aligns:
        if not (qs <= c0 < qe):
            continue
        ref = ts + (c0 - qs) if strand == "+" else te - (c0 - qs) - 1
        cand = (tname, ref, strand, primary, qe - qs)
        if best is None or (primary, qe - qs) > (best[3], best[4]):
            best = cand
    return best[:3] if best else None


def gap_ref_brackets(agp_lines: list, gaps: set, contig_idx: dict,
                     max_bracket: int = PLACE_MAX_BRACKET) -> dict:
    """Bracket each gap on the reference: {agp_idx: (ref_chr, r0, r1, flip)}.

    The gap's flanking components' gap-side ends are mapped through the
    tiered PAF.  A bracket is kept only when both ends land on the same
    reference chromosome with consistent strand and span at most
    *max_bracket*.  *flip* is True when scaffold-forward runs against
    the reference forward direction.
    """
    brackets: dict = {}
    for i, line in enumerate(agp_lines):
        if not isinstance(line, AGPGapLine) or i not in gaps:
            continue
        left = right = None
        for j in range(i - 1, -1, -1):
            if agp_lines[j].object_name != line.object_name:
                break
            if isinstance(agp_lines[j], AGPSeqLine):
                left = agp_lines[j]
                break
        for j in range(i + 1, len(agp_lines)):
            if agp_lines[j].object_name != line.object_name:
                break
            if isinstance(agp_lines[j], AGPSeqLine):
                right = agp_lines[j]
                break
        if left is None or right is None:
            continue
        # gap-side end of each flank, in 1-based contig coordinates
        lc = left.component_end if left.orientation == "+" \
            else left.component_beg
        rc = right.component_beg if right.orientation == "+" \
            else right.component_end
        lm = _map_to_ref(contig_idx.get(left.component_id, []), lc)
        rm = _map_to_ref(contig_idx.get(right.component_id, []), rc)
        if lm is None or rm is None:
            continue
        if lm[0] != rm[0] or lm[2] != rm[2]:
            continue
        r0, r1 = sorted((lm[1], rm[1]))
        if r1 - r0 <= 0 or r1 - r0 > max_bracket:
            continue
        brackets[i] = (lm[0], r0, r1, lm[2] == "-")
    return brackets


_CG_OP = re.compile(r"(\d+)([MIDNSHP=X])")


def assign_tx_to_genes(refmap_paf: str, ref_genes: dict,
                       min_ovlp: float = PLACE_MIN_OVLP) -> dict:
    """Assign transcripts to ref genes via their splice alignments.

    Only primary alignments are used; an alignment is assigned to the
    same-chromosome gene it overlaps best, provided
    max(ovlp/gene_len, ovlp/aln_span) >= *min_ovlp*.  Each transcript is
    assigned to at most one gene (its best), each gene keeps its most
    overlapping transcript.  Returns {gene_id: (tx_name, paf_fields)}.
    """
    per_tx: dict[str, tuple] = {}    # tx -> (ovlp, gene_id, fields)
    per_gene: dict[str, tuple] = {}  # gene_id -> (ovlp, tx, fields)
    with open(refmap_paf) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[2] == "*" or p[5] == "*":
                continue
            if "tp:A:P" not in p[12:]:
                continue
            ts, te = int(p[7]), int(p[8])
            best = None
            for gb, ge, gid, _st in ref_genes.get(p[5], []):
                if gb >= te:
                    break
                ov = min(te, ge) - max(ts, gb)
                if ov <= 0:
                    continue
                if max(ov / (ge - gb), ov / max(te - ts, 1)) < min_ovlp:
                    continue
                if best is None or ov > best[0]:
                    best = (ov, gid)
            if best is None:
                continue
            ov, gid = best
            if ov > per_tx.get(p[0], (0,))[0]:
                per_tx[p[0]] = (ov, gid, p)
    for tx, (ov, gid, p) in per_tx.items():
        if ov > per_gene.get(gid, (0,))[0]:
            per_gene[gid] = (ov, tx, p)
    return {gid: (tx, p) for gid, (_ov, tx, p) in per_gene.items()}


def build_gene_placement(tx_seq: str, paf_fields: list, r0: int, r1: int,
                         min_spacer: int = PLACE_MIN_SPACER,
                         max_spacer: int = PLACE_MAX_SPACER) -> list:
    """Exon/N-spacer blocks for one assigned transcript, REF-forward.

    The transcript's cg CIGAR (spliced alignment vs the reference) is
    walked with a query and a reference cursor: M/=/X chunks are exonic
    sequence cut from the transcript, N ops become intron placeholders
    of the reference-derived length (clamped to [min_spacer, max_spacer])
    — a pure (None, est_len) gap-row tuple when the intron is clipped by
    the bracket, otherwise a sequence block with the canonical splice
    consensus stamped on its edges (GT..AG, or CT..AC for minus-strand
    genes) so downstream spliced aligners can cross it; I bases are kept
    (appended to the current exon), D ops are skipped.  Exons are
    clipped to the bracket [r0, r1) so the placement never duplicates
    flanking sequence.  Returns a list mixing str blocks and
    (None, est_len) tuples; empty when no exonic base falls inside the
    bracket.

    CIGAR cursor convention (verified against minimap2 2.x PAF output):
    the reference cursor always ascends from ts; the query cursor
    ascends from qs for "+" alignments but DESCENDS from qe for "-"
    alignments (the CIGAR is written against the reverse-complemented
    query), so "-" chunks are revcomped individually and the CIGAR
    order already is the ref-forward order — no block reversal needed.
    """
    strand = paf_fields[4]
    qs0, qe0 = int(paf_fields[2]), int(paf_fields[3])
    cg = None
    for t in paf_fields[12:]:
        if t.startswith("cg:Z:"):
            cg = t[5:]
            break
    blocks: list = []
    if cg is None:                      # no CIGAR: one exon block
        seg = tx_seq[qs0:qe0]
        if strand == "-":
            seg = _revcomp(seg)
        blocks = [seg] if seg else []
    else:
        cur: list[str] = []
        qpos = qs0 if strand == "+" else qe0
        rpos = int(paf_fields[7])

        def flush() -> None:
            if cur:
                blocks.append("".join(cur))
                cur.clear()

        for n_s, op in _CG_OP.findall(cg):
            n = int(n_s)
            if op in "M=X":
                lo, hi = max(rpos, r0), min(rpos + n, r1)
                if lo < hi:
                    off = lo - rpos
                    ln = hi - lo
                    if strand == "+":
                        cur.append(tx_seq[qpos + off:qpos + off + ln])
                    else:
                        # block covers query [qpos-n, qpos); ref offset
                        # counts back from the block's query end
                        qb_hi = qpos - off
                        cur.append(_revcomp(tx_seq[qb_hi - ln:qb_hi]))
                qpos += n if strand == "+" else -n
                rpos += n
            elif op == "I":
                if r0 <= rpos <= r1:
                    if strand == "+":
                        cur.append(tx_seq[qpos:qpos + n])
                    else:
                        cur.append(_revcomp(tx_seq[qpos - n:qpos]))
                qpos += n if strand == "+" else -n
            elif op == "S":
                qpos += n if strand == "+" else -n
            elif op == "N":
                lo, hi = max(rpos, r0), min(rpos + n, r1)
                if hi - lo >= min_spacer:
                    flush()
                    sp = min(hi - lo, max_spacer)
                    if lo == rpos and hi == rpos + n and sp >= 4:
                        # unclipped intron: splicing is proven by the
                        # CIGAR — stamp the canonical splice consensus
                        # (strand-aware) onto the placeholder so miniprot
                        # / ab initio annotators can splice across it
                        donor, acc = ("GT", "AG") if strand == "+" \
                            else ("CT", "AC")
                        blocks.append(donor + "N" * (sp - 4) + acc)
                    else:
                        blocks.append((None, sp))
                rpos += n
            elif op == "D":
                rpos += n
            # H/P consume nothing on either cursor here
        flush()
    blocks = [b for b in blocks
              if (isinstance(b, str) and b) or isinstance(b, tuple)]
    return blocks


def plan_placements(brackets: dict, ref_genes: dict, gene_tx: dict,
                    tx_seqs: dict, max_genes: int = PLACE_MAX_GENES,
                    max_spacer: int = PLACE_MAX_SPACER,
                    min_edge: int = PLACE_MIN_EDGE) -> tuple:
    """Assemble per-gap placements in scaffold-forward order.

    *brackets* comes from gap_ref_brackets (restricted to the gaps still
    unfilled), *gene_tx* from assign_tx_to_genes, *tx_seqs* holds the
    assigned transcripts' sequences.  Genes fully inside the bracket are
    laid out in reference order; intergenic and bracket-edge spacers are
    estimated-N gap rows (clamped; edge spacers shorter than *min_edge*
    are dropped so a gene flush against the flank hugs it).  Returns
    (placements, detail): placements = {agp_idx: [blocks]},
    detail = {agp_idx: [(gene_id, tx_name, exon_bases)]}.
    """
    placements: dict = {}
    detail: dict = {}
    for gidx, (chr_, r0, r1, flip) in sorted(brackets.items()):
        genes = [g for g in ref_genes.get(chr_, [])
                 if g[0] >= r0 and g[1] <= r1 and g[2] in gene_tx]
        if not genes:
            continue
        if len(genes) > max_genes:
            logger.info("Placement skipped for gap idx %d (%s:%d-%d): "
                        "%d genes > --place-max-genes %d",
                        gidx, chr_, r0, r1, len(genes), max_genes)
            continue
        genes.sort()
        placed: list = []     # (gb, ge, gene_id, tx_name, blocks)
        for gb, ge, gid, _st in genes:
            txn, paf = gene_tx[gid]
            seq = tx_seqs.get(txn)
            if seq is None:
                continue
            gblocks = build_gene_placement(seq, paf, r0, r1,
                                           max_spacer=max_spacer)
            if gblocks:
                placed.append((gb, ge, gid, txn, gblocks))
        if not placed:
            continue
        blocks: list = []
        d = placed[0][0] - r0
        if d >= min_edge:
            blocks.append((None, min(d, max_spacer)))
        for k, (gb, ge, gid, txn, gblocks) in enumerate(placed):
            blocks.extend(gblocks)
            if k + 1 < len(placed):
                d = placed[k + 1][0] - ge
                blocks.append((None, min(max(d, min_edge), max_spacer)))
        d = r1 - placed[-1][1]
        if d >= min_edge:
            blocks.append((None, min(d, max_spacer)))
        if flip:
            blocks = [_revcomp(b) if isinstance(b, str) else b
                      for b in reversed(blocks)]
        placements[gidx] = blocks
        detail[gidx] = [(gid, txn,
                         sum(len(b) for b in gblocks if isinstance(b, str)))
                        for _gb, _ge, gid, txn, gblocks in placed]
    return placements, detail


def run_rnafill(agp_path: str, query_fasta: str, tiered_paf: str,
                output_dir: str, transcripts: list | None = None,
                fill_mode: str = "cdna",
                tx_preset: str = "splice:hq",
                threads: int = 4,
                fill_tiers: set = frozenset({"protein", "asm5"}),
                min_span_tx: int = 1,
                overlap_closure: bool = True,
                one_sided: bool = False,
                ext_min_tail: int = TX_MIN_EXT,
                ext_max_tail: int = TX_MAX_EXT,
                abut_window: int = TX_ABUT_WINDOW,
                max_depth_factor: float | None = 3.0,
                flank: int = TX_FLANK,
                internal_n: int = 0,
                reads: list | None = None,
                proteins: str | None = None,
                bait_pad: int = 2000,
                miniprot: str = "miniprot",
                minimap2: str = "minimap2",
                recruit_extra: list | None = None,
                ref: str | None = None,
                ref_gff: str | None = None,
                place_max_bracket: int = PLACE_MAX_BRACKET,
                place_max_genes: int = PLACE_MAX_GENES,
                place_min_ovlp: float = PLACE_MIN_OVLP,
                place_intact_cov: float = PLACE_INTACT_COV,
                place_max_spacer: int = PLACE_MAX_SPACER) -> dict:
    """Transcript-guided filling of genic gaps via flank recruitment.

    *transcripts* is a list of FASTA/FASTQ(.gz) files (Iso-Seq, ONT cDNA
    or assembled transcripts; assemble short-read RNA-seq first).
    See module docstring for the fill modes.  With *one_sided*, gaps
    that no transcript spans are still partially filled when a
    transcript covers one edge and extends >= *ext_min_tail* bp into
    the gap (the gap stays open; the known sequence is written next to
    its flank).  With *internal_n* > 0, N runs of at least that length
    inside components are exploded into gap rows first and filled by the
    same machinery (short-read assemblies carry N runs inside their
    scaffolds; HiFi assemblies have none).  With *ref* (reference genome
    FASTA), gaps still unfilled after all flank-based mechanisms get a
    ref-guided placement pass: the gap is bracketed on the reference via
    its flanking contigs' tiered-PAF alignments, reference genes lying
    fully inside the bracket (necessarily missing from the query — the
    region is all N) are matched to transcripts mapped against the
    reference, and the transcript sequence is written into the gap as
    exon components separated by estimated-N intron/intergenic spacers.
    Reference gene coordinates come from *ref_gff* when given, else from
    miniprot mapping *proteins* back onto *ref*.  With miniprot-derived
    loci, an anti-duplication filter applies: a ref gene whose protein
    already has an intact hit anywhere in the query (aligned coverage >=
    *place_intact_cov*, translation free of X) is not missing and is
    excluded from placement.  Returns the report dict.
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(agp_path) as f:
        agp_lines = AGPReader().parse(f.read())
    contig_seqs = read_fasta(query_fasta)
    if internal_n:
        agp_lines, n_int = explode_internal_ns(agp_lines, contig_seqs,
                                               min_len=internal_n)
        logger.info("Internal N runs exploded into fillable gap rows: %d "
                    "(min length %d bp)", n_int, internal_n)
    tiers = read_tiers(tiered_paf)
    eligible = classify_gaps(agp_lines, tiers, allowed=fill_tiers)
    n_gaps = sum(1 for l in agp_lines if isinstance(l, AGPGapLine))
    logger.info("%d gaps total, %d eligible for filling (flanks in %s)",
                n_gaps, len(eligible), "/".join(sorted(fill_tiers)))

    scaf_seqs = scaffold_sequences(agp_lines, contig_seqs)

    report: dict = {"gaps_total": n_gaps, "gaps_eligible": len(eligible)}
    if reads:
        if not proteins:
            raise ValueError("read recruitment requires reference proteins "
                             "(--proteins)")
        report["recruit"] = run_recruit(
            agp_lines, contig_seqs, scaf_seqs, eligible, output_dir,
            reads, proteins, flank=flank, bait_pad=bait_pad,
            threads=threads, minimap2=minimap2, miniprot=miniprot)
        if not transcripts:
            logger.info("Recruitment done. Assemble the recruited reads "
                        "with any transcript assembler (e.g. Trinity), "
                        "then re-run rna-fill with -T to fill.")
            return report
    if not transcripts:
        raise ValueError("rna-fill requires transcripts (-T) or reads "
                         "(--reads)")

    records, flank_len = build_flank_fasta(agp_lines, eligible, scaf_seqs,
                                           flank=flank)
    flanks_fa = os.path.join(output_dir, "rnafill.flanks.fa")
    _write_records_fa(records, flanks_fa)
    logger.info("Flank mini-reference: %d records -> %s",
                len(records), flanks_fa)

    paf = os.path.join(output_dir, "rnafill.flanks.paf")
    run_minimap2(flanks_fa, transcripts, paf, tx_preset, threads=threads,
                 extra=recruit_extra or ["-N", "20"], minimap2=minimap2)

    detail: dict = {}
    fills, stats = parse_tx_span_fills(
        paf, transcripts, scaf_seqs, fill_mode=fill_mode, min_tx=min_span_tx,
        abut_window=abut_window, max_depth_factor=max_depth_factor,
        detail=detail)
    logger.info("Transcript span closures: %d fills, %d abut, %d skipped "
                "as intronic (%d gaps had >1 candidate)",
                stats["span_filled"], stats["abut"],
                stats["intron_skipped"], stats["multi_candidate_gaps"])

    ovl: dict = {}
    if overlap_closure:
        tails = parse_tails(paf, flank_len, transcripts)
        prefix = os.path.join(output_dir, "rnafill")
        ovl = overlap_closure_fills(tails, scaf_seqs, prefix,
                                    threads=threads, min_ovlp=TX_MIN_OVLP,
                                    min_ident=TX_MIN_OVLP_IDENT,
                                    uniq_preset="splice",
                                    minimap2=minimap2)
        if fill_mode == "exon-only" and ovl:
            ovl = {g: s for g, s in ovl.items()
                   if not _canonical_splice_edges(
                       scaf_seqs[_gid_coords(g)[0]],
                       _gid_coords(g)[1], _gid_coords(g)[2])}
        logger.info("Tail-overlap closures (splice uniqueness): %d gaps",
                    len(ovl))

    merged = dict(ovl)
    merged.update(fills)   # span fills win on conflict
    for g in ovl:
        if g not in detail:
            detail[g] = ("overlap", "", len(ovl[g]))
    gid_to_idx = {f"{l.object_name}:{l.object_beg}-{l.object_end}": i
                  for i, l in enumerate(agp_lines)
                  if isinstance(l, AGPGapLine)}
    fills_idx = {gid_to_idx[g]: s for g, s in merged.items()
                 if g in gid_to_idx}

    # one-sided extensions for gaps no transcript spans
    ext_idx: dict = {}
    if one_sided:
        ext = parse_tx_extensions(paf, transcripts,
                                  min_tail=ext_min_tail,
                                  max_tail=ext_max_tail)
        for g, (lseq, rseq, qn) in ext.items():
            if g not in gid_to_idx or g in merged:
                continue
            ext_idx[gid_to_idx[g]] = (lseq, rseq)
            kind = "ext_LR" if lseq and rseq else \
                   "ext_L" if lseq else "ext_R"
            detail[g] = (kind, qn,
                         (len(lseq) if lseq else 0) +
                         (len(rseq) if rseq else 0))
        logger.info("One-sided extensions: %d gaps (%d bp written next to "
                    "flanks)", len(ext_idx),
                    sum(len(l or "") + len(r or "")
                        for l, r in ext_idx.values()))

    # ref-guided whole-gene placement for gaps nothing else could fill
    place_idx: dict = {}
    placed_rows: list = []
    if ref:
        remaining = set(eligible) - set(fills_idx) - set(ext_idx)
        cidx = contig_ref_index(tiered_paf)
        brackets = gap_ref_brackets(agp_lines, remaining, cidx,
                                    max_bracket=place_max_bracket)
        logger.info("Ref-guided placement: %d/%d unfilled gaps bracketed "
                    "on the reference", len(brackets), len(remaining))
        if brackets:
            if ref_gff:
                ref_genes = read_ref_genes(ref_gff)
            else:
                if not proteins:
                    raise ValueError(
                        "ref-guided placement without --ref-gff requires "
                        "--proteins (miniprot derives the ref gene loci)")
                ref_genes = derive_ref_genes(ref, proteins, output_dir,
                                             miniprot=miniprot,
                                             threads=threads)
            n_genes = sum(len(v) for v in ref_genes.values())
            logger.info("Reference gene loci: %d (%s)", n_genes,
                        ref_gff or "miniprot-derived")
            if ref_gff:
                logger.info("NOTE: --ref-gff loci skip the intact-protein "
                            "anti-duplication check (no protein-id mapping)")
            else:
                # anti-duplication: a ref gene whose protein already has an
                # intact hit anywhere in the query is not missing — do not
                # place a second copy (misplaced contigs / paralogs would
                # otherwise duplicate it)
                _gff_q, faa_q = _ensure_query_miniprot(
                    scaf_seqs, proteins, output_dir, miniprot=miniprot,
                    threads=threads)
                intact = find_intact_proteins(
                    faa_q, min_cov=place_intact_cov)
                n_before = sum(len(v) for v in ref_genes.values())
                ref_genes = {c: [g for g in gs if g[2] not in intact]
                             for c, gs in ref_genes.items()}
                logger.info("Intact-protein filter: %d/%d ref loci already "
                            "present in the query (cov>=%.2f, no X) — "
                            "excluded from placement",
                            n_before - sum(len(v) for v in
                                           ref_genes.values()),
                            n_before, place_intact_cov)
            refmap = os.path.join(output_dir, "rnafill.refmap.paf")
            if not os.path.exists(refmap):
                run_minimap2(ref, transcripts, refmap, tx_preset,
                             threads=threads, extra=["-N", "5", "-c"],
                             minimap2=minimap2)
            gene_tx = assign_tx_to_genes(refmap, ref_genes,
                                         min_ovlp=place_min_ovlp)
            logger.info("Transcripts assigned to ref genes: %d",
                        len(gene_tx))
            tx_seqs = _fetch_read_seqs(
                transcripts, {txn for txn, _p in gene_tx.values()})
            place_idx, pdetail = plan_placements(
                brackets, ref_genes, gene_tx, tx_seqs,
                max_genes=place_max_genes, max_spacer=place_max_spacer)
            idx_to_gid = {i: g for g, i in gid_to_idx.items()}
            for gidx, metas in pdetail.items():
                for gene_id, txn, nb in metas:
                    placed_rows.append((idx_to_gid[gidx], gene_id, txn, nb))
            logger.info("Ref-guided placements: %d gaps, %d genes (%d bp "
                        "of exonic sequence)", len(place_idx),
                        len(placed_rows),
                        sum(nb for _g, _i, _t, nb in placed_rows))

    out_agp = os.path.join(output_dir, "nearscaff.rnafill.agp")
    out_fasta = os.path.join(output_dir, "nearscaff.rnafill.scaffolds.fa")
    fill_report = write_outputs(agp_lines, eligible, fills_idx, contig_seqs,
                                out_agp, out_fasta, fill_prefix="gapfill_tx",
                                extensions=ext_idx, placements=place_idx)
    report.update(fill_report)
    report.update(stats)
    report["overlap_closed"] = len(ovl)
    report["genes_placed"] = len(placed_rows)
    report["fill_mode"] = fill_mode

    # per-closure manifest for downstream cross-validation
    closures_tsv = os.path.join(output_dir, "rnafill.closures.tsv")
    with open(closures_tsv, "w") as f:
        f.write("gid\tscaffold\tgap_beg\tgap_end\tkind\tfill_len"
                "\tsource_tx\n")
        for g in sorted(detail):
            kind, qn, mid = detail[g]
            scaf, gb, ge = _gid_coords(g)
            f.write(f"{g}\t{scaf}\t{gb + 1}\t{ge}\t{kind}\t{mid}\t{qn}\n")
        for g, gene_id, txn, nb in placed_rows:
            scaf, gb, ge = _gid_coords(g)
            f.write(f"{g}\t{scaf}\t{gb + 1}\t{ge}\tplaced:{gene_id}\t{nb}"
                    f"\t{txn}\n")
    logger.info("Closure manifest: %s", closures_tsv)
    logger.info("rna-fill complete (%s mode): %d/%d eligible gaps closed "
                "(%d bp filled)", fill_mode, report["gaps_closed"],
                report["gaps_eligible"], report["bases_filled"])
    logger.info("Gap-filled AGP: %s", out_agp)
    logger.info("Gap-filled FASTA: %s", out_fasta)
    if fill_mode == "cdna":
        logger.info("NOTE: fills are cDNA — introns inside closed gaps "
                    "are compressed out; use for annotation completeness")
    return report
