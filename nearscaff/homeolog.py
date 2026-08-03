"""Allo4D-style homeolog block-pair discovery (robust reimplementation).

Finds high-quality homeologous block pairs between a (mixed) polyploid
assembly and a diploid reference, without the upstream Allo4D scripts'
failure modes (parallel-write races, MCScanX row-orientation and
pandas-version bugs — see docs/roadmap.md).

Chain:
  1. ``annotate_genome`` — miniprot --gff + gffread, with CDS stop-trim
     (gffread includes the stop codon; downstream Ks/translation checks
     require len(cds) == 3 * len(pep)) and ID prefixing (miniprot names
     every run's genes MP000001..., which collides across genomes).
  2. ``find_homeolog_pairs`` — jcvi ``compara.catalog ortholog`` (LAST +
     built-in synteny scan); sp2 genes with exactly 2 distinct sp4
     anchor hits form 1v2 clusters.  Row orientation is normalized by
     BED membership, not by column position or ID prefix.
  3. ``cluster_side_blocks`` — per cluster, the sp41-side and sp42-side
     blocks as (scaffold, start, end) spans; sides spanning multiple
     scaffolds are dropped.

External tools (must be on PATH or passed explicitly): miniprot, gffread,
jcvi (with LAST; command prefix overridable via NEARSCAFF_JCVI,
e.g. ``/path/to/python -m jcvi.compara.catalog``).
"""

import logging
import os
import subprocess
from collections import defaultdict

logger = logging.getLogger("nearscaff.homeolog")


def _tool(name: str) -> str:
    """External tool path, overridable via NEARSCAFF_<NAME> env var."""
    return os.environ.get(f"NEARSCAFF_{name.upper()}", name)


def _run(cmd: list[str], what: str):
    logger.info("  %s ...", what)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{what} failed: {result.stderr.strip()[:500]}")
    return result


# ============================================================================
# 1. Annotation (miniprot --gff + gffread)
# ============================================================================

def gene_level_pep(ref_pep: str, out_path: str) -> str:
    """Keep the longest isoform per gene (strips .mRNAn / .tN suffixes).

    Transcript-level protein files (e.g. one entry per mRNA) inflate the
    collinearity hit counts downstream; Allo4D needs gene-level inputs.
    If fewer than 10% of IDs carry an isoform suffix the file is assumed
    gene-level already and is returned unmodified.
    """
    import re

    def read_fasta(path):
        seqs = {}
        name, parts = None, []
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(parts)
                    name = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line.strip())
        if name:
            seqs[name] = "".join(parts)
        return seqs

    seqs = read_fasta(ref_pep)
    iso = [n for n in seqs if re.search(r'\.(mRNA|t)\d+$', n)]
    if len(iso) < 0.1 * len(seqs):
        logger.info("  %s: %d sequences, no isoform suffixes — using as-is",
                    os.path.basename(ref_pep), len(seqs))
        return ref_pep
    best: dict[str, tuple[str, str]] = {}
    for n, s in seqs.items():
        gene = re.sub(r'\.(mRNA|t)\d+$', '', n)
        if gene not in best or len(s) > len(best[gene][1]):
            best[gene] = (n, s)
    with open(out_path, "w") as fh:
        for _gene, (n, s) in sorted(best.items()):
            fh.write(f">{n}\n{s}\n")
    logger.info("  %s: %d sequences -> %d genes (longest isoform)",
                os.path.basename(ref_pep), len(seqs), len(best))
    return out_path


def annotate_genome(genome_fa: str, ref_pep: str, work_dir: str, label: str,
                    threads: int = 4, prefix: str = "MP",
                    miniprot: str | None = None, gffread: str | None = None,
                    ) -> tuple[str, str, str, str]:
    """Annotate *genome_fa* with *ref_pep*; returns (gff3, pep, cds, bed).

    mRNA IDs get *prefix* (passed to miniprot -P) so annotations of
    different genomes never collide.  CDS is stop-trimmed and pep/bed are
    filtered to the IDs that survived.
    """
    miniprot = miniprot or _tool("miniprot")
    gffread = gffread or _tool("gffread")
    os.makedirs(work_dir, exist_ok=True)
    gff3 = os.path.join(work_dir, f"{label}.gff3")
    pep = os.path.join(work_dir, f"{label}.pep")
    cds = os.path.join(work_dir, f"{label}.cds")
    bed = os.path.join(work_dir, f"{label}.bed")
    if all(os.path.exists(p) and os.path.getsize(p) > 0
           for p in (pep, cds, bed)):
        logger.info("  %s annotation cached — skipping", label)
        return gff3, pep, cds, bed

    index = os.path.join(work_dir, f"{label}.mpi")
    if not os.path.exists(index):
        _run([miniprot, "-t", str(threads), "-d", index, genome_fa],
             f"miniprot index ({label})")

    raw_gff = os.path.join(work_dir, f"{label}.raw.gff3")
    if not os.path.exists(raw_gff):
        with open(raw_gff, "w") as fh:
            r = _run([miniprot, "-t", str(threads), "--gff", "-P", prefix,
                      index, ref_pep], f"miniprot map ({label})")
            fh.write(r.stdout)

    # gffread chokes on ##PAF comment lines
    with open(raw_gff) as f, open(gff3, "w") as out:
        for line in f:
            if not line.startswith("##PAF"):
                out.write(line)

    _run([gffread, gff3, "-g", genome_fa, "-x", cds, "-y", pep],
         f"gffread ({label})")

    # CDS stop-trim: gffread includes the stop codon; keep only records
    # with len(cds) == 3 * len(pep) after trimming one terminal codon.
    pep_len = {}
    name = None
    n = 0
    with open(pep) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    pep_len[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name:
        pep_len[name] = n

    kept = set()
    with open(cds) as f, open(cds + ".fixed", "w") as out:
        seq, name = [], None

        def flush():
            if name is None:
                return
            plen = pep_len.get(name)
            if plen:
                s = "".join(seq)
                if len(s) == 3 * plen + 3:
                    s = s[:-3]
                if len(s) == 3 * plen:
                    out.write(f">{name}\n{s}\n")
                    kept.add(name)

        for line in f:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        flush()
    os.replace(cds + ".fixed", cds)

    # filter pep + bed to kept IDs; build bed from mRNA lines
    with open(pep) as f, open(pep + ".tmp", "w") as out:
        write = False
        for line in f:
            if line.startswith(">"):
                write = line[1:].split()[0] in kept
            if write:
                out.write(line)
    os.replace(pep + ".tmp", pep)

    with open(gff3) as f, open(bed, "w") as out:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 9 and p[2] == "mRNA":
                mid = p[8].split(";")[0].split("=")[1]
                if mid in kept:
                    out.write(f"{p[0]}\t{mid}\t{p[3]}\t{p[4]}\n")

    logger.info("  %s: %d genes annotated (pep/cds/bed consistent)",
                label, len(kept))
    return gff3, pep, cds, bed


# ============================================================================
# 2. Homeolog pair discovery (jcvi compara.catalog ortholog)
# ============================================================================

def _read_bed(bed_path: str) -> dict[str, tuple[str, int, int]]:
    pos = {}
    with open(bed_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                pos[p[1]] = (p[0], int(p[2]), int(p[3]))
    return pos


def _write_jcvi_bed(bed_path: str, out_path: str):
    """Convert our 4-col bed (scaffold gene start end) to jcvi 6-col bed."""
    with open(bed_path) as f, open(out_path, "w") as out:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                out.write(f"{p[0]}\t{int(p[2]) - 1}\t{p[3]}\t{p[1]}\t0\t+\n")


def find_homeolog_pairs(
    sp4_cds: str, sp2_cds: str, sp4_bed: str, sp2_bed: str,
    work_dir: str, threads: int = 4, min_block_pairs: int = 5,
    cscore: float = 0.7, min_size: int = 5,
    jcvi: str | None = None,
) -> dict[str, list[tuple[str, str, str]]]:
    """jcvi ortholog anchors -> 1v2 homeolog clusters.

    Runs ``python -m jcvi.compara.catalog ortholog sp2 sp4`` (LAST +
    jcvi's built-in synteny scan — replaces the old blastp + MCScanX
    wiring).  Returns {cluster_id: [(sp2, sp41, sp42), ...]} where each
    sp2 gene has exactly 2 distinct sp4 anchor hits, and only anchor
    blocks with > *min_block_pairs* such pairs are kept.

    *jcvi* overrides the catalog command prefix (default
    ``NEARSCAFF_JCVI`` env var or ``python -m jcvi.compara.catalog``).
    """
    import shlex

    jcvi = jcvi or os.environ.get("NEARSCAFF_JCVI",
                                  "python -m jcvi.compara.catalog")
    # jcvi works in its own subdir so that its {name}.bed inputs never
    # collide with our 4-col beds (same basename!) in work_dir.
    jdir = os.path.join(work_dir, "jcvi")
    os.makedirs(jdir, exist_ok=True)
    anchors = os.path.join(jdir, "sp2.sp4.anchors")

    if not os.path.exists(anchors):
        # jcvi expects {name}.cds + {name}.bed (6-col) in the cwd
        for name, cds, bed in (("sp4", sp4_cds, sp4_bed),
                               ("sp2", sp2_cds, sp2_bed)):
            dst = os.path.join(jdir, f"{name}.cds")
            if not os.path.exists(dst):
                os.link(cds, dst)
            _write_jcvi_bed(bed, os.path.join(jdir, f"{name}.bed"))
        cmd = shlex.split(jcvi) + [
            "ortholog", "sp2", "sp4", "--no_strip_names",
            "--cscore", str(cscore), "-n", str(min_size)]
        logger.info("  %s ...", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True,
                                cwd=jdir)
        if result.returncode != 0:
            raise RuntimeError(
                f"jcvi ortholog failed: {result.stderr.strip()[-500:]}")
    if not os.path.exists(anchors):
        raise RuntimeError(f"jcvi produced no anchors at {anchors}")

    # ---- parse anchors; orientation normalized by BED membership ----
    sp4_genes = set(_read_bed(sp4_bed))
    block_rows: dict[int, list[tuple[str, str]]] = defaultdict(list)
    block_id = -1
    with open(anchors) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("###"):
                block_id += 1
                continue
            p = line.split("\t")
            row = {}
            for g in p[0:2]:  # anchors: geneA \t geneB \t score
                if g in sp4_genes:
                    row["sp4"] = g
                else:
                    row["sp2"] = g
            if "sp4" in row and "sp2" in row:
                block_rows[block_id].append((row["sp2"], row["sp4"]))

    # 1v2: sp2 genes with exactly 2 distinct sp4 hits
    sp2_hits: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for block, rows in block_rows.items():
        for sp2, sp4 in rows:
            if sp4 not in {h for h, _b in sp2_hits[sp2]}:
                sp2_hits[sp2].append((sp4, block))
    one_v_two = {sp2: v for sp2, v in sp2_hits.items() if len(v) == 2}
    logger.info("  anchors: %d blocks, %d sp2 genes with 1v2 pairs",
                len(block_rows), len(one_v_two))

    # block filter: anchor blocks with > min_block_pairs 1v2 pairs
    block_count: dict[int, int] = defaultdict(int)
    for sp2, v in one_v_two.items():
        block_count[v[0][1]] += 1
    good_blocks = {b for b, c in block_count.items() if c > min_block_pairs}

    clusters: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for sp2, (h1, h2) in sorted(one_v_two.items(),
                                key=lambda kv: kv[1][0][1]):
        if h1[1] not in good_blocks:
            continue
        cid = f"cluster_{h1[1]}"
        clusters[cid].append((sp2, h1[0], h2[0]))
    logger.info("  %d clusters, %d 1v2 triples (blocks > %d pairs)",
                len(clusters), sum(len(v) for v in clusters.values()),
                min_block_pairs)
    return dict(clusters)


# ============================================================================
# 3. Cluster-side blocks
# ============================================================================

def cluster_side_blocks(
    clusters: dict[str, list[tuple[str, str, str]]],
    sp4_bed: str,
) -> dict[str, tuple[str, dict[str, tuple[str, int, int]]]]:
    """Per cluster: {"A": (scaffold, start, end), "B": (...)} side spans.

    Sides whose genes scatter across multiple scaffolds are dropped.
    Returns {cluster_id: {"A": (scaf, s, e), "B": (scaf, s, e)}} with both
    sides present.
    """
    pos = _read_bed(sp4_bed)
    out = {}
    for cid, triples in clusters.items():
        sides = {}
        ok = True
        for idx, side in ((1, "A"), (2, "B")):
            scafs = set()
            s_min, e_max = None, 0
            for t in triples:
                g = t[idx]
                if g not in pos:
                    continue
                sc, s, e = pos[g]
                scafs.add(sc)
                s_min = s if s_min is None else min(s_min, s)
                e_max = max(e_max, e)
            if len(scafs) != 1:
                ok = False
                break
            sides[side] = (scafs.pop(), s_min, e_max)
        if ok and sides.get("A") and sides.get("B"):
            out[cid] = sides
    logger.info("  %d clusters with clean single-scaffold side blocks",
                len(out))
    return out
