"""Unit tests for long-read gap filling (method "lr").

parse_* tests use synthetic PAF/FASTQ (no external tools);
overlap_closure_fills uses real minimap2 (skipped when unavailable).
"""
import os
import shutil

import pytest

from nearscaff.agp import AGPSeqLine, AGPGapLine, AGPReader
from nearscaff.gapfill import (scaffold_sequences, build_flank_fasta,
                               parse_span_fills, parse_tails,
                               overlap_closure_fills, _revcomp)

HAS_MINIMAP2 = shutil.which("minimap2") is not None

AGP_TEXT = "\n".join([
    "s1\t1\t8\t1\tW\tc1\t1\t8\t+",
    "s1\t9\t108\t2\tU\t100\tscaffold\tyes\tna",
    "s1\t109\t116\t3\tW\tc2\t1\t8\t+",
    "s1\t117\t216\t4\tU\t100\tscaffold\tyes\tna",
    "s1\t217\t224\t5\tW\tc3\t1\t8\t+",
    "s1\t225\t324\t6\tU\t100\tscaffold\tyes\tna",
    "s1\t325\t332\t7\tW\tc4\t1\t8\t-",
]) + "\n"

CONTIGS = {"c1": "AAAAACGT", "c2": "CCCCGGGG", "c3": "TTTTAAAA",
           "c4": "AAACCCGG"}  # RC(c4) = CCGGGTTT


def _agp():
    return AGPReader().parse(AGP_TEXT)


def test_scaffold_sequences_orientation_and_gaps():
    scafs = scaffold_sequences(_agp(), CONTIGS)
    assert scafs == {"s1": "AAAAACGT" + "N" * 100 + "CCCCGGGG" + "N" * 100
                          + "TTTTAAAA" + "N" * 100 + "CCGGGTTT"}


def test_build_flank_fasta_coords():
    scafs = scaffold_sequences(_agp(), CONTIGS)
    records, flank_len = build_flank_fasta(_agp(), {1}, scafs, flank=8)
    d = dict(records)
    gid = "s1:9-108"
    # gap 1 spans [8,108) 0-based; L flank = [0,8), R flank = [108,116)
    assert d[f"{gid}|L"] == "AAAAACGT"
    assert d[f"{gid}|R"] == "CCCCGGGG"
    assert flank_len[f"{gid}|L"] == 8


def test_build_flank_fasta_drops_n_rich():
    scafs = scaffold_sequences(_agp(), CONTIGS)
    # gap 3's R flank (after c3, before c4) is short contig + long N-run:
    # with flank=200 the R flank = "TTTTAAAA" + 100 N -> mostly Ns, dropped
    records, _ = build_flank_fasta(_agp(), {3}, scafs, flank=200,
                                   max_n_frac=0.5)
    names = [n for n, _ in records]
    assert not any(n.endswith("|R") for n in names)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# ---- span / tail parsing with synthetic PAF + FASTQ -----------------------

# flank records: gid|L len 100, gid|R len 100
FLANK_LEN = {"g1|L": 100, "g1|R": 100}

READS_FQ = ("@r1\n" + "A" * 50 + "FILLSEQACGT" + "C" * 100 + "\n+\n" + "I" * 161 + "\n"
            "@r2\n" + _revcomp("A" * 50 + "REVSEQTTTT" + "C" * 100) + "\n+\n" + "I" * 160 + "\n"
            "@r3\n" + "C" * 100 + "TAILSEQACGT" + "\n+\n" + "I" * 112 + "\n")


def test_parse_span_fills_fwd(tmp_path):
    # r1 '+' (161 bp): L aln [0,50) on g1|L ending at 100; R aln [61,161) on g1|R
    paf = _write(tmp_path, "x.paf",
                 "r1\t161\t0\t50\t+\tg1|L\t100\t50\t100\t50\t50\t60\n"
                 "r1\t161\t61\t161\t+\tg1|R\t100\t0\t100\t100\t100\t60\n")
    fq = _write(tmp_path, "r.fq", READS_FQ)
    fills = parse_span_fills(paf, fq)
    assert fills == {"g1": "FILLSEQACGT"}


def test_parse_span_fills_rev(tmp_path):
    # r2 '-' (160 bp): R-flank aln [0,100), L-flank aln [110,160)
    paf = _write(tmp_path, "x.paf",
                 "r2\t160\t0\t100\t-\tg1|R\t100\t0\t100\t100\t100\t60\n"
                 "r2\t160\t110\t160\t-\tg1|L\t100\t50\t100\t50\t50\t60\n")
    fq = _write(tmp_path, "r.fq", READS_FQ)
    fills = parse_span_fills(paf, fq)
    assert fills == {"g1": "REVSEQTTTT"}


def test_parse_span_fills_min_reads(tmp_path):
    paf = _write(tmp_path, "x.paf",
                 "r1\t161\t0\t50\t+\tg1|L\t100\t50\t100\t50\t50\t60\n"
                 "r1\t161\t61\t161\t+\tg1|R\t100\t0\t100\t100\t100\t60\n")
    fq = _write(tmp_path, "r.fq", READS_FQ)
    assert parse_span_fills(paf, fq, min_reads=2) == {}


def test_parse_span_fills_edge_shortfall_and_sanity(tmp_path):
    """Alignments stopping >max_shortfall before the gap edge are rejected;
    coordinate sanity (gap end after gap start) is enforced."""
    fq = _write(tmp_path, "r.fq", READS_FQ)
    # shortfall 400 bp > default 300 -> rejected (flank len 500, aln ends 100)
    paf = _write(tmp_path, "x.paf",
                 "r1\t161\t0\t100\t+\tg1|L\t500\t400\t500\t100\t100\t60\n"
                 "r1\t161\t61\t161\t+\tg1|R\t500\t0\t100\t100\t100\t60\n")
    assert parse_span_fills(paf, fq) == {}
    # alignment ends at edge but gap end before gap start -> rejected
    paf = _write(tmp_path, "y.paf",
                 "r1\t161\t0\t50\t+\tg1|L\t100\t50\t100\t50\t50\t60\n"
                 "r1\t161\t61\t161\t+\tg1|R\t100\t50\t100\t50\t50\t60\n")
    assert parse_span_fills(paf, fq) == {}


def test_parse_span_fills_extrapolates_shortfall(tmp_path):
    """Alignment ending a bit before the edge: fill starts at the gap
    boundary (extrapolated into the unaligned flank remainder), not at
    the alignment end."""
    # read: A*50 | F*50 | FILLSEQACGT | C*100  (211 bp)
    # L aln covers read [0,50) on target [100,150) of a 200 bp flank
    # -> shortfall 50; gap boundary at read 50+50=100
    # R aln covers read [111,211) from target 0 -> boundary at read 111
    read = "A" * 50 + "F" * 50 + "FILLSEQACGT" + "C" * 100
    fq = _write(tmp_path, "r.fq",
                f"@rX\n{read}\n+\n{'I' * len(read)}\n")
    paf = _write(tmp_path, "x.paf",
                 "rX\t211\t0\t50\t+\tg1|L\t200\t100\t150\t50\t50\t60\n"
                 "rX\t211\t111\t211\t+\tg1|R\t200\t0\t100\t100\t100\t60\n")
    fills = parse_span_fills(paf, fq)
    assert fills == {"g1": "FILLSEQACGT"}


def test_parse_span_fills_depth_gate(tmp_path):
    """Over-collapsed gaps (depth >> median) are rejected; normal ones close."""
    lines = []
    reads = []
    # g_normal: 2 spanning reads; g_hot: 30 spanning reads -> rejected
    for gid, n in (("gN", 2), ("gH", 30)):
        for i in range(n):
            rn = f"{gid}_{i}"
            read = "A" * 50 + "FILLSEQACGT" + "C" * 100
            reads.append(f"@{rn}\n{read}\n+\n{'I' * len(read)}\n")
            lines.append(f"{rn}\t161\t0\t50\t+\t{gid}|L\t100\t50\t100\t50\t50\t60")
            lines.append(f"{rn}\t161\t61\t161\t+\t{gid}|R\t100\t0\t100\t100\t100\t60")
    paf = _write(tmp_path, "d.paf", "\n".join(lines) + "\n")
    fq = _write(tmp_path, "d.fq", "".join(reads))
    fills = parse_span_fills(paf, fq, max_depth_factor=3.0)
    assert "gN" in fills and "gH" not in fills
    # gate disabled -> both close
    fills = parse_span_fills(paf, fq, max_depth_factor=None)
    assert "gN" in fills and "gH" in fills


def test_parse_tails_directions(tmp_path):
    # Tails are stored in genomic-forward orientation: L tails run
    # flank -> deep gap, R tails run deep gap -> flank.
    # L, '+' read: right clip -> tail = seq[qe:]
    # R, '+' read: left clip is already genomic-forward -> tail = seq[:qs]
    # R, '-' read: right clip -> tail = rc(seq[qe:])
    paf = _write(tmp_path, "x.paf",
                 "r1\t161\t0\t50\t+\tg1|L\t100\t50\t100\t50\t50\t60\n"
                 "r1\t161\t61\t161\t+\tg1|R\t100\t0\t100\t100\t100\t60\n"
                 "r3\t112\t0\t100\t-\tg1|R\t100\t0\t100\t100\t100\t60\n")
    fq = _write(tmp_path, "r.fq", READS_FQ)
    tails = parse_tails(paf, FLANK_LEN, fq, min_clip=10)
    assert tails["g1"]["L"] == ["FILLSEQACGT" + "C" * 100]
    assert tails["g1"]["R"] == ["A" * 50 + "FILLSEQACGT",
                                _revcomp("TAILSEQACGT")]


# ---- overlap closure (real minimap2) --------------------------------------

@pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
def test_overlap_closure_unique_accepted(tmp_path):
    left = "GATTACA" * 100          # 700 bp unique-ish
    middle = "ACGTTGCA" * 100       # 800 bp (overlap)
    right = "TTGCAGTC" * 100        # 800 bp
    tails = {"g1": {"L": [left + middle], "R": [middle + right]}}
    scaf_seqs = {"s1": "N" * 50 + left + middle + right + "N" * 50}
    fills = overlap_closure_fills(tails, scaf_seqs, str(tmp_path / "w"),
                                  threads=2, min_ovlp=200, min_ident=0.8)
    assert fills.get("g1") == left + middle + right


@pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
def test_overlap_closure_repeat_rejected(tmp_path):
    rep = "ACGTTGCA" * 100
    tails = {"g1": {"L": ["GATTACA" * 100 + rep],
                    "R": [rep + "TTGCAGTC" * 100]}}
    # overlap segment appears TWICE in the scaffolds -> rejected
    scaf_seqs = {"s1": "N" * 50 + rep + "N" * 100 + rep + "N" * 50}
    fills = overlap_closure_fills(tails, scaf_seqs, str(tmp_path / "w"),
                                  threads=2, min_ovlp=200, min_ident=0.8)
    assert fills == {}
