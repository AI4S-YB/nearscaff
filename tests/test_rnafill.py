"""Unit tests for transcriptome-guided gap filling (rna-fill).

parse_* tests use synthetic PAF/FASTA (no external tools).
"""
import gzip

from nearscaff.gapfill import _fetch_read_seqs, _revcomp
from nearscaff.rnafill import (_canonical_splice_edges, _gid_coords,
                               parse_tx_span_fills, TX_ABUT_WINDOW)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# ---- sequence fetching -----------------------------------------------------

def test_fetch_read_seqs_fasta_and_fastq(tmp_path):
    fa = _write(tmp_path, "t.fa", ">tx1 some comment\nACGT\nACGT\n>tx2\n"
                + "A" * 100 + "\n")
    fq = _write(tmp_path, "t.fq", "@r1\nTTTT\n+\nIIII\n")
    gz = tmp_path / "t.fa.gz"
    with gzip.open(gz, "wt") as f:
        f.write(">tx3\n" + "G" * 50 + "\n")
    seqs = _fetch_read_seqs([fa, fq, str(gz)], {"tx1", "tx2", "r1", "tx3"})
    assert seqs == {"tx1": "ACGTACGT", "tx2": "A" * 100,
                    "r1": "TTTT", "tx3": "G" * 50}


# ---- splice-edge heuristic --------------------------------------------------

def test_canonical_splice_edges():
    # gap [10, 20) in a 30 bp scaffold; donor GT at left edge, AG at right
    scaf = "AAAAAAAA" + "GT" + "N" * 10 + "AG" + "CCCCCCCC"
    assert _canonical_splice_edges(scaf, 10, 20)
    # minus strand: CT .. AC
    scaf_m = "AAAAAAAA" + "CT" + "N" * 10 + "AC" + "CCCCCCCC"
    assert _canonical_splice_edges(scaf_m, 10, 20)
    # no signal
    scaf_n = "AAAAAAAA" + "AA" + "N" * 10 + "CC" + "CCCCCCCC"
    assert not _canonical_splice_edges(scaf_n, 10, 20)


def test_gid_coords():
    assert _gid_coords("Chr1:2001-2100") == ("Chr1", 2000, 2100)


# ---- transcript span fills --------------------------------------------------

# flank records: g1|L len 2000, g1|R len 2000
MID = "MIDSEQACGT" * 6          # 60 bp middle segment
TX1 = "A" * 500 + MID + "C" * 500           # spans with 60 bp middle
TX2 = "G" * 500 + "T" * 500                 # spans with 0 bp middle (abut)
TX3 = _revcomp("A" * 500 + "REVSEQTTTT" * 5 + "C" * 500)  # '-' strand

TX_FA = (f">tx1\n{TX1}\n>tx2\n{TX2}\n>tx3\n{TX3}\n")

SCAF_DUMMY = {"s1": "N" * 5000}
# g2 gap with canonical GT..AG edges: flank L ends ...GT, flank R starts AG...
SCAF_SPLICE = {"s2": "A" * 1998 + "GT" + "N" * 100 + "AG" + "C" * 1998}


def _paf(tmp_path, lines):
    return _write(tmp_path, "x.paf", "\n".join(lines) + "\n")


def _fa(tmp_path, text=TX_FA):
    return _write(tmp_path, "t.fa", text)


def test_tx_span_fill_exonic(tmp_path):
    """Transcript hitting both flanks with a 60 bp middle -> fill = middle."""
    paf = _paf(tmp_path, [
        f"tx1\t{len(TX1)}\t0\t500\t+\tg1|L\t2000\t1500\t2000\t500\t500\t60",
        f"tx1\t{len(TX1)}\t560\t1060\t+\tg1|R\t2000\t0\t500\t500\t500\t60",
    ])
    fills, stats = parse_tx_span_fills(paf, _fa(tmp_path), SCAF_DUMMY)
    assert fills == {"g1": MID}
    assert stats["span_filled"] == 1 and stats["abut"] == 0


def test_tx_span_fill_revstrand(tmp_path):
    """'-' transcript: middle extracted and revcomped back to genomic fwd."""
    mid = "REVSEQTTTT" * 5   # 50 bp
    paf = _paf(tmp_path, [
        f"tx3\t{len(TX3)}\t0\t500\t-\tg1|R\t2000\t0\t500\t500\t500\t60",
        f"tx3\t{len(TX3)}\t550\t1050\t-\tg1|L\t2000\t1500\t2000\t500\t500\t60",
    ])
    fills, stats = parse_tx_span_fills(paf, _fa(tmp_path), SCAF_DUMMY)
    assert fills == {"g1": mid}


def test_tx_span_abut_cdna_vs_exononly(tmp_path):
    """Zero-length middle (adjacent exons): abut in cdna, skipped in
    exon-only."""
    lines = [
        f"tx2\t{len(TX2)}\t0\t500\t+\tg1|L\t2000\t1500\t2000\t500\t500\t60",
        f"tx2\t{len(TX2)}\t500\t1000\t+\tg1|R\t2000\t0\t500\t500\t500\t60",
    ]
    paf = _paf(tmp_path, lines)
    fills, stats = parse_tx_span_fills(paf, _fa(tmp_path), SCAF_DUMMY,
                                       fill_mode="cdna")
    assert fills == {"g1": ""}
    assert stats["abut"] == 1
    fills, stats = parse_tx_span_fills(paf, _fa(tmp_path), SCAF_DUMMY,
                                       fill_mode="exon-only")
    assert fills == {}
    assert stats["intron_skipped"] == 1


def test_tx_span_exononly_splice_edge_rejected(tmp_path):
    """Middle > abut-window but gap edges are GT..AG -> intronic, skipped
    in exon-only mode, kept in cdna mode."""
    lines = [
        f"tx1\t{len(TX1)}\t0\t500\t+\ts2:2001-2100|L\t2000\t1500\t2000\t500\t500\t60",
        f"tx1\t{len(TX1)}\t560\t1060\t+\ts2:2001-2100|R\t2000\t0\t500\t500\t500\t60",
    ]
    paf = _paf(tmp_path, lines)
    fills, stats = parse_tx_span_fills(paf, _fa(tmp_path), SCAF_SPLICE,
                                       fill_mode="exon-only")
    assert fills == {}
    assert stats["intron_skipped"] == 1
    fills, stats = parse_tx_span_fills(paf, _fa(tmp_path), SCAF_SPLICE,
                                       fill_mode="cdna")
    assert fills == {"s2:2001-2100": MID}


def test_tx_span_median_candidate(tmp_path):
    """Multiple isoforms per gap: the median-length candidate wins."""
    mids = ["M" * 60, "G" * 70, "Q" * 80]
    reads = []
    lines = []
    for i, mid in enumerate(mids):
        rn = f"iso{i}"
        read = "A" * 500 + mid + "C" * 500
        reads.append(f">{rn}\n{read}\n")
        r_end = 500 + len(mid)
        lines.append(f"{rn}\t{len(read)}\t0\t500\t+\tg1|L\t2000\t1500\t2000\t500\t500\t60")
        lines.append(f"{rn}\t{len(read)}\t{r_end}\t{r_end + 500}\t+\tg1|R\t2000\t0\t500\t500\t500\t60")
    paf = _paf(tmp_path, lines)
    fills, stats = parse_tx_span_fills(paf, _fa(tmp_path, "".join(reads)),
                                       SCAF_DUMMY)
    assert fills == {"g1": "G" * 70}
    assert stats["multi_candidate_gaps"] == 1


def test_tx_span_depth_gate(tmp_path):
    """Over-collapsed gaps (transcript depth >> median) are rejected."""
    lines = []
    reads = []
    for gid, n in (("gN", 2), ("gH", 30)):
        for i in range(n):
            rn = f"{gid}_{i}"
            read = "A" * 500 + MID + "C" * 500
            reads.append(f">{rn}\n{read}\n")
            lines.append(f"{rn}\t{len(read)}\t0\t500\t+\t{gid}|L\t2000\t1500\t2000\t500\t500\t60")
            lines.append(f"{rn}\t{len(read)}\t560\t1060\t+\t{gid}|R\t2000\t0\t500\t500\t500\t60")
    paf = _paf(tmp_path, lines)
    fills, stats = parse_tx_span_fills(paf, _fa(tmp_path, "".join(reads)),
                                       SCAF_DUMMY, max_depth_factor=3.0)
    assert "gN" in fills and "gH" not in fills
    assert stats["depth_rejected"] == 1


def test_tx_span_prefers_primary_per_side(tmp_path):
    """A secondary hit recorded first must not block the primary hit of
    the same transcript on the same flank (cross-locus guard)."""
    fa = _fa(tmp_path)
    # tx1's L flank: a secondary hit recorded first would extrapolate to
    # a spurious abut; the primary hit (arriving later) must replace it
    paf = _paf(tmp_path, [
        f"tx1\t{len(TX1)}\t0\t400\t+\tg1|L\t2000\t1400\t1800\t400\t400\t60\ttp:A:S",
        f"tx1\t{len(TX1)}\t0\t500\t+\tg1|L\t2000\t1500\t2000\t500\t500\t60\ttp:A:P",
        f"tx1\t{len(TX1)}\t560\t1060\t+\tg1|R\t2000\t0\t500\t500\t500\t60\ttp:A:P",
    ])
    fills, _ = parse_tx_span_fills(paf, fa, SCAF_DUMMY)
    assert fills == {"g1": TX1[500:560]}


def test_tx_span_shortfall_and_sanity(tmp_path):
    """Alignments stopping >max_shortfall before the edge are rejected;
    a tolerable shortfall is extrapolated (fill starts at the gap
    boundary, not at the alignment end)."""
    fa = _fa(tmp_path)
    # L alignment stops 400 bp short of the edge (> default 300)
    paf = _paf(tmp_path, [
        f"tx1\t{len(TX1)}\t0\t400\t+\tg1|L\t2000\t1200\t1600\t400\t400\t60",
        f"tx1\t{len(TX1)}\t560\t1060\t+\tg1|R\t2000\t0\t500\t500\t500\t60",
    ])
    fills, _ = parse_tx_span_fills(paf, fa, SCAF_DUMMY)
    assert fills == {}
    # shortfall 100 <= 300: boundary extrapolated 100 bp into the
    # transcript, so the fill excludes the 100 bp of flank sequence
    paf = _paf(tmp_path, [
        f"tx1\t{len(TX1)}\t0\t400\t+\tg1|L\t2000\t1500\t1900\t400\t400\t60",
        f"tx1\t{len(TX1)}\t560\t1060\t+\tg1|R\t2000\t0\t500\t500\t500\t60",
    ])
    fills, _ = parse_tx_span_fills(paf, fa, SCAF_DUMMY)
    assert fills == {"g1": TX1[500:560]}


# ---- one-sided extensions ---------------------------------------------------

from nearscaff.rnafill import parse_tx_extensions  # noqa: E402
from nearscaff.agp import AGPReader, AGPSeqLine, AGPGapLine  # noqa: E402
from nearscaff.gapfill import write_outputs  # noqa: E402

TXL = "A" * 500 + "E" * 800          # L flank hit + 800 bp 3' clip
TXR = "R" * 700 + "C" * 500          # 700 bp 5' clip + R flank hit
EXT_FA = f">txL\n{TXL}\n>txR\n{TXR}\n"


def test_tx_extensions_both_sides(tmp_path):
    fa = _write(tmp_path, "t.fa", EXT_FA)
    paf = _paf(tmp_path, [
        f"txL\t{len(TXL)}\t0\t500\t+\tg1|L\t2000\t1500\t2000\t500\t500\t60\ttp:A:P",
        f"txR\t{len(TXR)}\t700\t1200\t+\tg1|R\t2000\t0\t500\t500\t500\t60\ttp:A:P",
    ])
    ext = parse_tx_extensions(paf, fa, min_tail=500)
    assert ext == {"g1": ("E" * 800, "R" * 700, "txL")}


def test_tx_extensions_primary_only(tmp_path):
    fa = _write(tmp_path, "t.fa", EXT_FA)
    # same L hit but flagged secondary -> rejected
    paf = _paf(tmp_path, [
        f"txL\t{len(TXL)}\t0\t500\t+\tg1|L\t2000\t1500\t2000\t500\t500\t60\ttp:A:S",
    ])
    assert parse_tx_extensions(paf, fa, min_tail=500) == {}


def test_tx_extensions_min_tail_and_truncation(tmp_path):
    fa = _write(tmp_path, "t.fa", EXT_FA)
    paf = _paf(tmp_path, [
        f"txL\t{len(TXL)}\t0\t500\t+\tg1|L\t2000\t1500\t2000\t500\t500\t60\ttp:A:P",
    ])
    # min_tail above the clip length -> rejected
    assert parse_tx_extensions(paf, fa, min_tail=1000) == {}
    # max_tail caps the deep end (L tail keeps the flank-proximal part)
    ext = parse_tx_extensions(paf, fa, min_tail=500, max_tail=600)
    assert ext["g1"][0] == "E" * 600


def test_write_outputs_extensions(tmp_path):
    """One-sided extension: gap stays open, known sequence written next
    to its flank; part numbers stay consecutive."""
    contigs = {"c1": "AAAAACGT", "c2": "CCCCGGGG"}
    agp_text = "\n".join([
        "s1\t1\t8\t1\tW\tc1\t1\t8\t+",
        "s1\t9\t108\t2\tU\t100\tscaffold\tyes\tna",
        "s1\t109\t116\t3\tW\tc2\t1\t8\t+",
    ]) + "\n"
    agp_lines = AGPReader().parse(agp_text)
    out_agp = str(tmp_path / "o.agp")
    out_fa = str(tmp_path / "o.fa")
    report = write_outputs(agp_lines, {1}, {}, contigs, out_agp, out_fa,
                           extensions={1: ("LLL", "RRR")})
    assert report["gaps_closed"] == 0
    assert report["gaps_extended"] == 1
    assert report["bases_filled"] == 6
    lines = AGPReader().parse(open(out_agp).read())
    assert [l.part_number for l in lines] == list(range(1, len(lines) + 1))
    w = [l for l in lines if isinstance(l, AGPSeqLine)]
    g = [l for l in lines if isinstance(l, AGPGapLine)]
    assert len(w) == 4 and len(g) == 1        # c1, LLL, gap, RRR, c2
    assert (w[1].object_beg, w[1].object_end) == (9, 11)
    assert (g[0].object_beg, g[0].object_end) == (12, 111)
    assert (w[2].object_beg, w[2].object_end) == (112, 114)
    from nearscaff.gapfill import read_fasta as _rf
    fa = _rf(out_fa)
    assert fa == {"s1": "AAAAACGT" + "LLL" + "N" * 100 + "RRR" + "CCCCGGGG"}


# ---- internal N-run explosion ------------------------------------------------

from nearscaff.rnafill import explode_internal_ns  # noqa: E402
from nearscaff.gapfill import scaffold_sequences  # noqa: E402


def _explode_test_contig(orient="+"):
    # 1000 bp component with a 60 bp N run at oriented [400, 460)
    left = "A" * 400
    right = "C" * 540
    if orient == "+":
        seq = left + "N" * 60 + right
    else:
        seq = _revcomp(left + "N" * 60 + right)
    return seq


def test_explode_internal_ns_plus():
    contigs = {"c1": _explode_test_contig("+")}
    agp_lines = AGPReader().parse("s1\t1\t1000\t1\tW\tc1\t1\t1000\t+\n")
    out, n = explode_internal_ns(agp_lines, contigs, min_len=20, min_edge=50)
    assert n == 1
    assert len(out) == 3
    w = [l for l in out if isinstance(l, AGPSeqLine)]
    g = [l for l in out if isinstance(l, AGPGapLine)]
    assert len(w) == 2 and len(g) == 1
    # scaffold coords
    assert (w[0].object_beg, w[0].object_end) == (1, 400)
    assert (g[0].object_beg, g[0].object_end) == (401, 460)
    assert (w[1].object_beg, w[1].object_end) == (461, 1000)
    assert g[0].gap_length == 60
    # component coords ('+')
    assert (w[0].component_beg, w[0].component_end) == (1, 400)
    assert (w[1].component_beg, w[1].component_end) == (461, 1000)
    assert [l.part_number for l in out] == [1, 2, 3]
    # the exploded AGP must rebuild the identical scaffold sequence
    assert scaffold_sequences(out, contigs)["s1"] == contigs["c1"]


def test_explode_internal_ns_minus_strand():
    contigs = {"c1": _explode_test_contig("-")}
    agp_lines = AGPReader().parse("s1\t1\t1000\t1\tW\tc1\t1\t1000\t-\n")
    out, n = explode_internal_ns(agp_lines, contigs, min_len=20, min_edge=50)
    assert n == 1
    w = [l for l in out if isinstance(l, AGPSeqLine)]
    g = [l for l in out if isinstance(l, AGPGapLine)]
    # scaffold coords are orientation-independent
    assert (g[0].object_beg, g[0].object_end) == (401, 460)
    # component coords are mirrored: oriented [0,400) <- comp [601,1000]
    assert (w[0].component_beg, w[0].component_end) == (601, 1000)
    assert (w[1].component_beg, w[1].component_end) == (1, 540)
    assert scaffold_sequences(out, contigs)["s1"] == _revcomp(
        contigs["c1"])


def test_explode_internal_ns_filters():
    # run shorter than min_len: untouched
    contigs = {"c1": "A" * 400 + "N" * 10 + "C" * 590}
    agp_lines = AGPReader().parse("s1\t1\t1000\t1\tW\tc1\t1\t1000\t+\n")
    out, n = explode_internal_ns(agp_lines, contigs, min_len=20)
    assert n == 0 and len(out) == 1
    # run too close to the component edge (< min_edge on one side)
    contigs = {"c1": "A" * 30 + "N" * 60 + "C" * 910}
    out, n = explode_internal_ns(agp_lines, contigs, min_len=20, min_edge=50)
    assert n == 0 and len(out) == 1


def test_explode_internal_ns_multiple_runs():
    contigs = {"c1": "A" * 200 + "N" * 30 + "G" * 300 + "N" * 40 +
               "C" * 460}
    agp_lines = AGPReader().parse("s1\t1\t1030\t1\tW\tc1\t1\t1030\t+\n")
    out, n = explode_internal_ns(agp_lines, contigs, min_len=20, min_edge=50)
    assert n == 2
    w = [l for l in out if isinstance(l, AGPSeqLine)]
    g = [l for l in out if isinstance(l, AGPGapLine)]
    assert len(w) == 3 and len(g) == 2
    assert (g[0].object_beg, g[0].object_end) == (201, 230)
    assert (g[1].object_beg, g[1].object_end) == (531, 570)
    assert (w[1].component_beg, w[1].component_end) == (231, 530)
    assert [l.part_number for l in out] == [1, 2, 3, 4, 5]
    assert scaffold_sequences(out, contigs)["s1"] == contigs["c1"]


# ---- read recruitment (broken-gene bait) ------------------------------------

from nearscaff.rnafill import find_broken_loci, extract_read_pairs  # noqa: E402


def test_find_broken_loci(tmp_path):
    # two alignments: p1 translation contains X, p2 clean; p1 has a second
    # (secondary) record that must not create a duplicate locus
    trans = _write(tmp_path, "t.faa",
                   "p1\t100\t0\t100\t+\ts1\t1000\t100\t200\t300\t300\t0\n"
                   "##STA\tMMAAXXXAA*\n"
                   "p2\t100\t0\t100\t+\ts1\t1000\t500\t600\t300\t300\t0\n"
                   "##STA\tMMAAAAAAA*\n"
                   "p1\t100\t0\t100\t+\ts1\t1000\t700\t800\t300\t300\t0\n"
                   "##STA\tMMAAXXXAA*\n")
    gff = _write(tmp_path, "a.gff",
                 "s1\tmp\tmRNA\t100\t200\t.\t+\t.\tID=m1;Target=p1 1 100\n"
                 "s1\tmp\tmRNA\t500\t600\t.\t+\t.\tID=m2;Target=p2 1 100\n"
                 "s1\tmp\tmRNA\t700\t800\t.\t+\t.\tID=m3;Target=p1 1 100\n")
    loci = find_broken_loci(gff, trans, {"s1": 5000}, pad=50)
    # m1 (50-250 after pad) and m3 (650-850) — both from p1, not merged
    assert loci == {"s1": [(50, 250), (650, 850)]}


def test_find_broken_loci_merge_and_clamp(tmp_path):
    trans = _write(tmp_path, "t.faa",
                   "p1\t10\t0\t10\t+\ts1\t100\t1\t10\t30\t30\t0\n"
                   "##STA\tXX*\n")
    gff = _write(tmp_path, "a.gff",
                 "s1\tmp\tmRNA\t5\t100\t.\t+\t.\tID=m1;Target=p1 1 10\n"
                 "s1\tmp\tmRNA\t90\t200\t.\t+\t.\tID=m2;Target=p1 1 10\n")
    # pad 50: [-45..150] and [40..250] merge and clamp to [1, 250]
    loci = find_broken_loci(gff, trans, {"s1": 300}, pad=50)
    assert loci == {"s1": [(1, 250)]}


def test_extract_read_pairs(tmp_path):
    paf = _write(tmp_path, "r.paf",
                 "r1\t100\t0\t100\t+\tbait\t1000\t0\t100\t100\t100\t60\n"
                 "r3\t100\t0\t50\t+\tbait\t1000\t0\t50\t50\t50\t60\n"
                 "unmapped\t*\n")
    fq1 = _write(tmp_path, "r1.fq",
                 "@r1 1:N:0\nAAAA\n+\nIIII\n@r2 1:N:0\nCCCC\n+\nIIII\n"
                 "@r3 1:N:0\nGGGG\n+\nIIII\n")
    fq2 = _write(tmp_path, "r2.fq",
                 "@r1 2:N:0\nTTTT\n+\nIIII\n@r2 2:N:0\nAAAA\n+\nIIII\n"
                 "@r3 2:N:0\nCCCC\n+\nIIII\n")
    o1 = str(tmp_path / "o1.fq")
    o2 = str(tmp_path / "o2.fq")
    n = extract_read_pairs(paf, [fq1, fq2], [o1, o2])
    assert n == 2                       # r1, r3 (unmapped line skipped)
    assert open(o1).read().count("@") == 2
    assert "@r2" not in open(o1).read()
    assert "@r3" in open(o2).read()


# ---- ref-guided whole-gene placement -----------------------------------------

import random  # noqa: E402
import shutil  # noqa: E402

import pytest  # noqa: E402

from nearscaff.rnafill import (  # noqa: E402
    read_ref_genes, contig_ref_index, gap_ref_brackets,
    assign_tx_to_genes, build_gene_placement, plan_placements,
    run_rnafill)

HAS_MINIMAP2 = shutil.which("minimap2") is not None


def _refgenes_text():
    return ("chrR\tsrc\tmRNA\t2001\t2400\t.\t+\t.\tID=g1\n"
            "chrR\tsrc\tmRNA\t6001\t6500\t.\t-\t.\tID=g2\n")


def test_read_ref_genes_mrna_and_fallback(tmp_path):
    gff = _write(tmp_path, "r.gff", _refgenes_text())
    genes = read_ref_genes(gff)
    assert genes == {"chrR": [(2000, 2400, "g1", "+"),
                              (6000, 6500, "g2", "-")]}
    # no mRNA lines -> fall back to gene lines
    gff2 = _write(tmp_path, "r2.gff",
                  "chrR\tsrc\tgene\t101\t200\t.\t+\t.\tID=gx\n")
    assert read_ref_genes(gff2) == {"chrR": [(100, 200, "gx", "+")]}


def _tiered_paf(tmp_path, lines):
    return _write(tmp_path, "tiered.paf", "\n".join(lines) + "\n")


def test_contig_ref_index_and_brackets(tmp_path):
    # cA -> chrR [10000,15000) '+', cB -> chrR [15700,20700) '+'
    paf = _tiered_paf(tmp_path, [
        "cA\t5000\t0\t5000\t+\tchrR\t30000\t10000\t15000\t5000\t5000\t60"
        "\ttp:A:P\tnc:Z:protein",
        "cB\t5000\t0\t5000\t+\tchrR\t30000\t15700\t20700\t5000\t5000\t60"
        "\ttp:A:P\tnc:Z:protein",
    ])
    idx = contig_ref_index(paf)
    assert idx["cA"][0][:6] == (0, 5000, "+", "chrR", 10000, 15000)
    agp_lines = AGPReader().parse(
        "s1\t1\t5000\t1\tW\tcA\t1\t5000\t+\n"
        "s1\t5001\t5100\t2\tU\t100\tscaffold\tyes\tna\n"
        "s1\t5101\t10100\t3\tW\tcB\t1\t5000\t+\n")
    br = gap_ref_brackets(agp_lines, {1}, idx)
    assert br == {1: ("chrR", 14999, 15700, False)}


def test_gap_ref_brackets_minus_strand(tmp_path):
    # both flanks on '-' strand: bracket identical, flip=True
    paf = _tiered_paf(tmp_path, [
        "cA\t5000\t0\t5000\t-\tchrR\t30000\t10000\t15000\t5000\t5000\t60"
        "\ttp:A:P",
        "cB\t5000\t0\t5000\t-\tchrR\t30000\t15700\t20700\t5000\t5000\t60"
        "\ttp:A:P",
    ])
    idx = contig_ref_index(paf)
    agp_lines = AGPReader().parse(
        "s1\t1\t5000\t1\tW\tcA\t1\t5000\t-\n"
        "s1\t5001\t5100\t2\tU\t100\tscaffold\tyes\tna\n"
        "s1\t5101\t10100\t3\tW\tcB\t1\t5000\t-\n")
    br = gap_ref_brackets(agp_lines, {1}, idx)
    assert br == {1: ("chrR", 14999, 15700, True)}


def test_gap_ref_brackets_rejects(tmp_path):
    agp_lines = AGPReader().parse(
        "s1\t1\t5000\t1\tW\tcA\t1\t5000\t+\n"
        "s1\t5001\t5100\t2\tU\t100\tscaffold\tyes\tna\n"
        "s1\t5101\t10100\t3\tW\tcB\t1\t5000\t+\n")
    # different chromosomes -> no bracket
    paf = _tiered_paf(tmp_path, [
        "cA\t5000\t0\t5000\t+\tchrR\t30000\t10000\t15000\t5000\t5000\t60"
        "\ttp:A:P",
        "cB\t5000\t0\t5000\t+\tchrX\t30000\t15700\t20700\t5000\t5000\t60"
        "\ttp:A:P",
    ])
    assert gap_ref_brackets(agp_lines, {1}, contig_ref_index(paf)) == {}
    # bracket too large -> no bracket
    paf2 = _tiered_paf(tmp_path, [
        "cA\t5000\t0\t5000\t+\tchrR\t9000000\t10000\t15000\t5000\t5000\t60"
        "\ttp:A:P",
        "cB\t5000\t0\t5000\t+\tchrR\t9000000\t15700\t20700\t5000\t5000\t60"
        "\ttp:A:P",
    ])
    assert gap_ref_brackets(agp_lines, {1}, contig_ref_index(paf2),
                            max_bracket=500) == {}
    # secondary alignment used when no primary covers the coordinate
    paf3 = _tiered_paf(tmp_path, [
        "cA\t5000\t0\t5000\t+\tchrR\t30000\t10000\t15000\t5000\t5000\t60"
        "\ttp:A:S",
        "cB\t5000\t0\t5000\t+\tchrR\t30000\t15700\t20700\t5000\t5000\t60"
        "\ttp:A:S",
    ])
    br = gap_ref_brackets(agp_lines, {1}, contig_ref_index(paf3))
    assert br == {1: ("chrR", 14999, 15700, False)}


def _tx_paf_fields(ts, te, cg, strand="+", qlen=300, qs=0, qe=300,
                   tx="tx"):
    return [tx, str(qlen), str(qs), str(qe), strand, "chrR", "30000",
            str(ts), str(te), str(qe - qs), str(te - ts), "60",
            "tp:A:P", f"cg:Z:{cg}"]


def test_build_gene_placement_plus(tmp_path):
    tx = "A" * 100 + "C" * 200
    p = _tx_paf_fields(1000, 1350, "100M50N200M")
    blocks = build_gene_placement(tx, p, 900, 2000)
    # unclipped intron gets the canonical splice consensus on its edges
    assert blocks == ["A" * 100, "GT" + "N" * 46 + "AG", "C" * 200]


def test_build_gene_placement_max_spacer():
    tx = "A" * 100 + "C" * 200
    p = _tx_paf_fields(1000, 1350, "100M50N200M")
    # unclipped intron longer than max_spacer is truncated, splice
    # consensus edges are kept
    blocks = build_gene_placement(tx, p, 900, 2000, max_spacer=20)
    assert blocks == ["A" * 100, "GT" + "N" * 16 + "AG", "C" * 200]
    # default cap leaves the 50 bp intron untouched
    assert build_gene_placement(tx, p, 900, 2000)[1] == "GT" + "N" * 46 + "AG"


def test_build_gene_placement_clipped_to_bracket():
    tx = "A" * 100 + "C" * 200
    p = _tx_paf_fields(1000, 1350, "100M50N200M")
    # bracket cuts 49 bp off exon1's head and 150 bp off exon2's tail;
    # the intron itself is unclipped -> splice edges stamped
    blocks = build_gene_placement(tx, p, 1049, 1200)
    assert blocks == ["A" * 51, "GT" + "N" * 46 + "AG", "C" * 50]
    # bracket cuts THROUGH the intron -> plain estimated-N spacer tuple
    blocks = build_gene_placement(tx, p, 1049, 1130)
    assert blocks == ["A" * 51, (None, 30)]
    # no exonic base inside the bracket -> empty
    assert build_gene_placement(tx, p, 500, 900) == []


def test_build_gene_placement_small_intron_merged():
    tx = "A" * 100 + "C" * 100
    p = _tx_paf_fields(1000, 1210, "100M10N100M", qlen=200, qe=200)
    # 10 bp intron < min_spacer -> exons stay joined
    blocks = build_gene_placement(tx, p, 900, 2000)
    assert blocks == ["A" * 100 + "C" * 100]


def test_build_gene_placement_minus_strand():
    # '-' CIGAR walks the query backward from qe; chunks are revcomped
    # individually and CIGAR order is already ref-forward.
    e1 = "ACGT" * 25                      # 100 bp
    e2 = "TGCA" * 50                      # 200 bp
    tx = _revcomp(e1 + e2)                # revcomp: e2' then e1'
    p = _tx_paf_fields(1000, 1350, "100M50N200M", strand="-")
    blocks = build_gene_placement(tx, p, 900, 2000)
    # minus-strand gene: intron placeholder carries CT..AC (revcomp
    # of GT..AG)
    assert blocks == [e1, "CT" + "N" * 46 + "AC", e2]


def test_assign_tx_to_genes(tmp_path):
    ref_genes = {"chrR": [(1000, 1400, "g1", "+"), (5000, 5400, "g2", "+")]}
    lines = [
        # full overlap, primary -> assigned to g1
        "txA\t400\t0\t400\t+\tchrR\t30000\t1000\t1400\t400\t400\t60"
        "\ttp:A:P\tcg:Z:400M",
        # secondary -> ignored
        "txB\t400\t0\t400\t+\tchrR\t30000\t1000\t1400\t400\t400\t60"
        "\ttp:A:S\tcg:Z:400M",
        # partial overlap (200/400 = 0.5, passes) but loses g1 to txA
        "txC\t400\t0\t400\t+\tchrR\t30000\t1200\t1600\t400\t400\t60"
        "\ttp:A:P\tcg:Z:400M",
        # overlap too small (<0.5 of both gene and alignment) -> none
        "txD\t600\t0\t600\t+\tchrR\t30000\t5300\t5900\t600\t600\t60"
        "\ttp:A:P\tcg:Z:600M",
    ]
    paf = _paf(tmp_path, lines)
    got = assign_tx_to_genes(paf, ref_genes)
    assert set(got) == {"g1"}
    assert got["g1"][0] == "txA"


def test_plan_placements_multi_gene():
    brackets = {1: ("chrR", 1000, 9000, False)}
    ref_genes = {"chrR": [(2000, 2400, "g1", "+"), (6000, 6500, "g2", "+"),
                          (7000, 7100, "g3", "+")]}  # g3: no transcript
    gene_tx = {"g1": ("tx1", _tx_paf_fields(2000, 2400, "400M",
                                            qlen=400, qe=400, tx="tx1")),
               "g2": ("tx2", _tx_paf_fields(6000, 6500, "500M",
                                            qlen=500, qe=500, tx="tx2"))}
    tx_seqs = {"tx1": "A" * 400, "tx2": "C" * 500}
    pl, detail = plan_placements(brackets, ref_genes, gene_tx, tx_seqs)
    assert pl[1] == [(None, 1000), "A" * 400, (None, 3600),
                     "C" * 500, (None, 2500)]
    assert detail[1] == [("g1", "tx1", 400), ("g2", "tx2", 500)]


def test_plan_placements_flip():
    brackets = {1: ("chrR", 1000, 9000, True)}
    ref_genes = {"chrR": [(2000, 2400, "g1", "+"), (6000, 6500, "g2", "+")]}
    gene_tx = {"g1": ("tx1", _tx_paf_fields(2000, 2400, "400M",
                                            qlen=400, qe=400, tx="tx1")),
               "g2": ("tx2", _tx_paf_fields(6000, 6500, "500M",
                                            qlen=500, qe=500, tx="tx2"))}
    tx_seqs = {"tx1": "A" * 400, "tx2": "C" * 500}
    pl, _d = plan_placements(brackets, ref_genes, gene_tx, tx_seqs)
    # scaffold-forward = reversed ref order, sequences revcomped
    assert pl[1] == [(None, 2500), "G" * 500, (None, 3600),
                     "T" * 400, (None, 1000)]


def test_plan_placements_max_genes():
    brackets = {1: ("chrR", 1000, 9000, False)}
    ref_genes = {"chrR": [(2000, 2400, "g1", "+"), (3000, 3400, "g2", "+"),
                          (4000, 4400, "g3", "+")]}
    gene_tx = {g: (f"tx{g[-1]}", _tx_paf_fields(b, e, "400M",
                                                qlen=400, qe=400,
                                                tx=f"tx{g[-1]}"))
               for g, (b, e) in
               {"g1": (2000, 2400), "g2": (3000, 3400),
                "g3": (4000, 4400)}.items()}
    tx_seqs = {f"tx{i}": "A" * 400 for i in "123"}
    pl, detail = plan_placements(brackets, ref_genes, gene_tx, tx_seqs,
                                 max_genes=2)
    assert pl == {} and detail == {}


def test_write_outputs_placements(tmp_path):
    """Ref-guided placement: gap replaced by exon components + estimated-N
    spacer rows; placements win over extensions; parts stay consecutive."""
    contigs = {"c1": "AAAAACGT", "c2": "CCCCGGGG"}
    agp_text = "\n".join([
        "s1\t1\t8\t1\tW\tc1\t1\t8\t+",
        "s1\t9\t108\t2\tU\t100\tscaffold\tyes\tna",
        "s1\t109\t116\t3\tW\tc2\t1\t8\t+",
    ]) + "\n"
    agp_lines = AGPReader().parse(agp_text)
    out_agp = str(tmp_path / "o.agp")
    out_fa = str(tmp_path / "o.fa")
    report = write_outputs(
        agp_lines, {1}, {}, contigs, out_agp, out_fa,
        extensions={1: ("LLL", "RRR")},          # must lose to placements
        placements={1: ["AAA", (None, 50), "CCC"]})
    assert report["gaps_placed"] == 1
    assert report["gaps_extended"] == 0
    assert report["bases_filled"] == 6
    lines = AGPReader().parse(open(out_agp).read())
    assert [l.part_number for l in lines] == list(range(1, len(lines) + 1))
    w = [l for l in lines if isinstance(l, AGPSeqLine)]
    g = [l for l in lines if isinstance(l, AGPGapLine)]
    assert len(w) == 4 and len(g) == 1      # c1, AAA, gap, CCC, c2
    assert g[0].gap_length == 50 and g[0].linkage == "no"
    assert (w[1].object_beg, w[1].object_end) == (9, 11)
    assert (g[0].object_beg, g[0].object_end) == (12, 61)
    assert (w[2].object_beg, w[2].object_end) == (62, 64)
    assert (w[3].object_beg, w[3].object_end) == (65, 72)
    from nearscaff.gapfill import read_fasta as _rf
    fa = _rf(out_fa)
    assert fa == {"s1": "AAAAACGT" + "AAA" + "N" * 50 + "CCC" + "CCCCGGGG"}


@pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
def test_run_rnafill_ref_guided_end_to_end(tmp_path):
    """Synthetic locus: a gene whose region is entirely a gap in the query
    gets placed from ref coordinates, exons written, intron as N."""
    rng = random.Random(7)
    ref = "".join(rng.choice("ACGT") for _ in range(30000))
    # gene g1 on chrR: exons [15000,15200) + [15350,15700), intron 150 bp
    e1, e2 = ref[15000:15200], ref[15350:15700]
    tx = e1 + e2
    # query: cA = ref[10000:15000), cB = ref[15700:20700), gap in between
    query_fa = _write(tmp_path, "q.fa",
                      f">cA\n{ref[10000:15000]}\n>cB\n{ref[15700:20700]}\n")
    ref_fa = _write(tmp_path, "ref.fa", f">chrR\n{ref}\n")
    tx_fa = _write(tmp_path, "tx.fa", f">txG\n{tx}\n")
    gff = _write(tmp_path, "ref.gff",
                 "chrR\tsrc\tmRNA\t15001\t15700\t.\t+\t.\tID=g1\n")
    agp = _write(tmp_path, "in.agp",
                 "s1\t1\t5000\t1\tW\tcA\t1\t5000\t+\n"
                 "s1\t5001\t5100\t2\tU\t100\tscaffold\tyes\tna\n"
                 "s1\t5101\t10100\t3\tW\tcB\t1\t5000\t+\n")
    tiered = _tiered_paf(tmp_path, [
        f"cA\t5000\t0\t5000\t+\tchrR\t30000\t10000\t15000\t5000\t5000\t60"
        "\ttp:A:P\tnc:Z:protein",
        f"cB\t5000\t0\t5000\t+\tchrR\t30000\t15700\t20700\t5000\t5000\t60"
        "\ttp:A:P\tnc:Z:protein",
    ])
    report = run_rnafill(agp, query_fa, tiered, str(tmp_path / "out"),
                         transcripts=[tx_fa], ref=ref_fa, ref_gff=gff,
                         overlap_closure=False)
    assert report["genes_placed"] == 1
    assert report["gaps_placed"] == 1
    from nearscaff.gapfill import read_fasta as _rf
    fa = _rf(str(tmp_path / "out" / "nearscaff.rnafill.scaffolds.fa"))
    scaf = fa["s1"]
    # exons (allowing a few bp of alignment fuzz at the edges) are in,
    # in scaffold-forward order, separated by an N run
    i1 = scaf.find(e1[5:-5])
    i2 = scaf.find(e2[5:-5])
    assert 5000 < i1 < i2
    between = scaf[i1 + len(e1[5:-5]):i2]
    assert set(between) <= {"N"} or "N" * 50 in between
    # closures manifest records the placement
    manifest = open(tmp_path / "out" / "rnafill.closures.tsv").read()
    assert "placed:g1" in manifest and "txG" in manifest


def test_find_intact_proteins(tmp_path):
    """Intact = aligned coverage >= min_cov AND translation without X."""
    from nearscaff.rnafill import find_intact_proteins
    trans = _write(tmp_path, "t.faa", "\n".join([
        # intact: cov 0.95, no X
        "P1\t100\t0\t95\t+\ts1\t1000\t0\t300\t95\t95\t60",
        "##STA\t" + "M" * 95,
        # X-containing: full cov but crosses an N run
        "P2\t100\t0\t100\t+\ts1\t1000\t400\t700\t100\t100\t60",
        "##STA\t" + "M" * 50 + "XX" + "M" * 48,
        # truncated: cov 0.5, no X
        "P3\t100\t0\t50\t+\ts1\t1000\t800\t950\t50\t50\t60",
        "##STA\t" + "M" * 50,
        # two hits: one truncated with X, one intact -> intact
        "P4\t100\t0\t60\t+\ts1\t1000\t0\t180\t60\t60\t60",
        "##STA\t" + "M" * 58 + "XX",
        "P4\t100\t5\t100\t-\ts2\t1000\t0\t285\t95\t95\t60",
        "##STA\t" + "M" * 95,
    ]) + "\n")
    assert find_intact_proteins(trans) == {"P1", "P4"}
    assert find_intact_proteins(trans, min_cov=0.4) == {"P1", "P3", "P4"}
