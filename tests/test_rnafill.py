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
