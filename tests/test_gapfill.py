"""Unit tests for the shared gapfill machinery (tiers, classify, write-back)."""
import os

from nearscaff.agp import AGPSeqLine, AGPGapLine, AGPReader
from nearscaff.gapfill import (read_tiers, classify_gaps, write_outputs,
                               read_fasta, _revcomp)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# c1(protein) gap c2(asm5) gap c3(asm20) gap c4(asm5,-)
CONTIGS = {"c1": "AAAAACGT", "c2": "CCCCGGGG", "c3": "TTTTAAAA",
           "c4": "AAACCCGG"}  # RC(c4) = CCGGGTTT

AGP_TEXT = "\n".join([
    "s1\t1\t8\t1\tW\tc1\t1\t8\t+",
    "s1\t9\t108\t2\tU\t100\tscaffold\tyes\tna",
    "s1\t109\t116\t3\tW\tc2\t1\t8\t+",
    "s1\t117\t216\t4\tU\t100\tscaffold\tyes\tna",
    "s1\t217\t224\t5\tW\tc3\t1\t8\t+",
    "s1\t225\t324\t6\tU\t100\tscaffold\tyes\tna",
    "s1\t325\t332\t7\tW\tc4\t1\t8\t-",
]) + "\n"

TIERED_PAF = "\n".join([
    "c1\t8\t0\t8\t+\tChr1\t100\t0\t8\t8\t8\t60\tnc:Z:protein",
    "c2\t8\t0\t8\t+\tChr1\t100\t10\t18\t8\t8\t60\tnc:Z:asm5",
    "c3\t8\t0\t8\t+\tChr1\t100\t20\t28\t8\t8\t60\tnc:Z:asm20",
    "c4\t8\t0\t8\t-\tChr1\t100\t30\t38\t8\t8\t60\tnc:Z:asm5",
]) + "\n"


def _agp_lines():
    return AGPReader().parse(AGP_TEXT)


def _tiers(tmp_path):
    return read_tiers(_write(tmp_path, "tiered.paf", TIERED_PAF))


def test_read_tiers(tmp_path):
    tiers = _tiers(tmp_path)
    assert tiers == {"c1": "protein", "c2": "asm5", "c3": "asm20", "c4": "asm5"}


def test_classify_gaps_requires_both_flanks_allowed(tmp_path):
    tiers = _tiers(tmp_path)
    eligible = classify_gaps(_agp_lines(), tiers)
    # line indices: 1 (c1|c2) eligible; 3 (c2|c3) and 5 (c3|c4) touch asm20
    assert eligible == {1}


def test_classify_gaps_custom_tiers(tmp_path):
    tiers = _tiers(tmp_path)
    eligible = classify_gaps(_agp_lines(), tiers,
                             allowed={"protein", "asm5", "asm20"})
    assert eligible == {1, 3, 5}


def test_write_outputs_end_to_end(tmp_path):
    agp_lines = _agp_lines()
    eligible = {1}
    fills = {1: "ACGTAC"}
    out_agp = str(tmp_path / "out.agp")
    out_fa = str(tmp_path / "out.fa")
    report = write_outputs(agp_lines, eligible, fills, CONTIGS,
                           out_agp, out_fa)
    assert report == {"gaps_total": 3, "gaps_eligible": 1,
                      "gaps_closed": 1, "gaps_resized": 0,
                      "gaps_endjoined": 0, "gaps_extended": 0,
                      "gaps_placed": 0, "bases_filled": 6}

    lines = AGPReader().parse(open(out_agp).read())
    # gap 1 replaced by a W component of length 6; coordinates re-numbered
    w = [l for l in lines if isinstance(l, AGPSeqLine)]
    g = [l for l in lines if isinstance(l, AGPGapLine)]
    assert len(w) == 5 and len(g) == 2
    assert w[1].component_id == "s1_gapfill2"
    assert (w[1].object_beg, w[1].object_end) == (9, 14)
    assert (w[2].object_beg, w[2].object_end) == (15, 22)
    # part numbers are consecutive
    assert [l.part_number for l in lines] == list(range(1, len(lines) + 1))

    fa = read_fasta(out_fa)
    assert fa == {"s1": "AAAAACGT" + "ACGTAC" + "CCCCGGGG" + "N" * 100
                       + "TTTTAAAA" + "N" * 100 + "CCGGGTTT"}


def test_write_outputs_collapsed_gap(tmp_path):
    agp_lines = _agp_lines()
    fills = {1: ""}  # collapsed: gap line dropped
    out_agp = str(tmp_path / "out.agp")
    out_fa = str(tmp_path / "out.fa")
    report = write_outputs(agp_lines, {1}, fills, CONTIGS, out_agp, out_fa)
    assert report["gaps_closed"] == 1
    assert report["bases_filled"] == 0
    lines = AGPReader().parse(open(out_agp).read())
    # first gap gone; coordinates/part numbers stay consecutive
    assert [l.part_number for l in lines] == list(range(1, len(lines) + 1))
    assert lines[0].object_end == 8
    assert (lines[1].object_beg, lines[1].object_end) == (9, 16)
    fa = read_fasta(out_fa)
    assert fa["s1"].startswith("AAAAACGT" + "CCCCGGGG")


def test_write_outputs_endjoin_trim(tmp_path):
    """End-join: gap dropped AND the preceding component shortened."""
    agp_lines = _agp_lines()
    out_agp = str(tmp_path / "out.agp")
    out_fa = str(tmp_path / "out.fa")
    report = write_outputs(agp_lines, {1}, {}, CONTIGS, out_agp, out_fa,
                           trims={1: 3})  # trim 3 bp off c1's 3' end
    assert report["gaps_endjoined"] == 1
    lines = AGPReader().parse(open(out_agp).read())
    # first gap gone; c1 shortened from 8 to 5 bp
    assert lines[0].component_id == "c1"
    assert lines[0].component_end == 5
    assert (lines[0].object_beg, lines[0].object_end) == (1, 5)
    assert (lines[1].object_beg, lines[1].object_end) == (6, 13)
    fa = read_fasta(out_fa)
    assert fa["s1"].startswith("AAAAA" + "CCCCGGGG")


def test_write_outputs_endjoin_trim_minus_strand(tmp_path):
    """End-join trim of a "-" component must cut the contig's 5' side
    (component_beg), which is its scaffold-side 3' end."""
    # s2: cA(10 bp, '-') - gap - cB(8 bp, '+'); trim 4 bp at the join
    contigs = {"cA": "AAAACCCCGG", "cB": "CCCCGGGG"}
    agp_text = "\n".join([
        "s2\t1\t10\t1\tW\tcA\t1\t10\t-",
        "s2\t11\t15\t2\tU\t5\tscaffold\tyes\tna",
        "s2\t16\t23\t3\tW\tcB\t1\t8\t+",
    ]) + "\n"
    agp_lines = AGPReader().parse(agp_text)
    out_agp = str(tmp_path / "out.agp")
    out_fa = str(tmp_path / "out.fa")
    report = write_outputs(agp_lines, set(), {}, contigs, out_agp, out_fa,
                           trims={1: 4})
    assert report["gaps_endjoined"] == 1
    lines = AGPReader().parse(open(out_agp).read())
    # cA keeps its 3' end in contig coords: beg 1+4=5, end unchanged
    assert lines[0].component_id == "cA"
    assert (lines[0].component_beg, lines[0].component_end) == (5, 10)
    assert (lines[0].object_beg, lines[0].object_end) == (1, 6)
    fa = read_fasta(out_fa)
    # scaffold keeps rc(cA[4:10]) = the flank-adjacent end
    assert fa["s2"] == _revcomp("CCCCGG") + "CCCCGGGG"
