"""Unit tests for homeolog.py — collinear homeolog pair discovery."""

import os

import pytest

from nearscaff.homeolog import (
    _read_bed, cluster_side_blocks, find_homeolog_pairs, gene_level_pep,
)


@pytest.fixture()
def beds(tmp_path):
    sp4_bed = tmp_path / "sp4.bed"
    sp2_bed = tmp_path / "sp2.bed"
    # sp4: two subgenome copies of 6 loci on two scaffolds each
    sp4_lines = []
    for i in range(1, 7):
        sp4_lines.append(f"ctgA\tMP{i:06d}\t{i * 1000}\t{i * 1000 + 300}")
        sp4_lines.append(f"ctgB\tMP{i:06d}x\t{i * 1000}\t{i * 1000 + 300}")
    sp4_bed.write_text("\n".join(sp4_lines) + "\n")
    sp2_bed.write_text(
        "\n".join(f"chr1\tDP{i:06d}\t{i * 1000}\t{i * 1000 + 300}"
                  for i in range(1, 7)) + "\n")
    return str(sp4_bed), str(sp2_bed)


def _write_anchors(path, rows, block=0):
    """rows: list of (gene1, gene2) — orientation mixed on purpose."""
    with open(path, "w") as f:
        f.write("###\n")
        for g1, g2 in rows:
            f.write(f"{g1}\t{g2}\t100\n")


class TestPairExtraction:
    def test_orientation_normalized_1v2(self, tmp_path, beds):
        """Anchor rows come in both orientations; classification must use
        bed membership, not column position."""
        sp4_bed, sp2_bed = beds
        work = tmp_path / "work"
        work.mkdir(); (work / "jcvi").mkdir()
        rows = []
        for i in range(1, 7):
            # normal orientation: sp2 first
            rows.append((f"DP{i:06d}", f"MP{i:06d}"))
            # reversed orientation: sp4 first
            rows.append((f"MP{i:06d}x", f"DP{i:06d}"))
        _write_anchors(work / "jcvi" / "sp2.sp4.anchors", rows)
        clusters = find_homeolog_pairs(
            "x", "y", sp4_bed, sp2_bed, str(work), min_block_pairs=5)
        triples = [t for v in clusters.values() for t in v]
        assert len(triples) == 6
        for sp2, sp41, sp42 in triples:
            assert sp2.startswith("DP")
            assert sp41.startswith("MP") and sp42.startswith("MP")
            assert {sp41, sp42} == {f"MP{sp2[2:]}", f"MP{sp2[2:]}x"}

    def test_sp2_with_one_or_three_hits_excluded(self, tmp_path, beds):
        sp4_bed, sp2_bed = beds
        work = tmp_path / "work"
        work.mkdir(); (work / "jcvi").mkdir()
        rows = [(f"DP{i:06d}", f"MP{i:06d}") for i in range(1, 7)]
        rows.append(("DP000001", "MP000002x"))  # DP1 gets a 2nd hit
        rows.append(("DP000001", "MP000003x"))  # and a 3rd
        _write_anchors(work / "jcvi" / "sp2.sp4.anchors", rows)
        clusters = find_homeolog_pairs(
            "x", "y", sp4_bed, sp2_bed, str(work), min_block_pairs=0)
        sp2s = {t[0] for v in clusters.values() for t in v}
        assert "DP000001" not in sp2s  # 3 hits -> excluded

    def test_small_blocks_filtered(self, tmp_path, beds):
        sp4_bed, sp2_bed = beds
        work = tmp_path / "work"
        work.mkdir(); (work / "jcvi").mkdir()
        _write_anchors(work / "jcvi" / "sp2.sp4.anchors",
                       [(f"DP{i:06d}", f"MP{i:06d}") for i in range(1, 7)]
                       + [(f"MP{i:06d}x", f"DP{i:06d}")
                          for i in range(1, 7)])
        # block holds 6 1v2 pairs: kept when threshold is 5 (6 > 5) ...
        clusters = find_homeolog_pairs(
            "x", "y", sp4_bed, sp2_bed, str(work), min_block_pairs=5)
        assert len([t for v in clusters.values() for t in v]) == 6
        # ... filtered when the threshold is 6 (6 > 6 is False)
        clusters2 = find_homeolog_pairs(
            "x", "y", sp4_bed, sp2_bed, str(work), min_block_pairs=6)
        assert clusters2 == {}


class TestClusterSideBlocks:
    def test_spans_and_multiscaffold_drop(self, tmp_path):
        sp4_bed = tmp_path / "sp4.bed"
        sp4_bed.write_text(
            "ctgA\tg1\t100\t200\nctgA\tg2\t300\t400\n"
            "ctgB\th1\t150\t250\nctgB\th2\t350\t450\n"
            "ctgC\tx1\t100\t200\n")
        clusters = {
            "cluster_0": [("d1", "g1", "h1"), ("d2", "g2", "h2")],
            "cluster_1": [("d3", "g1", "x1")],  # B side fine, A side fine
            "cluster_2": [("d4", "g1", "h1"), ("d5", "x1", "h2")],  # A spans 2
        }
        blocks = cluster_side_blocks(clusters, str(sp4_bed))
        assert blocks["cluster_0"]["A"] == ("ctgA", 100, 400)
        assert blocks["cluster_0"]["B"] == ("ctgB", 150, 450)
        assert "cluster_2" not in blocks


class TestGeneLevelPep:
    def test_longest_isoform(self, tmp_path):
        pep = tmp_path / "ref.pep"
        pep.write_text(">g1.mRNA1\nAAAA\n>g1.mRNA2\nAAAAAA\n>g2.mRNA1\nCC\n")
        out = tmp_path / "gene.pep"
        result = gene_level_pep(str(pep), str(out))
        text = open(result).read()
        assert text.count(">") == 2
        assert ">g1.mRNA2\nAAAAAA\n" in text

    def test_gene_level_passthrough(self, tmp_path):
        pep = tmp_path / "ref.pep"
        pep.write_text(">KAF1.1\nAA\n>KAF2.1\nCC\n")
        out = tmp_path / "gene.pep"
        assert gene_level_pep(str(pep), str(out)) == str(pep)
        assert not out.exists()
