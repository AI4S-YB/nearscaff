"""Unit tests for nearscaff synteny-consistent contig reordering."""
import os
import tempfile

import networkx as nx

from nearscaff.config import NearscaffConfig
from nearscaff.pipeline import _write_anchors_tsv, _load_anchor_coords, _extract_agp_paths
from nearscaff.agp import AGPSeqLine


class _A:
    """Minimal anchor stand-in matching the Anchor attribute names."""
    def __init__(self, qc, qg, rc, rg, rs, re_, strand, score, ident):
        self.query_contig = qc
        self.query_gene = qg
        self.ref_chr = rc
        self.ref_gene = rg
        self.r_start = rs
        self.r_end = re_
        self.strand = strand
        self.score = score
        self.identity = ident


def test_write_then_load_anchor_coords_roundtrip():
    anchors = [
        _A("ctg1", "g1", "chr1", "RG1", 1000, 2000, "+", 500, 0.9),
        _A("ctg1", "g2", "chr1", "RG2", 5000, 6000, "+", 400, 0.8),
        _A("ctg2", "g3", "chr2", "RG3", 300, 400, "-", 300, 0.7),
    ]
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "gene_anchors.tsv")
        _write_anchors_tsv(anchors, path)
        coords = _load_anchor_coords(path)
    assert sorted(coords.keys()) == ["ctg1", "ctg2"]
    assert {c for (c, rs, st) in coords["ctg1"]} == {"chr1"}
    assert sorted(rs for (c, rs, st) in coords["ctg1"]) == [1000, 5000]
    assert coords["ctg2"][0] == ("chr2", 300, "-")


def test_load_anchor_coords_degrades_on_old_format():
    # old 7-column format (no ref_start/ref_end) -> empty dict (caller falls back)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "old.tsv")
        with open(path, "w") as fh:
            fh.write("query_contig\tquery_gene\tref_chr\tref_gene\tstrand\tscore\tidentity\n")
            fh.write("ctg1\tg1\tchr1\tRG1\t+\t500\t0.9\n")
        assert _load_anchor_coords(path) == {}


def test_scaffold_config_synteny_reorder_defaults():
    cfg = NearscaffConfig()
    assert cfg.scaffold.synteny_reorder is True
    assert cfg.scaffold.suspect_anchor_span == 5_000_000
    assert cfg.scaffold.suspect_divergence == 10_000_000


def _cover_path(bases):
    """Cover graph: a single path through these contig bases."""
    g = nx.Graph()
    for b in bases:
        g.add_edge(b + "_b", b + "_e", weight=float("inf"))
    for i in range(len(bases) - 1):
        g.add_edge(bases[i] + "_e", bases[i + 1] + "_b", weight=0.5)
    return g


def _order_of(agp_lines):
    return [l.component_id for l in agp_lines if isinstance(l, AGPSeqLine)]


def test_synteny_reorder_uses_protein_coordinate_over_noisy_midpoint():
    # protein coords: A=100, B=200, C=300 ; nucleotide midpoints (noisy): A=200,B=100,C=300
    cover = _cover_path(["B", "A", "C"])  # graph/pre-sort order ≠ protein order → relocated>0
    lengths = {"A": 1000, "B": 1000, "C": 1000}
    contig_ref = {"A": ("chr1", 150, 250), "B": ("chr1", 50, 150), "C": ("chr1", 250, 350)}
    anchor_coords = {
        "A": [("chr1", 100, "+")], "B": [("chr1", 200, "+")], "C": [("chr1", 300, "+")],
    }
    report = []
    lines = _extract_agp_paths(
        cover, lengths, contig_ref=contig_ref,
        anchor_coords=anchor_coords, synteny_reorder=True,
        suspect_anchor_span=5_000_000, suspect_divergence=10_000_000,
        report=report,
    )
    assert _order_of(lines) == ["A", "B", "C"]  # protein order, not noisy midpoint
    assert report[0][3] >= 2  # n_relocated: midpoint would have been B,A,C


def test_synteny_reorder_falls_back_to_midpoint_without_anchors():
    cover = _cover_path(["A", "B", "C"])
    lengths = {"A": 1000, "B": 1000, "C": 1000}
    contig_ref = {"A": ("chr1", 50, 150), "B": ("chr1", 150, 250), "C": ("chr1", 250, 350)}
    lines = _extract_agp_paths(cover, lengths, contig_ref=contig_ref,
                               anchor_coords={}, synteny_reorder=True, report=[])
    assert _order_of(lines) == ["A", "B", "C"]  # midpoint order


def test_synteny_reorder_flag_suspect_for_huge_anchor_span():
    cover = _cover_path(["A", "X", "B"])
    lengths = {"A": 1000, "X": 1000, "B": 1000}
    contig_ref = {"A": ("chr1", 1, 100), "X": ("chr1", 100, 200), "B": ("chr1", 200, 300)}
    anchor_coords = {"A": [("chr1", 100, "+")],
                     "X": [("chr1", 100, "+"), ("chr1", 8_000_100, "+")],  # span 8Mb
                     "B": [("chr1", 9000, "+")]}
    report = []
    _extract_agp_paths(cover, lengths, contig_ref=contig_ref,
                       anchor_coords=anchor_coords, synteny_reorder=True,
                       suspect_anchor_span=5_000_000, suspect_divergence=10_000_000,
                       report=report)
    assert report[0][4] == 1  # n_suspect: X flagged


def test_synteny_reorder_disabled_matches_midpoint():
    cover = _cover_path(["A", "B", "C"])
    lengths = {"A": 1000, "B": 1000, "C": 1000}
    # midpoints A<B<C, but anchors would reorder
    contig_ref = {"A": ("chr1", 50, 150), "B": ("chr1", 150, 250), "C": ("chr1", 250, 350)}
    anchor_coords = {"A": [("chr1", 300, "+")], "B": [("chr1", 200, "+")], "C": [("chr1", 100, "+")]}
    lines = _extract_agp_paths(cover, lengths, contig_ref=contig_ref,
                               anchor_coords=anchor_coords, synteny_reorder=False, report=[])
    assert _order_of(lines) == ["A", "B", "C"]  # midpoint order (disabled)


def test_synteny_reorder_uses_anchor_chr_not_contig_ref_chr():
    # Production scenario: anchors are in a DIFFERENT coordinate space
    # (ref-protein chr "KAH7") than contig_ref (scaffolding chr "CM040171").
    # The reorder must key off the ANCHOR chr (dominant = KAH7), not the
    # contig_ref chr, or it silently no-ops (the original bug).
    cover = _cover_path(["B", "A", "C"])  # graph order != protein order
    lengths = {"A": 1000, "B": 1000, "C": 1000}
    contig_ref = {"A": ("CM040171", 50, 150), "B": ("CM040171", 150, 250),
                  "C": ("CM040171", 250, 350)}
    anchor_coords = {"A": [("KAH7", 100, "+")], "B": [("KAH7", 200, "+")],
                     "C": [("KAH7", 300, "+")]}  # KAH7 protein order A<B<C
    report = []
    lines = _extract_agp_paths(cover, lengths, contig_ref=contig_ref,
                               anchor_coords=anchor_coords, synteny_reorder=True,
                               suspect_anchor_span=5_000_000, suspect_divergence=10_000_000,
                               report=report)
    assert _order_of(lines) == ["A", "B", "C"]  # ordered by KAH7 protein coords
    assert report[0][2] == 3  # n_with_anchors: all 3 matched (anchor chr KAH7)


def test_synteny_reorder_ignores_foreign_chromosome_anchors():
    # scaffold majority chr = chr1 (A,B,C all ref chr1). X has anchors ONLY on chr2.
    # X must fall back to its chr1 nucleotide midpoint, not use chr2 protein coord.
    cover = _cover_path(["A", "X", "B"])
    lengths = {"A": 1000, "X": 1000, "B": 1000}
    contig_ref = {"A": ("chr1", 1, 100), "X": ("chr1", 150, 250), "B": ("chr1", 300, 400)}
    anchor_coords = {"A": [("chr1", 50, "+")],
                     "X": [("chr2", 9_999_999, "+")],   # foreign chr -> ignored
                     "B": [("chr1", 350, "+")]}
    report = []
    lines = _extract_agp_paths(cover, lengths, contig_ref=contig_ref,
                               anchor_coords=anchor_coords, synteny_reorder=True,
                               suspect_anchor_span=5_000_000, suspect_divergence=10_000_000,
                               report=report)
    # X placed by its chr1 midpoint (150-250 -> 200), between A(50) and B(350)
    assert _order_of(lines) == ["A", "X", "B"]
    # X did NOT count as having anchors (its only anchor is foreign-chr)
    assert report[0][2] == 2  # n_with_anchors = A and B only


def test_build_config_synteny_reorder_defaults():
    from argparse import Namespace
    from nearscaff.cli import _build_config
    cfg = _build_config(Namespace(threads=4, output="o"))
    assert cfg.scaffold.synteny_reorder is True
    assert cfg.scaffold.suspect_anchor_span == 5_000_000
    assert cfg.scaffold.suspect_divergence == 10_000_000


def test_build_config_no_synteny_reorder():
    from argparse import Namespace
    from nearscaff.cli import _build_config
    cfg = _build_config(Namespace(threads=4, output="o",
                                  no_synteny_reorder=True,
                                  suspect_anchor_span=3_000_000,
                                  suspect_divergence=7_000_000))
    assert cfg.scaffold.synteny_reorder is False
    assert cfg.scaffold.suspect_anchor_span == 3_000_000
    assert cfg.scaffold.suspect_divergence == 7_000_000
