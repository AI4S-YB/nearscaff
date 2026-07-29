"""Unit tests for nearscaff.pipeline._extract_agp_paths."""
import networkx as nx

from nearscaff.pipeline import _extract_agp_paths, MIN_ORIENT_MAPQ
from nearscaff.agp import AGPSeqLine


def _two_contig_cover():
    """Cover graph with one path: ctg1_b - ctg1_e - ctg2_b - ctg2_e."""
    g = nx.Graph()
    g.add_edge("ctg1_b", "ctg1_e")
    g.add_edge("ctg1_e", "ctg2_b")
    g.add_edge("ctg2_b", "ctg2_e")
    return g


def _seq_lines(lines):
    return [l for l in lines if isinstance(l, AGPSeqLine)]


def test_components_reordered_by_reference_midpoint():
    cover = _two_contig_cover()
    lengths = {"ctg1": 100, "ctg2": 100}
    # Graph order is ctg1 -> ctg2, but reference midpoints say the opposite.
    contig_ref = {"ctg1": ("chr1", 9000, 9100), "ctg2": ("chr1", 1000, 1100)}
    lines = _extract_agp_paths(cover, lengths, contig_ref=contig_ref)
    order = [l.component_id for l in _seq_lines(lines)]
    assert order == ["ctg2", "ctg1"]


def test_orientation_marked_unknown_on_low_mapq():
    cover = _two_contig_cover()
    lengths = {"ctg1": 100, "ctg2": 100}
    contig_strand = {"ctg1": "-", "ctg2": "-"}
    contig_mapq = {"ctg1": 60, "ctg2": MIN_ORIENT_MAPQ - 1}
    lines = _extract_agp_paths(cover, lengths,
                               contig_strand=contig_strand,
                               contig_mapq=contig_mapq)
    orient = {l.component_id: l.orientation for l in _seq_lines(lines)}
    assert orient["ctg1"] == "-"
    assert orient["ctg2"] == "?"


def test_orientation_from_strand_when_mapq_sufficient():
    cover = _two_contig_cover()
    lengths = {"ctg1": 100, "ctg2": 100}
    contig_strand = {"ctg1": "+", "ctg2": "-"}
    contig_mapq = {"ctg1": MIN_ORIENT_MAPQ, "ctg2": 60}
    lines = _extract_agp_paths(cover, lengths,
                               contig_strand=contig_strand,
                               contig_mapq=contig_mapq)
    orient = {l.component_id: l.orientation for l in _seq_lines(lines)}
    assert orient == {"ctg1": "+", "ctg2": "-"}
