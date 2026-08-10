"""Unit tests for nearscaff chromosome-purity enforcement."""
import networkx as nx
import pytest

from nearscaff.config import NearscaffConfig
from nearscaff.pipeline import _enforce_chromosome_purity


def test_scaffold_config_chr_purity_defaults():
    cfg = NearscaffConfig()
    assert cfg.scaffold.enforce_chr_purity is True
    assert cfg.scaffold.min_chr_share == 0.20
    assert cfg.scaffold.min_chr_len == 1_000_000


def _cover(edges):
    """Build a cover graph. Each edge = (u, v, weight); inf = intra-contig."""
    g = nx.Graph()
    for u, v, w in edges:
        g.add_edge(u, v, weight=w)
    return g


def test_impure_two_chromosome_component_is_split():
    cover = _cover([
        ("A_b", "A_e", float("inf")),
        ("A_e", "B_b", 0.5),     # cross-chromosome matching edge
        ("B_b", "B_e", float("inf")),
    ])
    contig_ref = {"A": ("chr1", 1, 2_000_000), "B": ("chr2", 1, 2_000_000)}
    lengths = {"A": 2_000_000, "B": 2_000_000}

    n_cut = _enforce_chromosome_purity(cover, contig_ref, lengths,
                                       min_share=0.20, min_len=1_000_000)

    assert n_cut == 1
    comps = [set(c) for c in nx.connected_components(cover)]
    assert len(comps) == 2
    assert {"A_b", "A_e"} in comps
    assert {"B_b", "B_e"} in comps


def test_pure_component_is_left_untouched():
    cover = _cover([
        ("A_b", "A_e", float("inf")),
        ("A_e", "B_b", 0.5),
        ("B_b", "B_e", float("inf")),
    ])
    contig_ref = {"A": ("chr1", 1, 2_000_000),
                  "B": ("chr1", 5_000_000, 7_000_000)}
    lengths = {"A": 2_000_000, "B": 2_000_000}

    n_cut = _enforce_chromosome_purity(cover, contig_ref, lengths,
                                       min_share=0.20, min_len=1_000_000)
    assert n_cut == 0
    assert nx.number_connected_components(cover) == 1


def test_minor_contig_absorbed_into_dominant_neighbour():
    # A(chr1,big) --strong-- D(chr3,<1Mb minor) --weak-- B(chr2,big)
    cover = _cover([
        ("A_b", "A_e", float("inf")),
        ("A_e", "D_b", 0.9),
        ("D_b", "D_e", float("inf")),
        ("D_e", "B_b", 0.1),
        ("B_b", "B_e", float("inf")),
    ])
    contig_ref = {"A": ("chr1", 1, 2_000_000),
                  "D": ("chr3", 1, 500_000),
                  "B": ("chr2", 1, 2_000_000)}
    lengths = {"A": 2_000_000, "D": 500_000, "B": 2_000_000}

    n_cut = _enforce_chromosome_purity(cover, contig_ref, lengths,
                                       min_share=0.20, min_len=1_000_000)
    assert n_cut == 1
    comps = [set(c) for c in nx.connected_components(cover)]
    assert len(comps) == 2
    chr1_comp = next(c for c in comps if "A_b" in c)
    assert chr1_comp == {"A_b", "A_e", "D_b", "D_e"}  # D absorbed into chr1
    assert {"B_b", "B_e"} in comps


def test_unknown_contig_bridge_is_severed():
    # A(chr1) --strong-- U(no ref) --weak-- B(chr2)
    cover = _cover([
        ("A_b", "A_e", float("inf")),
        ("A_e", "U_b", 0.9),
        ("U_b", "U_e", float("inf")),
        ("U_e", "B_b", 0.1),
        ("B_b", "B_e", float("inf")),
    ])
    contig_ref = {"A": ("chr1", 1, 2_000_000),
                  "B": ("chr2", 1, 2_000_000)}  # U absent -> unknown
    lengths = {"A": 2_000_000, "U": 1_000_000, "B": 2_000_000}

    n_cut = _enforce_chromosome_purity(cover, contig_ref, lengths,
                                       min_share=0.20, min_len=1_000_000)
    assert n_cut == 1
    comps = [set(c) for c in nx.connected_components(cover)]
    assert len(comps) == 2
    assert {"B_b", "B_e"} in comps  # B isolated; U stays with A


def test_single_significant_chromosome_not_split():
    # A(chr1,big) - M(chr2,<1Mb minor): only chr1 is significant -> skip.
    cover = _cover([
        ("A_b", "A_e", float("inf")),
        ("A_e", "M_b", 0.5),
        ("M_b", "M_e", float("inf")),
    ])
    contig_ref = {"A": ("chr1", 1, 5_000_000), "M": ("chr2", 1, 100_000)}
    lengths = {"A": 5_000_000, "M": 100_000}

    n_cut = _enforce_chromosome_purity(cover, contig_ref, lengths,
                                       min_share=0.20, min_len=1_000_000)
    assert n_cut == 0
    assert nx.number_connected_components(cover) == 1
