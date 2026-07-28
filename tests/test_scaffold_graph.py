"""Tests for Fused Scaffold Graph."""
from nearscaff.scaffold_graph import FusedScaffoldGraph
from nearscaff.types import EdgeType


def test_add_and_fuse_edges():
    sg = FusedScaffoldGraph()
    sg.add_node("ctg1_b"); sg.add_node("ctg1_e")
    sg.add_node("ctg2_b"); sg.add_node("ctg2_e")

    sg.add_fused_edge("ctg1_e", "ctg2_b", weight=0.9,
                      edge_type=EdgeType.PROTEIN_SYNTENY,
                      gap_size=15000, gap_type="scaffold")
    sg.add_fused_edge("ctg1_e", "ctg2_b", weight=0.6,
                      edge_type=EdgeType.NUCLEOTIDE_CHAIN,
                      gap_size=14500, gap_type="scaffold")

    fused = sg.fuse_weights()
    assert ("ctg1_e", "ctg2_b") in fused
    # 0.9*1.0 + 0.6*0.85 = 0.9 + 0.51 = 1.41
    assert abs(fused[("ctg1_e", "ctg2_b")] - 1.41) < 0.01


def test_node_naming_validation():
    sg = FusedScaffoldGraph()
    sg.add_node("ctg1_b")
    sg.add_node("ctg1_e")
    try:
        sg.add_node("ctg1")  # no _b/_e suffix
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_best_buddy_scale():
    sg = FusedScaffoldGraph()
    # Three contigs: ctg1, ctg2, ctg3
    for ctg in ["ctg1", "ctg2", "ctg3"]:
        sg.add_node(ctg + "_b"); sg.add_node(ctg + "_e")

    weights = {
        ("ctg1_e", "ctg2_b"): 10.0,   # strong
        ("ctg2_e", "ctg3_b"): 8.0,    # strong
        ("ctg1_e", "ctg3_b"): 3.0,    # weak alternative
    }
    scaled = sg.best_buddy_scale(weights)
    # Reciprocal best hit: 10 vs next-best-alt 3 → 10/3 ≈ 3.33
    assert abs(scaled[("ctg1_e", "ctg2_b")] - 3.333) < 0.01
    # Non-reciprocal: denominator = max(max_incident_u, max_incident_v) = 10
    assert abs(scaled[("ctg1_e", "ctg3_b")] - 0.3) < 0.01


def test_max_weight_matching():
    sg = FusedScaffoldGraph()
    for ctg in ["a", "b", "c", "d"]:
        sg.add_node(ctg + "_b"); sg.add_node(ctg + "_e")

    weights = {
        ("a_e", "b_b"): 10.0,
        ("b_e", "c_b"): 8.0,
        ("c_e", "d_b"): 9.0,
        ("a_e", "c_b"): 3.0,  # weak alternative
    }
    matching = sg.get_max_weight_matching(weights)
    assert len(matching) >= 2  # at least 2 edges selected


def test_cover_graph_no_cycle():
    sg = FusedScaffoldGraph()
    for ctg in ["a", "b"]:
        sg.add_node(ctg + "_b"); sg.add_node(ctg + "_e")

    weights = {("a_e", "b_b"): 10.0}
    matching = {("a_e", "b_b")}
    cover = sg.cover_graph(matching, weights)
    # Should have 3 edges: a_b-a_e (inf), b_b-b_e (inf), a_e-b_b (10.0)
    assert cover.number_of_edges() == 3
