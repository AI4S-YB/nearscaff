"""Tests for the Block Tree builder (nearscaff.blocktree)."""

import pytest

from nearscaff.blocktree import (
    detect_interleaved,
    build_conflict_graph,
    assign_subgenomes,
    build_block_tree,
    build_chromosome_layer,
    build_local_clusters,
    block_tree_to_json,
    block_tree_from_json,
)
from nearscaff.types import SyntenyBlock, BlockTreeNode


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def make_block(
    bid: str,
    ref_chr: str,
    ref_start: int,
    ref_end: int,
    query_contigs: set[str],
    orientation: str = "+",
    score: float = 100,
    count: int = 10,
) -> SyntenyBlock:
    """Create a minimal :class:`SyntenyBlock` for testing."""
    return SyntenyBlock(
        block_id=bid,
        ref_chr=ref_chr,
        ref_start=ref_start,
        ref_end=ref_end,
        query_contigs=set(query_contigs),
        gene_pairs=[(f"q{i}", f"r{i}", 10) for i in range(count)],
        orientation=orientation,
        score=score,
        anchor_count=count,
    )


# ---------------------------------------------------------------------------
# Tests — detect_interleaved
# ---------------------------------------------------------------------------

class TestDetectInterleaved:
    def test_basic_interleaved(self):
        b1 = make_block("b1", "Chr1", 1000, 5000, {"ctgA", "ctgB"}, "+")
        b2 = make_block("b2", "Chr1", 2000, 6000, {"ctgC", "ctgD"}, "+")
        assert detect_interleaved(b1, b2, 0.5) == "INTERLEAVED"

    def test_no_overlap_large_distance(self):
        b1 = make_block("b1", "Chr1", 1000, 5000, {"ctgA", "ctgB"}, "+")
        b3 = make_block("b3", "Chr1", 10000, 15000, {"ctgE"})
        assert detect_interleaved(b1, b3, 0.5) is None

    def test_shared_contig_not_interleaved(self):
        b1 = make_block("b1", "Chr1", 1000, 5000, {"ctgA", "ctgB"}, "+")
        b4 = make_block("b4", "Chr1", 1000, 5000, {"ctgA", "ctgX"})
        assert detect_interleaved(b1, b4, 0.5) != "INTERLEAVED"

    def test_different_chromosomes(self):
        b1 = make_block("b1", "Chr1", 1000, 5000, {"ctgA"}, "+")
        b2 = make_block("b2", "Chr2", 2000, 6000, {"ctgB"}, "+")
        assert detect_interleaved(b1, b2) is None

    def test_different_orientation(self):
        b1 = make_block("b1", "Chr1", 1000, 5000, {"ctgA"}, "+")
        b2 = make_block("b2", "Chr1", 2000, 6000, {"ctgB"}, "-")
        assert detect_interleaved(b1, b2) is None

    def test_insufficient_overlap(self):
        b1 = make_block("b1", "Chr1", 1000, 5000, {"ctgA"}, "+")
        b2 = make_block("b2", "Chr1", 4500, 9500, {"ctgB"}, "+")
        # overlap = 5000-4500 = 500, shorter = 4000 (b1), ratio = 0.125 < 0.5
        assert detect_interleaved(b1, b2, 0.5) is None

    def test_zero_length_block(self):
        b1 = make_block("b1", "Chr1", 1000, 1000, {"ctgA"}, "+")  # zero length
        b2 = make_block("b2", "Chr1", 1000, 5000, {"ctgB"}, "+")
        assert detect_interleaved(b1, b2) is None

    def test_custom_threshold(self):
        b1 = make_block("b1", "Chr1", 1000, 5000, {"ctgA"}, "+")
        b2 = make_block("b2", "Chr1", 4500, 9500, {"ctgB"}, "+")
        # overlap = 500, ratio = 0.125 — should pass at threshold=0.1
        assert detect_interleaved(b1, b2, 0.1) == "INTERLEAVED"


# ---------------------------------------------------------------------------
# Tests — build_conflict_graph
# ---------------------------------------------------------------------------

class TestBuildConflictGraph:
    def test_two_interleaved_pairs(self):
        blocks = [
            make_block("bA1", "Chr1", 1000, 5000, {"ctgA1"}, "+"),
            make_block("bA2", "Chr1", 5000, 10000, {"ctgA2"}, "+"),
            make_block("bB1", "Chr1", 1500, 4500, {"ctgB1"}, "+"),
            make_block("bB2", "Chr1", 4500, 9500, {"ctgB2"}, "+"),
        ]
        graph = build_conflict_graph(blocks, 0.5)
        assert "Chr1" in graph
        chr1 = graph["Chr1"]
        assert "bA1" in chr1
        assert "bB1" in chr1["bA1"]
        assert "bA2" in chr1
        assert "bB2" in chr1["bA2"]

    def test_no_conflicts(self):
        blocks = [
            make_block("b1", "Chr1", 1000, 5000, {"ctgA"}, "+"),
            make_block("b2", "Chr1", 10000, 15000, {"ctgB"}, "+"),
        ]
        graph = build_conflict_graph(blocks)
        assert graph == {} or all(len(adj) == 0 for adj in graph.values())

    def test_empty_blocks_list(self):
        graph = build_conflict_graph([], 0.5)
        assert graph == {}


# ---------------------------------------------------------------------------
# Tests — assign_subgenomes
# ---------------------------------------------------------------------------

class TestAssignSubgenomes:
    def test_two_subgenomes_from_conflicts(self):
        blocks = [
            make_block("bA1", "Chr1", 1000, 5000, {"ctgA1"}, "+"),
            make_block("bA2", "Chr1", 5000, 10000, {"ctgA2"}, "+"),
            make_block("bB1", "Chr1", 1500, 4500, {"ctgB1"}, "+"),
            make_block("bB2", "Chr1", 4500, 9500, {"ctgB2"}, "+"),
        ]
        conflict_graph = build_conflict_graph(blocks, 0.5)
        sg_nodes = assign_subgenomes(blocks, conflict_graph, ploidy_hint=2)
        assert len(sg_nodes) == 2
        for sg in sg_nodes:
            assert sg.level == "subgenome"

    def test_empty_blocks(self):
        sg_nodes = assign_subgenomes([], {})
        assert sg_nodes == []

    def test_no_conflicts_single_subgenome(self):
        blocks = [
            make_block("b1", "Chr1", 1000, 5000, {"ctgA"}, "+"),
            make_block("b2", "Chr1", 10000, 15000, {"ctgB"}, "+"),
        ]
        conflict_graph = build_conflict_graph(blocks)
        sg_nodes = assign_subgenomes(blocks, conflict_graph)
        assert len(sg_nodes) == 1


# ---------------------------------------------------------------------------
# Tests — build_block_tree (integration)
# ---------------------------------------------------------------------------

class TestBuildBlockTree:
    def test_two_subgenomes(self):
        blocks = [
            make_block("bA1", "Chr1", 1000, 5000, {"ctgA1"}, "+", 100),
            make_block("bA2", "Chr1", 5000, 10000, {"ctgA2"}, "+", 80),
            make_block("bB1", "Chr1", 1500, 4500, {"ctgB1"}, "+", 90),
            make_block("bB2", "Chr1", 4500, 9500, {"ctgB2"}, "+", 85),
        ]
        root = build_block_tree(blocks, ["Chr1"], ploidy_hint=2)
        subgenomes = list(root.iter_level("subgenome"))
        assert len(subgenomes) == 2

        # Verify tree structure
        for sg in subgenomes:
            assert len(sg.children) >= 1  # at least one chromosome
            for chr_node in sg.children:
                assert chr_node.level == "chromosome"
                assert chr_node.ref_chr == "Chr1"
                assert len(chr_node.children) >= 1  # at least one local cluster

    def test_empty_blocks_graceful(self):
        root = build_block_tree([], ["Chr1"])
        assert root.level == "root"
        assert root.children == []

    def test_single_block(self):
        blocks = [make_block("b1", "Chr1", 1000, 5000, {"ctg1"}, "+", 100)]
        root = build_block_tree(blocks, ["Chr1"])
        assert len(list(root.iter_level("subgenome"))) >= 1

    def test_multiple_chromosomes(self):
        blocks = [
            make_block("b1", "Chr1", 1000, 5000, {"ctgA"}, "+", 100),
            make_block("b2", "Chr2", 1000, 5000, {"ctgB"}, "+", 90),
        ]
        root = build_block_tree(blocks, ["Chr1", "Chr2"])
        subgenomes = list(root.iter_level("subgenome"))
        assert len(subgenomes) == 1
        sg = subgenomes[0]
        chr_names = {c.ref_chr for c in sg.children}
        assert chr_names == {"Chr1", "Chr2"}

    def test_local_clusters_merge(self):
        """Adjacent same-orientation blocks close together should merge."""
        blocks = [
            make_block("b1", "Chr1", 1000, 2000, {"ctgA"}, "+", 100),
            make_block("b2", "Chr1", 2100, 3000, {"ctgA"}, "+", 90),
        ]
        root = build_block_tree(blocks, ["Chr1"])
        locals = list(root.iter_level("local"))
        # gap = 2100 - 2000 = 100 < 500000, same orientation → merged
        assert len(locals) == 1
        assert locals[0].ref_start == 1000
        assert locals[0].ref_end == 3000

    def test_local_clusters_split_orientation(self):
        """Adjacent blocks with different orientations should NOT merge."""
        blocks = [
            make_block("b1", "Chr1", 1000, 2000, {"ctgA"}, "+", 100),
            make_block("b2", "Chr1", 2100, 3000, {"ctgA"}, "-", 90),
        ]
        root = build_block_tree(blocks, ["Chr1"])
        locals = list(root.iter_level("local"))
        assert len(locals) == 2


# ---------------------------------------------------------------------------
# Tests — JSON roundtrip
# ---------------------------------------------------------------------------

class TestBlockTreeJson:
    def test_roundtrip(self):
        root = BlockTreeNode(
            node_id="root",
            level="subgenome",
        )
        child = BlockTreeNode(
            node_id="chr1",
            level="chromosome",
            ref_chr="Chr1",
            ref_start=0,
            ref_end=1000000,
            confidence=0.9,
        )
        root.add_child(child)

        json_str = block_tree_to_json(root)
        restored = block_tree_from_json(json_str)

        assert restored.node_id == "root"
        assert restored.level == "subgenome"
        assert restored.children[0].node_id == "chr1"
        assert restored.children[0].confidence == 0.9

    def test_roundtrip_preserves_flags(self):
        root = BlockTreeNode(
            node_id="root",
            level="root",
            flags={"preliminary", "low_confidence"},
        )
        json_str = block_tree_to_json(root)
        restored = block_tree_from_json(json_str)
        assert restored.flags == {"preliminary", "low_confidence"}

    def test_roundtrip_nested_tree(self):
        """Full 3-level tree roundtrip."""
        blocks = [
            make_block("bA1", "Chr1", 1000, 5000, {"ctgA1"}, "+", 100, count=3),
            make_block("bB1", "Chr1", 1500, 4500, {"ctgB1"}, "+", 90, count=4),
        ]
        root = build_block_tree(blocks, ["Chr1"], ploidy_hint=2)
        json_str = block_tree_to_json(root)
        restored = block_tree_from_json(json_str)

        assert restored.node_id == "root"
        assert restored.level == "root"
        # Should have subgenomes
        restored_sgs = list(restored.iter_level("subgenome"))
        assert len(restored_sgs) == len(list(root.iter_level("subgenome")))


# ---------------------------------------------------------------------------
# Tests — layer construction functions
# ---------------------------------------------------------------------------

class TestLayerFunctions:
    def test_build_chromosome_layer(self):
        blocks = [
            make_block("b1", "Chr1", 1000, 5000, {"ctgA"}, "+", 100),
            make_block("b2", "Chr1", 6000, 10000, {"ctgA"}, "+", 80),
            make_block("b3", "Chr2", 1000, 4000, {"ctgB"}, "+", 90),
        ]
        sg_node = BlockTreeNode(
            node_id="sg_0",
            level="subgenome",
            query_contigs=["ctgA", "ctgB"],
            gene_pairs=[gp for b in blocks for gp in b.gene_pairs],
            synteny_score=90.0,
            orientation="+",
        )
        sg_node._block_ids = ["b1", "b2", "b3"]

        chr_nodes = build_chromosome_layer(sg_node, blocks)

        assert len(chr_nodes) == 2
        chr_ids = {n.node_id for n in chr_nodes}
        assert "sg_0_Chr1" in chr_ids
        assert "sg_0_Chr2" in chr_ids
        assert len(sg_node.children) == 2

    def test_build_local_clusters(self):
        blocks = [
            make_block("b1", "Chr1", 1000, 2000, {"ctgA"}, "+", 100),
            make_block("b2", "Chr1", 2100, 3000, {"ctgA"}, "+", 90),
            make_block("b3", "Chr1", 100000, 101000, {"ctgA"}, "+", 80),
        ]
        chr_node = BlockTreeNode(
            node_id="sg_0_Chr1",
            level="chromosome",
            ref_chr="Chr1",
            ref_start=1000,
            ref_end=101000,
            query_contigs=["ctgA"],
            gene_pairs=[gp for b in blocks for gp in b.gene_pairs],
            synteny_score=90.0,
            orientation="+",
        )
        chr_node._block_ids = ["b1", "b2", "b3"]

        local_nodes = build_local_clusters(chr_node, blocks, max_inter_block_gap=50_000)

        # b1 and b2 should merge (gap=100), b3 is far away
        assert len(local_nodes) == 2
        assert len(chr_node.children) == 2
        assert local_nodes[0].ref_start == 1000
        assert local_nodes[0].ref_end == 3000
