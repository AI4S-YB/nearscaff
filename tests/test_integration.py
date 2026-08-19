"""End-to-end integration tests for the nearscaff pipeline.

These tests require miniprot and minimap2 to be installed.
"""
import os
import sys
import json
import subprocess

import pytest
import networkx as nx

from nearscaff.config import NearscaffConfig
from nearscaff.pipeline import run_stage0, run_stage1, run_full
from nearscaff.blocktree import block_tree_from_json
from nearscaff.scaffold_graph import FusedScaffoldGraph
from nearscaff.types import EdgeType, SyntenyBlock
from nearscaff.blocktree import build_block_tree, block_tree_to_json
from nearscaff.agp import AGPReader

# Package root (parent of tests/) so CLI subprocesses can import nearscaff
PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run_cli(*cli_args, timeout=120):
    env = dict(os.environ)
    env["PYTHONPATH"] = PKG_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", "nearscaff", *cli_args],
        capture_output=True, text=True, timeout=timeout, env=env,
    )


# ---------------------------------------------------------------------------
# External tool checks
# ---------------------------------------------------------------------------

def _has_tool(name):
    return subprocess.run(["which", name], capture_output=True).returncode == 0


HAS_MINIPROT = _has_tool("miniprot")
HAS_MINIMAP2 = _has_tool("minimap2")
HAS_SAMTOOLS = _has_tool("samtools")


# ---------------------------------------------------------------------------
# Unit-level pipeline tests (no external tools needed)
# ---------------------------------------------------------------------------

class TestPipelineUnit:
    """Pipeline logic tests that don't need external tools."""

    def test_build_block_tree_from_synthetic_blocks(self):
        """Block Tree builds successfully from synthetic synteny blocks."""
        blocks = [
            SyntenyBlock(
                block_id="b1", ref_chr="Chr1",
                ref_start=1000, ref_end=5000,
                query_contigs={"ctgA"}, gene_pairs=[("q1", "r1", 10)],
                orientation="+", score=100, anchor_count=10,
            ),
            SyntenyBlock(
                block_id="b2", ref_chr="Chr1",
                ref_start=5100, ref_end=10000,
                query_contigs={"ctgA"}, gene_pairs=[("q2", "r2", 10)],
                orientation="+", score=90, anchor_count=10,
            ),
            SyntenyBlock(
                block_id="b3", ref_chr="Chr1",
                ref_start=15000, ref_end=20000,
                query_contigs={"ctgB"}, gene_pairs=[("q3", "r3", 10)],
                orientation="+", score=80, anchor_count=10,
            ),
        ]
        root = build_block_tree(blocks, ["Chr1"], ploidy_hint=1)
        assert root.level == "root"

        subgenomes = list(root.iter_level("subgenome"))
        assert len(subgenomes) == 1

        local_clusters = list(root.iter_level("local"))
        assert len(local_clusters) >= 1

        # JSON roundtrip
        json_str = block_tree_to_json(root)
        restored = block_tree_from_json(json_str)
        assert restored.level == root.level
        assert len(restored.children) == len(root.children)

    def test_scaffold_graph_full_cycle(self):
        """Scaffold Graph: add nodes -> add edges -> fuse -> match -> cover."""
        sg = FusedScaffoldGraph()
        for ctg in ["a", "b", "c"]:
            sg.add_node(ctg + "_b")
            sg.add_node(ctg + "_e")

        sg.add_fused_edge("a_e", "b_b", weight=0.9, edge_type=EdgeType.PROTEIN_SYNTENY,
                          gap_size=1000, source="synteny")
        sg.add_fused_edge("b_e", "c_b", weight=0.7, edge_type=EdgeType.NUCLEOTIDE_CHAIN,
                          gap_size=800, source="nucleotide")

        fused = sg.fuse_weights()
        assert len(fused) == 2

        scaled = sg.best_buddy_scale(fused)
        matching = sg.get_max_weight_matching(scaled)
        cover = sg.cover_graph(matching, scaled)

        # Should have a-e_b, b_b-b_e (intra), b_e-c_b
        assert isinstance(cover, nx.Graph)

    def test_agp_roundtrip(self):
        """AGP parsing and formatting roundtrip."""
        from nearscaff.agp import AGPSeqLine, AGPGapLine, AGPReader, AGPWriter

        lines = [
            AGPSeqLine("scaf1", 1, 1000, 1, "W", "ctg1", 1, 1000, "+"),
            AGPGapLine("scaf1", 1001, 1100, 2, "N", 100, "scaffold", "yes", "align_genus"),
            AGPSeqLine("scaf1", 1101, 2000, 3, "W", "ctg2", 1, 900, "+"),
        ]
        writer = AGPWriter()
        text = writer.format(lines)

        reader = AGPReader()
        parsed = reader.parse(text)
        assert len(parsed) == 3
        assert parsed[0].component_id == "ctg1"
        assert parsed[1].gap_length == 100

    def test_append_unplaced_singletons(self):
        """--keep-unplaced: contigs absent from the AGP become singleton
        scaffolds named after the contig; placed contigs untouched."""
        from nearscaff.agp import AGPSeqLine, AGPGapLine
        from nearscaff.pipeline import _append_unplaced_singletons

        lines = [
            AGPSeqLine("nearscaff_0001", 1, 1000, 1, "W", "ctg1", 1, 1000, "+"),
            AGPGapLine("nearscaff_0001", 1001, 1100, 2, "N", 100,
                       "scaffold", "yes", "align_genus"),
            AGPSeqLine("nearscaff_0001", 1101, 2000, 3, "W", "ctg2", 1, 900, "+"),
        ]
        lengths = {"ctg1": 1000, "ctg2": 900, "ctg3": 500, "ctg4": 50}
        out = _append_unplaced_singletons(lines, lengths, 0)
        assert len(out) == 5          # ctg3 + ctg4 appended
        singletons = out[3:]
        assert [l.object_name for l in singletons] == ["ctg3", "ctg4"]
        assert all(l.component_id == l.object_name and l.orientation == "+"
                   for l in singletons)
        # min_len filters small contigs
        out = _append_unplaced_singletons(lines, lengths, 100)
        assert [l.object_name for l in out[3:]] == ["ctg3"]


# ---------------------------------------------------------------------------
# Integration tests (require external tools)
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestPipelineIntegration:
    """Full pipeline tests requiring miniprot + minimap2."""

    @pytest.mark.skipif(not HAS_MINIPROT, reason="miniprot not installed")
    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    def test_stage0_anchors_and_blocktree(self, tiny_ref_fa, tiny_gff3, tiny_query_fa, tmp_path):
        """Stage 0: protein anchoring + synteny -> Block Tree."""
        config = NearscaffConfig()
        config.threads = 2
        config.keep_intermediate = True
        config.synteny.min_cluster_size = 2
        config.output_dir = str(tmp_path)

        root = run_stage0(config, tiny_ref_fa, tiny_gff3, tiny_query_fa,
                          str(tmp_path))

        assert root is not None, "Stage 0 should produce a Block Tree"
        tree_path = os.path.join(str(tmp_path), "block_tree.json")
        assert os.path.exists(tree_path), f"Block Tree JSON should exist at {tree_path}"

        # Verify JSON is valid
        with open(tree_path) as f:
            data = json.load(f)
        assert data["level"] == "root"

        # Verify tree structure
        local_clusters = list(root.iter_level("local"))
        assert len(local_clusters) >= 1, "Should have at least 1 local cluster"

    @pytest.mark.skipif(not HAS_MINIPROT, reason="miniprot not installed")
    @pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    def test_stage0_writes_alignment_cache(self, tiny_ref_fa, tiny_gff3, tiny_query_fa, tmp_path):
        """Stage0 should persist per-contig alignments for stage1 refine reuse."""
        config = NearscaffConfig()
        config.threads = 2
        config.synteny.min_cluster_size = 2
        config.output_dir = str(tmp_path)
        run_stage0(config, tiny_ref_fa, tiny_gff3, tiny_query_fa, str(tmp_path))
        cache_path = os.path.join(str(tmp_path), "intermediate", "contig_alignments.tsv")
        assert os.path.exists(cache_path), "stage0 should write the alignment cache"
        from nearscaff.nucleotide import read_align_cache
        cache = read_align_cache(cache_path)
        assert len(cache) >= 1
        entry = next(iter(cache.values()))
        assert set(entry) == {"chr", "r_start", "r_end", "strand",
                              "mapq", "hitlen", "identity"}

    @pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    def test_stage1_scaffold(self, tiny_ref_fa, tiny_query_fa, tiny_block_tree, tmp_path):
        """Stage 1: nucleotide alignment + scaffold graph -> AGP + FASTA."""
        config = NearscaffConfig()
        config.threads = 2
        config.nucleotide.region_margin = 1000
        config.nucleotide.preset = "asm5"
        config.output_dir = str(tmp_path)

        agp_path = run_stage1(config, tiny_block_tree, tiny_ref_fa, tiny_query_fa,
                              str(tmp_path))

        assert agp_path is not None, "Stage 1 should produce an AGP"
        assert os.path.exists(agp_path), f"AGP should exist at {agp_path}"

        # Verify AGP is valid
        with open(agp_path) as f:
            content = f.read()
        assert "nearscaff_" in content

        reader = AGPReader()
        lines = reader.parse(content)
        assert len(lines) >= 1, "AGP should contain at least one line"

        # Verify scaffold FASTA
        fasta_path = os.path.join(str(tmp_path), "nearscaff.scaffolds.fa")
        assert os.path.exists(fasta_path), "Stage 1 should produce scaffold FASTA"
        with open(fasta_path) as f:
            fasta_text = f.read()
        assert fasta_text.startswith(">"), "scaffold FASTA should not be empty"
        # Both query contigs must appear in the output (placed or appended)
        total_len = sum(len(l.strip()) for l in fasta_text.splitlines()
                        if not l.startswith(">"))
        # 80000 bp of query sequence plus any scaffold gap Ns
        assert total_len >= 80000, \
            f"scaffold FASTA should cover all query sequence (got {total_len})"

    @pytest.mark.skipif(not HAS_MINIPROT, reason="miniprot not installed")
    @pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    @pytest.mark.slow
    def test_full_pipeline(self, tiny_ref_fa, tiny_gff3, tiny_query_fa, tmp_path):
        """End-to-end pipeline run with synthetic data."""
        config = NearscaffConfig()
        config.threads = 2
        config.keep_intermediate = True
        config.synteny.min_cluster_size = 2  # easier to form clusters
        config.output_dir = str(tmp_path)

        agp_path = run_full(config, tiny_ref_fa, tiny_gff3, tiny_query_fa,
                            str(tmp_path))

        assert agp_path is not None, "Pipeline should produce an AGP file"
        assert os.path.exists(agp_path)

        # Check output files
        assert os.path.exists(os.path.join(str(tmp_path), "block_tree.json"))
        assert os.path.exists(os.path.join(str(tmp_path), "nearscaff.log"))
        assert os.path.exists(os.path.join(str(tmp_path), "nearscaff.agp"))
        assert os.path.exists(os.path.join(str(tmp_path), "nearscaff.scaffolds.fa"))

        # Verify AGP validity
        with open(agp_path) as f:
            content = f.read()
        reader = AGPReader()
        lines = reader.parse(content)
        assert len(lines) >= 1

        # Verify scaffold FASTA covers the full query (placed + appended)
        with open(os.path.join(str(tmp_path), "nearscaff.scaffolds.fa")) as f:
            fasta_text = f.read()
        total_len = sum(len(l.strip()) for l in fasta_text.splitlines()
                        if not l.startswith(">"))
        assert total_len >= 80000

    @pytest.mark.skipif(not HAS_MINIPROT, reason="miniprot not installed")
    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    def test_cli_run(self, tiny_ref_fa, tiny_gff3, tiny_query_fa, tmp_path):
        """CLI: nearscaff run subcommand."""
        result = _run_cli(
            "run",
            "-r", tiny_ref_fa, "-g", tiny_gff3, "-q", tiny_query_fa,
            "-o", str(tmp_path), "-t", "2", "--cscore", "0.3",
            "--min-cluster-size", "2", "--nucleotide-passes", "asm5",
        )
        assert result.returncode == 0, f"CLI run failed:\n{result.stderr}\n{result.stdout}"
        assert os.path.exists(os.path.join(str(tmp_path), "block_tree.json"))
        assert os.path.exists(os.path.join(str(tmp_path), "nearscaff.agp"))
        assert os.path.exists(os.path.join(str(tmp_path), "nearscaff.scaffolds.fa"))

    @pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    def test_cli_scaffold(self, tiny_ref_fa, tiny_query_fa, tiny_block_tree, tmp_path):
        """CLI: nearscaff scaffold subcommand."""
        result = _run_cli(
            "scaffold",
            "-b", tiny_block_tree, "-r", tiny_ref_fa, "-q", tiny_query_fa,
            "-o", str(tmp_path), "-t", "2", "--margin", "1000",
        )
        assert result.returncode == 0, f"CLI scaffold failed:\n{result.stderr}\n{result.stdout}"
        assert os.path.exists(os.path.join(str(tmp_path), "nearscaff.agp"))
        assert os.path.exists(os.path.join(str(tmp_path), "nearscaff.scaffolds.fa"))

    def test_cli_help(self):
        """CLI help works without external tools."""
        result = _run_cli("--help")
        assert result.returncode == 0
        assert "run" in result.stdout
        assert "scaffold" in result.stdout

    @pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    def test_stage1_writes_align_cache(self, tiny_ref_fa, tiny_query_fa, tiny_block_tree, tmp_path):
        config = NearscaffConfig()
        config.threads = 2
        config.nucleotide.region_margin = 1000
        config.nucleotide.preset = "asm5"
        config.nucleotide.reuse_ref_index = True
        config.output_dir = str(tmp_path)
        agp_path = run_stage1(config, tiny_block_tree, tiny_ref_fa, tiny_query_fa,
                              str(tmp_path))
        assert agp_path is not None and os.path.exists(agp_path)
        cache_path = os.path.join(str(tmp_path), "intermediate", "contig_alignments.tsv")
        assert os.path.exists(cache_path), f"align cache should exist at {cache_path}"

    @pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
    def test_build_ref_index_idempotent(self, tiny_ref_fa, tmp_path):
        from nearscaff.nucleotide import build_ref_index, index_path_for
        idx = build_ref_index(tiny_ref_fa, "asm5", str(tmp_path), threads=2)
        assert os.path.exists(idx)
        assert idx == index_path_for(tiny_ref_fa, "asm5", str(tmp_path))
        # idempotent: second call reuses, still succeeds
        idx2 = build_ref_index(tiny_ref_fa, "asm5", str(tmp_path), threads=2)
        assert idx == idx2

    @pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
    def test_align_to_full_reference_returns_paf(self, tiny_ref_fa, tiny_query_fa, tmp_path):
        from nearscaff.nucleotide import build_ref_index, align_to_full_reference, parse_nucleotide_paf
        idx = build_ref_index(tiny_ref_fa, "asm5", str(tmp_path), threads=2)
        paf = align_to_full_reference(idx, tiny_query_fa, preset="asm5",
                                      secondary=5, with_cigar=False, threads=2)
        entries = parse_nucleotide_paf(paf)
        assert len(entries) >= 1
        assert all(e["query"] in ("ctg1", "ctg2") for e in entries)

    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    def test_extract_contigs_faidx_path(self, tmp_path):
        from nearscaff.nucleotide import ensure_query_faid, _extract_contigs
        src = tmp_path / "q.fa"
        src.write_text(">ctg1\nACGTACGT\n>ctg2\nTTTTGGGG\n")
        qpath = str(src)
        assert ensure_query_faid(qpath) is True
        out = str(tmp_path / "out.fa")
        _extract_contigs(qpath, ["ctg2"], out)
        txt = open(out).read()
        assert ">ctg2" in txt and ">ctg1" not in txt

    @pytest.mark.skipif(not HAS_MINIMAP2, reason="minimap2 not installed")
    @pytest.mark.skipif(not HAS_SAMTOOLS, reason="samtools not installed")
    def test_align_unplaced_whole_ref_vs_per_region(self, tiny_ref_fa, tiny_query_fa, tmp_path):
        """Whole-reference path and per-region fallback must find the same queries."""
        from nearscaff.nucleotide import build_ref_index, ensure_query_faid
        from nearscaff.pipeline import _align_unplaced_to_scaffolds
        from types import SimpleNamespace
        ensure_query_faid(tiny_query_fa)
        idx = build_ref_index(tiny_ref_fa, "asm5", str(tmp_path), threads=2)
        region = SimpleNamespace(scaffold_idx=0, ref_chr="Chr1",
                                 ref_start=5000, ref_end=95000,
                                 contigs=["ctg1"])
        unplaced = {"ctg2"}
        whole = _align_unplaced_to_scaffolds(
            unplaced, [region], tiny_ref_fa, tiny_query_fa,
            preset="asm5", margin=1000, threads=2,
            ref_index=idx, secondary=5)
        per_region = _align_unplaced_to_scaffolds(
            unplaced, [region], tiny_ref_fa, tiny_query_fa,
            preset="asm5", margin=1000, threads=2,
            ref_index=None, secondary=5)
        assert {e["query"] for e in whole} == {e["query"] for e in per_region}
        assert "ctg2" in {e["query"] for e in whole}
