"""Unit tests for the Stage 1 nucleotide-extension helpers (no external tools)."""
from nearscaff.config import NucleotideConfig


def test_nucleotide_config_defaults():
    cfg = NucleotideConfig()
    assert cfg.secondary_alignments == 5
    assert cfg.reuse_ref_index is True
    assert cfg.index_dir is None
    # asm10 is redundant (looser asm20 subsumes it); default drops it for speed.
    assert cfg.nucleotide_passes == ["asm5", "asm20"]


from nearscaff.nucleotide import write_align_cache, read_align_cache


def test_align_cache_roundtrip(tmp_path):
    path = str(tmp_path / "cache.tsv")
    cache = {
        "ctg1": {"chr": "Chr1", "r_start": 100, "r_end": 5000, "strand": "+",
                 "mapq": 60, "hitlen": 4900, "identity": 0.987},
        "ctg2": {"chr": "Chr2", "r_start": 200, "r_end": 3000, "strand": "-",
                 "mapq": 40, "hitlen": 2800, "identity": 0.91},
    }
    write_align_cache(path, cache)
    loaded = read_align_cache(path)
    assert set(loaded) == {"ctg1", "ctg2"}
    assert loaded["ctg1"]["chr"] == "Chr1"
    assert loaded["ctg1"]["r_end"] == 5000
    assert loaded["ctg2"]["strand"] == "-"
    assert abs(loaded["ctg1"]["identity"] - 0.987) < 1e-6


def test_read_align_cache_missing_returns_empty(tmp_path):
    assert read_align_cache(str(tmp_path / "nope.tsv")) == {}


import os
from nearscaff.nucleotide import index_path_for


def test_index_path_for(tmp_path):
    p = index_path_for("/data/genomes/ref.fa", "asm5", str(tmp_path))
    assert os.path.basename(p) == "ref.fa.asm5.mmi"
    assert p.startswith(str(tmp_path))


from nearscaff.nucleotide import _extract_contigs


def test_extract_contigs_linear_fallback(tmp_path):
    src = tmp_path / "q.fa"
    src.write_text(">ctg1\nACGTACGT\n>ctg2\nTTTTGGGG\n")
    out = str(tmp_path / "out.fa")
    # base 环境无 samtools → _faidx_available 返回 False → 走线性扫描
    _extract_contigs(str(src), ["ctg2"], out)
    txt = open(out).read()
    assert ">ctg2" in txt and ">ctg1" not in txt


from types import SimpleNamespace
from nearscaff.scaffold_graph import FusedScaffoldGraph
from nearscaff.pipeline import _add_extension_edges


def _entry(query, rchr, rs, re_, nmatch, hitlen, mapq):
    return {"query": query, "qlen": 3000, "qstart": 0, "qend": 3000,
            "strand": "+", "ref_chr": rchr, "rlen": 100000,
            "rstart": rs, "rend": re_, "nmatch": nmatch,
            "hitlen": hitlen, "mapq": mapq}


def test_add_extension_edges_uses_secondary_near_scaffold():
    """主比对在 Chr2、次级在 Chr1 的 scaffold 端点附近 → 仍应建立延伸边。"""
    sg = FusedScaffoldGraph()
    for c in ("scaf_a", "ctg_new"):
        sg.add_node(c + "_b")
        sg.add_node(c + "_e")
    region = SimpleNamespace(scaffold_idx=0, ref_chr="Chr1",
                             ref_start=1000, ref_end=5000, contigs=["scaf_a"])
    entries = [
        # 强主比对，但在 Chr2（远离 scaffold）
        _entry("ctg_new", "Chr2", 50000, 53000, 2900, 3000, 60),
        # 弱次级比对，紧邻 Chr1 scaffold 右端（gap=200）
        _entry("ctg_new", "Chr1", 5200, 8000, 800, 2800, 20),
    ]
    n = _add_extension_edges(sg, entries, [region],
                             {"scaf_a": 10000, "ctg_new": 3000},
                             gap_min=0, gap_max=500000)
    assert n == 1  # 旧行为(best_per_query)会因主比对在 Chr2 而漏掉 → 返回 0


def test_add_extension_edges_rejects_beyond_gap_max():
    """距 scaffold 端点超过 gap_max 的比对不应建边（包容循环的过滤闸门）。"""
    sg = FusedScaffoldGraph()
    for c in ("scaf_a", "ctg_far"):
        sg.add_node(c + "_b")
        sg.add_node(c + "_e")
    region = SimpleNamespace(scaffold_idx=0, ref_chr="Chr1",
                             ref_start=1000, ref_end=5000, contigs=["scaf_a"])
    # ctg_far 在 Chr1 但距 scaffold 右端 ~2Mb，超过 gap_max=500000
    entries = [_entry("ctg_far", "Chr1", 2010000, 2013000, 2900, 3000, 60)]
    n = _add_extension_edges(sg, entries, [region],
                             {"scaf_a": 10000, "ctg_far": 3000},
                             gap_min=0, gap_max=500000)
    assert n == 0


from nearscaff.pipeline import _refine_contig_coordinates
from nearscaff.config import NearscaffConfig


def test_refine_uses_cache_without_realign(tmp_path, monkeypatch):
    """有缓存时，_refine_contig_coordinates 不应调用 _extract_contigs/minimap2。"""
    config = NearscaffConfig()
    contig_ref = {"ctg1": ("Chr1", 0, 1000), "ctg2": ("Chr1", 0, 1000)}
    cache = {
        "ctg1": {"chr": "Chr1", "r_start": 100, "r_end": 4900, "strand": "+",
                 "mapq": 60, "hitlen": 4800, "identity": 0.99},
        "ctg2": {"chr": "Chr1", "r_start": 5100, "r_end": 9000, "strand": "-",
                 "mapq": 55, "hitlen": 3900, "identity": 0.97},
    }
    called = {"n": 0}

    def _boom(*a, **k):
        called["n"] += 1
        raise AssertionError("should not align when cache covers all contigs")

    # _extract_contigs is imported locally from nearscaff.nucleotide inside the
    # function, so patch the source module (picked up by the local import).
    monkeypatch.setattr("nearscaff.nucleotide._extract_contigs", _boom)
    strand, mapq = _refine_contig_coordinates(
        contig_ref, {"ctg1": 5000, "ctg2": 4000},
        "ref.fa", "query.fa", config, threads=1,
        ref_index=None, align_cache=cache)
    assert called["n"] == 0
    assert strand == {"ctg1": "+", "ctg2": "-"}
    assert contig_ref["ctg1"] == ("Chr1", 100, 4900)


import subprocess


def test_cli_scaffold_help_lists_new_options():
    result = subprocess.run(["nearscaff", "scaffold", "--help"],
                            capture_output=True, text=True)
    assert result.returncode == 0
    assert "--secondary-alignments" in result.stdout
    assert "--no-reuse-index" in result.stdout


from nearscaff.nucleotide import _align_full_ref_cmd
from nearscaff.pipeline import _best_contig_alignments


def test_best_contig_alignments_picks_best_by_nmatch_with_full_fields():
    """每 contig 取 nmatch 最佳，且保留 strand/mapq/hitlen 供 stage1 复用。"""
    entries = [
        {"query": "ctg1", "ref_chr": "Chr1", "rstart": 100, "rend": 5000, "strand": "+",
         "nmatch": 4000, "hitlen": 4500, "mapq": 60},
        {"query": "ctg1", "ref_chr": "Chr1", "rstart": 60000, "rend": 65000, "strand": "-",
         "nmatch": 2000, "hitlen": 2500, "mapq": 30},  # weaker -> ignored
        {"query": "ctg2", "ref_chr": "Chr2", "rstart": 200, "rend": 3000, "strand": "-",
         "nmatch": 2500, "hitlen": 2800, "mapq": 55},
    ]
    best = _best_contig_alignments(entries)
    assert set(best) == {"ctg1", "ctg2"}
    c1 = best["ctg1"]  # (chr, r_start, r_end, identity, strand, mapq, hitlen)
    assert (c1[0], c1[1], c1[2]) == ("Chr1", 100, 5000)   # best by nmatch (4000)
    assert (c1[4], c1[5], c1[6]) == ("+", 60, 4500)        # strand, mapq, hitlen
    assert abs(c1[3] - 4000 / 4500) < 1e-6                  # identity = nmatch/hitlen
    assert best["ctg2"][4] == "-"                           # strand preserved for ctg2


def test_align_full_ref_cmd_secondary_zero_disables_secondaries():
    cmd = _align_full_ref_cmd("ref.mmi", "q.fa", "asm5",
                              secondary=0, with_cigar=False, threads=2)
    assert "--secondary=no" in cmd
    assert "-N" not in cmd


def test_align_full_ref_cmd_secondary_positive_uses_N():
    cmd = _align_full_ref_cmd("ref.mmi", "q.fa", "asm5",
                              secondary=5, with_cigar=False, threads=2)
    assert "-N" in cmd
    assert cmd[cmd.index("-N") + 1] == "5"
    assert "--secondary=no" not in cmd


def test_align_full_ref_cmd_with_cigar_adds_c():
    cmd = _align_full_ref_cmd("ref.mmi", "q.fa", "asm5",
                              secondary=5, with_cigar=True, threads=2)
    assert "-c" in cmd


def test_refine_missing_path_disables_secondaries(tmp_path, monkeypatch):
    """无缓存 contig 走精修时，应禁用次级比对（secondary=0）——精修只取最佳比对。"""
    from nearscaff.pipeline import _refine_contig_coordinates
    from nearscaff.config import NearscaffConfig
    config = NearscaffConfig()
    contig_ref = {"ctg1": ("Chr1", 0, 1000)}  # not in cache -> missing -> realign
    captured = {}

    def _fake_align(ref_index, query_fa, preset, secondary, with_cigar, threads):
        captured["secondary"] = secondary
        captured["with_cigar"] = with_cigar
        return ""  # empty PAF -> no coordinate update

    monkeypatch.setattr("nearscaff.nucleotide.align_to_full_reference", _fake_align)
    monkeypatch.setattr("nearscaff.nucleotide._extract_contigs", lambda *a, **k: None)
    _refine_contig_coordinates(contig_ref, {"ctg1": 1000}, "ref.fa", "query.fa",
                               config, threads=1, ref_index="ref.mmi",
                               align_cache={})
    assert captured.get("secondary") == 0, "refine should disable secondaries"
    assert captured.get("with_cigar") is True


def test_extract_contigs_faidx_e2big_falls_back_to_linear(tmp_path, monkeypatch):
    """faidx 命令行过长(E2BIG, e.g. 50万 contig)时应回退线性扫描，不能崩。"""
    import nearscaff.nucleotide as nuc
    import subprocess
    src = tmp_path / "q.fa"
    src.write_text(">ctg1\nACGTACGT\n>ctg2\nTTTTGGGG\n")
    out = str(tmp_path / "out.fa")
    monkeypatch.setattr(nuc, "_faidx_available", lambda p: True)

    def _e2big(*a, **k):
        raise OSError(7, "Argument list too long")

    monkeypatch.setattr(subprocess, "run", _e2big)
    nuc._extract_contigs(str(src), ["ctg2"], out)  # must not raise
    txt = open(out).read()
    assert ">ctg2" in txt and ">ctg1" not in txt
