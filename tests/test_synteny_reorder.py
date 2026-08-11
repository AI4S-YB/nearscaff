"""Unit tests for nearscaff synteny-consistent contig reordering."""
import os
import tempfile

from nearscaff.config import NearscaffConfig
from nearscaff.pipeline import _write_anchors_tsv, _load_anchor_coords


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
