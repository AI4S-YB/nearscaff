"""Tests for the synteny engine (C-score filtering)."""
from nearscaff.synteny import filter_cscore
from nearscaff.types import GeneAnchor


def make_anchor(qc, qg, rc, rg, score, q_start=0, r_start=0):
    return GeneAnchor(
        query_contig=qc, query_gene=qg, ref_chr=rc, ref_gene=rg,
        q_start=q_start, q_end=q_start+100, r_start=r_start, r_end=r_start+100,
        strand='+', score=float(score), identity=1.0,
        n_exons=1, n_frameshifts=0, n_stop_codons=0,
    )


def test_filter_cscore():
    a1 = make_anchor("ctg1", "gA", "Chr1", "gX", 800)
    a2 = make_anchor("ctg1", "gA", "Chr2", "gY", 600)
    a3 = make_anchor("ctg1", "gB", "Chr1", "gZ", 400)
    # gA best=800, gB best=400
    # a1 cscore=800/800=1.0, a2 cscore=600/800=0.75, a3 cscore=400/400=1.0
    filtered = filter_cscore([a1, a2, a3], cscore=0.7)
    assert len(filtered) == 3
    filtered_strict = filter_cscore([a1, a2, a3], cscore=0.8)
    assert len(filtered_strict) == 2  # a2 excluded (0.75 < 0.8)
