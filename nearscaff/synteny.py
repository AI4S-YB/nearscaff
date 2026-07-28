"""Synteny engine — C-score filtering of gene anchors.

Adapted from JCVI jcvi/compara/blastfilter.py.
"""
from collections import defaultdict
from nearscaff.types import GeneAnchor


def filter_cscore(anchors: list[GeneAnchor], cscore: float = 0.7):
    """Filter anchors by C-score.

    C-score(a,b) = score(a,b) / max(best_score(a), best_score(b))
    Reference: Putnam et al. (2008) Science
    """
    best_query = defaultdict(float)
    best_subject = defaultdict(float)
    for a in anchors:
        key_q = (a.query_contig, a.query_gene)
        key_s = (a.ref_chr, a.ref_gene)
        if a.score > best_query[key_q]:
            best_query[key_q] = a.score
        if a.score > best_subject[key_s]:
            best_subject[key_s] = a.score

    filtered = []
    for a in anchors:
        key_q = (a.query_contig, a.query_gene)
        key_s = (a.ref_chr, a.ref_gene)
        best = max(best_query[key_q], best_subject[key_s])
        cs = a.score / best if best > 0 else 0
        a._cscore = cs
        if cs >= cscore:
            filtered.append(a)
    return filtered
