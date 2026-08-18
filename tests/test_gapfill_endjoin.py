"""Unit tests for end-join site detection (overlapping contig ends)."""

from nearscaff.gapfill import parse_endjoins, _end_identity


def _paf(tmp_path, rows):
    """rows: [(read, gid|side, strand, qs, qe, ts, te, tlen)]"""
    p = tmp_path / "e.paf"
    p.write_text("\n".join(
        f"{r}\t1000\t{qs}\t{qe}\t{s}\t{t}\t{tl}\t{ts}\t{te}\t100\t100\t60"
        for r, t, s, qs, qe, ts, te, tl in rows) + "\n")
    return str(p)


def test_parse_endjoins_abut(tmp_path):
    """Reads hitting both flanks with extrapolated gap ~0 -> abut."""
    # flank L len 1500: aln ends at 1500 (edge). R: starts at 0.
    # est = (qs_r - 0) - (qe_l + 0) = 500 - 500 = 0
    paf = _paf(tmp_path, [
        ("r1", "g1|L", "+", 0, 500, 1000, 1500, 1500),
        ("r1", "g1|R", "+", 500, 1000, 0, 500, 1500),
        ("r2", "g1|L", "+", 0, 500, 1000, 1500, 1500),
        ("r2", "g1|R", "+", 505, 1000, 0, 495, 1500),
    ])
    assert parse_endjoins(paf) == {"g1": 0}


def test_parse_endjoins_overlap_trim(tmp_path):
    """Negative extrapolated gap -> overlap, trim from left contig end."""
    # est = (400 - 0) - (600 + 0) = -200 -> overlap 200 bp
    paf = _paf(tmp_path, [
        ("r1", "g1|L", "+", 0, 600, 900, 1500, 1500),
        ("r1", "g1|R", "+", 400, 1000, 0, 600, 1500),
        ("r2", "g1|L", "+", 0, 600, 900, 1500, 1500),
        ("r2", "g1|R", "+", 400, 1000, 0, 600, 1500),
    ])
    assert parse_endjoins(paf) == {"g1": 200}


def test_parse_endjoins_min_reads(tmp_path):
    paf = _paf(tmp_path, [
        ("r1", "g1|L", "+", 0, 500, 1000, 1500, 1500),
        ("r1", "g1|R", "+", 500, 1000, 0, 500, 1500),
    ])
    assert parse_endjoins(paf, min_reads=2) == {}


def test_end_identity():
    scaf = "AAAACCCCGGGG" + "TTTT" + "AAAACCCC"
    gb, ge, ov = 12, 12, 8
    assert _end_identity(scaf, gb, ge, ov) < 0.5
    # M + M: right 20-mer fully contained in left 20-mer -> 1.0
    m = "AAAACCCCTTTTGGGGACCC"
    scaf2 = "GGGG" + m + m
    assert _end_identity(scaf2, 24, 24, 20) >= 0.8
