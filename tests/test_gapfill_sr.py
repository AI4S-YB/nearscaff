"""Unit tests for short-read PE-span gap size estimation."""

from nearscaff.gapfill import parse_pe_spans


def _pe_paf(tmp_path, name, pairs):
    """pairs: [(read, te1_on_L, ts2_on_R)]; flank len 1500, reads 150 bp."""
    lines = []
    for rn, te1, ts2 in pairs:
        # R1 '+' on g1|L ending at te1 (150 bp aln)
        lines.append(f"{rn}\t150\t0\t150\t+\tg1|L\t1500\t{te1 - 150}\t{te1}"
                     f"\t150\t150\t60\ttp:A:P")
        # R2 '-' on g1|R starting at ts2
        lines.append(f"{rn}\t150\t0\t150\t-\tg1|R\t1500\t{ts2}\t{ts2 + 150}"
                     f"\t150\t150\t60\ttp:A:P")
    p1 = tmp_path / name
    p1.write_text("\n".join(l for i, l in enumerate(lines) if i % 2 == 0) + "\n")
    p2 = tmp_path / (name + ".2")
    p2.write_text("\n".join(l for i, l in enumerate(lines) if i % 2 == 1) + "\n")
    return str(p1), str(p2)


def test_pe_span_estimates_size(tmp_path):
    # two consistent pairs: shortfall 50 each side, 150 bp alignments
    # portion = 150 + 50 = 200 per mate -> est = 500 - 400 = 100
    p1, p2 = _pe_paf(tmp_path, "a.paf",
                     [("r1", 1450, 50), ("r2", 1440, 55)])
    est = parse_pe_spans(p1, p2, insert=500)
    assert "g1" in est
    assert 95 <= est["g1"] <= 110


def test_pe_span_requires_two_pairs(tmp_path):
    p1, p2 = _pe_paf(tmp_path, "a.paf", [("r1", 1450, 50)])
    assert parse_pe_spans(p1, p2, insert=500) == {}


def test_pe_span_rejects_wrong_orientation(tmp_path):
    # both reads '+' -> not a proper FR span
    lines = [
        "r1\t150\t0\t150\t+\tg1|L\t1500\t1300\t1450\t150\t150\t60\ttp:A:P",
        "r1\t150\t0\t150\t+\tg1|R\t1500\t50\t200\t150\t150\t60\ttp:A:P",
        "r2\t150\t0\t150\t+\tg1|L\t1500\t1300\t1450\t150\t150\t60\ttp:A:P",
        "r2\t150\t0\t150\t+\tg1|R\t1500\t50\t200\t150\t150\t60\ttp:A:P",
    ]
    p1 = tmp_path / "a.paf"
    p1.write_text(lines[0] + "\n" + lines[2] + "\n")
    p2 = tmp_path / "a.paf.2"
    p2.write_text(lines[1] + "\n" + lines[3] + "\n")
    assert parse_pe_spans(str(p1), str(p2), insert=500) == {}


def test_pe_span_rejects_low_mapq_and_secondary(tmp_path):
    p1, p2 = _pe_paf(tmp_path, "a.paf",
                     [("r1", 1450, 50), ("r2", 1450, 50)])
    # strip the tp tag and drop mapq -> all rejected
    p1x = tmp_path / "b.paf"
    p1x.write_text(open(p1).read().replace("\ttp:A:P", "").replace("\t60\n",
                                                                   "\t5\n"))
    p2x = tmp_path / "b.paf.2"
    p2x.write_text(open(p2).read().replace("\ttp:A:P", "").replace("\t60\n",
                                                                   "\t5\n"))
    assert parse_pe_spans(str(p1x), str(p2x), insert=500) == {}
