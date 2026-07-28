"""Tests for AGP v2.1 module."""
from nearscaff.agp import AGPReader, AGPWriter, AGPSeqLine, AGPGapLine


SAMPLE_AGP = (
    "scaf1\t1\t1000\t1\tW\tctg1\t1\t1000\t+\n"
    "scaf1\t1001\t1100\t2\tN\t100\tfragment\tyes\tmap\n"
    "scaf1\t1101\t2000\t3\tW\tctg2\t1\t900\t+\n"
)


def test_read_agp():
    reader = AGPReader()
    lines = reader.parse(SAMPLE_AGP)
    assert len(lines) == 3
    assert isinstance(lines[0], AGPSeqLine)
    assert lines[0].component_id == "ctg1"
    assert lines[0].orientation == "+"
    assert isinstance(lines[1], AGPGapLine)
    assert lines[1].gap_length == 100
    assert lines[1].gap_type == "fragment"
    assert isinstance(lines[2], AGPSeqLine)
    assert lines[2].component_id == "ctg2"


def test_write_agp():
    lines = [
        AGPSeqLine("scaf1", 1, 1000, 1, "W", "ctg1", 1, 1000, "+"),
        AGPGapLine("scaf1", 1001, 1100, 2, "N", 100, "fragment", "yes", "map"),
    ]
    writer = AGPWriter()
    output = writer.format(lines)
    expected = "scaf1\t1\t1000\t1\tW\tctg1\t1\t1000\t+\nscaf1\t1001\t1100\t2\tN\t100\tfragment\tyes\tmap\n"
    assert output == expected


def test_roundtrip():
    reader = AGPReader()
    writer = AGPWriter()
    original = reader.parse(SAMPLE_AGP)
    output = writer.format(original)
    parsed_again = reader.parse(output)
    assert len(parsed_again) == len(original)
    for a, b in zip(original, parsed_again):
        assert a.format() == b.format()


def test_parse_skips_comments():
    text = "# comment line\nscaf1\t1\t100\t1\tW\tctg1\t1\t100\t+\n"
    reader = AGPReader()
    lines = reader.parse(text)
    assert len(lines) == 1


def test_parse_skips_short_lines():
    text = "scaf1\t1\t100\nscaf1\t1\t100\t1\tW\tctg1\t1\t100\t+\n"
    reader = AGPReader()
    lines = reader.parse(text)
    assert len(lines) == 1


def test_all_component_types():
    """All valid sequence component types should parse."""
    for ct in ['A', 'D', 'F', 'G', 'O', 'P', 'W']:
        text = f"scaf\t1\t100\t1\t{ct}\tctg1\t1\t100\t+\n"
        lines = AGPReader().parse(text)
        assert len(lines) == 1
        assert lines[0].component_type == ct
