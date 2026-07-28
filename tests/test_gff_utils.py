"""Tests for GFF3 utilities."""
from nearscaff.gff_utils import extract_cds, translate_cds

TINY_GFF = """##gff-version 3
Chr1\t.\tgene\t100\t300\t.\t+\t.\tID=g1
Chr1\t.\tmRNA\t100\t300\t.\t+\t.\tID=g1.m1;Parent=g1
Chr1\t.\tCDS\t101\t150\t.\t+\t0\tID=g1.cds1;Parent=g1.m1
Chr1\t.\tCDS\t251\t300\t.\t+\t0\tID=g1.cds2;Parent=g1.m1
"""

# 300bp: CDS1=[100..150)=50bp (GFF 1-based 101-150), CDS2=[250..300)=50bp (GFF 251-300), total 100bp
# Using all A's for simplicity
TINY_FA = ">Chr1\n" + "A" * 300


def test_extract_cds():
    cds_seqs = extract_cds(TINY_GFF, TINY_FA)
    assert "g1.m1" in cds_seqs
    # CDS1: 50bp + CDS2: 50bp = 100bp total
    assert len(cds_seqs["g1.m1"]) == 100


def test_translate_cds():
    # AAA = K (Lysine)
    seq = "AAAAAA"  # 6bp = 2 codons -> KK
    aa = translate_cds(seq)
    assert aa == "KK"


def test_translate_stop():
    # TAA = stop
    seq = "AAATAA"  # K*
    aa = translate_cds(seq)
    assert aa == "K*"
