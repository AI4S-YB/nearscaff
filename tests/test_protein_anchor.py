"""Tests for protein anchoring module."""
import os, tempfile
from nearscaff.protein_anchor import (
    extract_proteins_from_gff, parse_protein_paf
)

SIMPLE_GFF = """##gff-version 3
Chr1\t.\tgene\t1\t96\t.\t+\t.\tID=g1
Chr1\t.\tmRNA\t1\t96\t.\t+\t.\tID=g1.m1;Parent=g1
Chr1\t.\tCDS\t1\t96\t.\t+\t0\tID=g1.cds;Parent=g1.m1
"""


def test_extract_proteins_from_gff():
    """CDS extraction + translation produces valid protein FASTA."""
    # 96bp CDS = 32aa: "MKW" + 29 x "K"
    # ATG=M, AAA=K, TGG=W, then 87 A's = 29 AAA codons = 29K
    fa_text = ">Chr1\n" + "ATGAAATGG" + "A" * 87
    with tempfile.NamedTemporaryFile(suffix='.gff3', mode='w', delete=False) as gf, \
         tempfile.NamedTemporaryFile(suffix='.fa', mode='w', delete=False) as ff, \
         tempfile.NamedTemporaryFile(suffix='.faa', mode='w', delete=False) as pf:
        try:
            gf.write(SIMPLE_GFF); gf.flush(); gf.close()
            ff.write(fa_text); ff.flush(); ff.close()
            proteins = extract_proteins_from_gff(gf.name, ff.name, pf.name)
            assert 'g1.m1' in proteins
            assert proteins['g1.m1'] == "MKW" + "K" * 29
            # Check protein FASTA file
            with open(pf.name) as f:
                content = f.read()
            assert '>g1.m1' in content
            assert 'MKW' in content
        finally:
            for p in [gf.name, ff.name, pf.name]:
                if os.path.exists(p):
                    os.unlink(p)


def test_extract_proteins_short_orfs_filtered():
    """Proteins shorter than 30aa should be filtered out."""
    # Only 3bp CDS = 1aa, too short
    gff = """##gff-version 3
Chr1\t.\tgene\t1\t3\t.\t+\t.\tID=g_short
Chr1\t.\tmRNA\t1\t3\t.\t+\t.\tID=g_short.m1;Parent=g_short
Chr1\t.\tCDS\t1\t3\t.\t+\t0\tID=g_short.cds;Parent=g_short.m1
"""
    fa = ">Chr1\nATGNNN"
    with tempfile.NamedTemporaryFile(suffix='.gff3', mode='w', delete=False) as gf, \
         tempfile.NamedTemporaryFile(suffix='.fa', mode='w', delete=False) as ff, \
         tempfile.NamedTemporaryFile(suffix='.faa', mode='w', delete=False) as pf:
        try:
            gf.write(gff); gf.flush(); gf.close()
            ff.write(fa); ff.flush(); ff.close()
            proteins = extract_proteins_from_gff(gf.name, ff.name, pf.name)
            assert 'g_short.m1' not in proteins  # filtered
        finally:
            for p in [gf.name, ff.name, pf.name]:
                if os.path.exists(p):
                    os.unlink(p)


def test_parse_protein_paf():
    """Parse a simplified PAF line into GeneAnchor."""
    paf_line = (
        "g1.m1\t150\t10\t30\t+\tChr1\t300000\t100500\t100560\t"
        "50\t55\t60\tAS:i:120\tms:i:100\tnp:i:40\tfs:i:0\tst:i:0\t"
        "cg:Z:150M"
    )
    anchors = parse_protein_paf(paf_line, ref_genes={'g1.m1': 'AT1G01010'})
    assert len(anchors) == 1
    a = anchors[0]
    assert a.ref_gene == 'AT1G01010'
    assert a.strand == '+'
    assert a.score == 120
    assert a.n_frameshifts == 0
    assert a.n_stop_codons == 0
    assert a.query_gene == 'g1.m1'


def test_parse_protein_paf_no_ref_genes():
    """Without ref_genes mapping, use query name as ref_gene."""
    paf_line = (
        "AT1G01010\t150\t10\t30\t+\tChr1\t300000\t100500\t100560\t"
        "50\t55\t60\tAS:i:80\tms:i:90\tnp:i:30\tfs:i:1\tst:i:0\t"
        "cg:Z:150M"
    )
    anchors = parse_protein_paf(paf_line)
    assert len(anchors) == 1
    assert anchors[0].ref_gene == 'AT1G01010'
    assert anchors[0].n_frameshifts == 1
