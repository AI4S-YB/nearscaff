"""Shared test fixtures for nearscaff tests."""
import os
import random
import pytest
import tempfile


def _random_cds(seed: int, length: int = 999) -> str:
    """Generate a pseudo-random CDS that miniprot can meaningfully align.

    Uses non-repetitive codons so each gene produces distinct, alignable
    protein sequence.  Seeded for reproducibility.
    """
    rng = random.Random(seed)
    bases = ["A", "C", "G", "T"]
    # Build sequence avoiding homopolymer runs > 4
    seq = []
    run = 1
    prev = ""
    while len(seq) < length:
        b = rng.choice(bases)
        if b == prev:
            run += 1
            if run > 4:
                continue
        else:
            run = 1
            prev = b
        seq.append(b)
    return "".join(seq)


@pytest.fixture
def tmp_fasta_dir():
    """Create temp directory that persists for the test duration."""
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def tiny_ref_fa():
    """Synthetic reference: one 100kb chromosome with 5 genes.

    Gene model: 999 bp CDS (333 aa), spaced across the chromosome.
    Each gene has a distinct coding sequence for unique miniprot anchoring.
    """
    # Gene positions: g1=10k, g2=30k, g3=60k, g4=75k, g5=90k
    gene_positions = [
        ("g1", 10000, "ATG" + _random_cds(1) + "TGA"),
        ("g2", 30000, "ATG" + _random_cds(2) + "TGA"),
        ("g3", 60000, "ATG" + _random_cds(3) + "TGA"),
        ("g4", 75000, "ATG" + _random_cds(4) + "TGA"),
        ("g5", 90000, "ATG" + _random_cds(5) + "TGA"),
    ]

    # Build chromosome: N-filled backbone with gene sequences patched in
    backbone = list("N" * 100000)
    for name, start, cds in gene_positions:
        for i, base in enumerate(cds):
            pos = start + i
            if pos < 100000:
                backbone[pos] = base

    seq = "".join(backbone)

    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        f.write(">Chr1\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")
        path = f.name

    yield path
    os.unlink(path)


@pytest.fixture
def tiny_gff3():
    """GFF3 with 5 genes matching tiny_ref_fa."""
    gff = """##gff-version 3
Chr1\t.\tgene\t10000\t10998\t.\t+\t.\tID=g1
Chr1\t.\tmRNA\t10000\t10998\t.\t+\t.\tID=g1.m1;Parent=g1
Chr1\t.\tCDS\t10000\t10998\t.\t+\t0\tID=g1.cds;Parent=g1.m1
Chr1\t.\tgene\t30000\t30998\t.\t+\t.\tID=g2
Chr1\t.\tmRNA\t30000\t30998\t.\t+\t.\tID=g2.m1;Parent=g2
Chr1\t.\tCDS\t30000\t30998\t.\t+\t0\tID=g2.cds;Parent=g2.m1
Chr1\t.\tgene\t60000\t60998\t.\t+\t.\tID=g3
Chr1\t.\tmRNA\t60000\t60998\t.\t+\t.\tID=g3.m1;Parent=g3
Chr1\t.\tCDS\t60000\t60998\t.\t+\t0\tID=g3.cds;Parent=g3.m1
Chr1\t.\tgene\t75000\t75998\t.\t+\t.\tID=g4
Chr1\t.\tmRNA\t75000\t75998\t.\t+\t.\tID=g4.m1;Parent=g4
Chr1\t.\tCDS\t75000\t75998\t.\t+\t0\tID=g4.cds;Parent=g4.m1
Chr1\t.\tgene\t90000\t90994\t.\t+\t.\tID=g5
Chr1\t.\tmRNA\t90000\t90994\t.\t+\t.\tID=g5.m1;Parent=g5
Chr1\t.\tCDS\t90000\t90994\t.\t+\t0\tID=g5.cds;Parent=g5.m1
"""
    with tempfile.NamedTemporaryFile(suffix=".gff3", mode="w", delete=False) as f:
        f.write(gff)
        path = f.name

    yield path
    os.unlink(path)


@pytest.fixture
def tiny_query_fa():
    """Query genome: 2 contigs from Chr1 with a gap.

    ctg1 = Chr1[0:40000] (contains g1 and g2)
    ctg2 = Chr1[50000:90000] (contains g3, g4)
    Purposeful 10kb gap between 40k-50k.
    """
    gene_positions = [
        ("g1", 10000, "ATG" + _random_cds(1) + "TGA"),
        ("g2", 30000, "ATG" + _random_cds(2) + "TGA"),
        ("g3", 60000, "ATG" + _random_cds(3) + "TGA"),
        ("g4", 75000, "ATG" + _random_cds(4) + "TGA"),
        ("g5", 90000, "ATG" + _random_cds(5) + "TGA"),
    ]

    backbone = list("N" * 100000)
    for name, start, cds in gene_positions:
        for i, base in enumerate(cds):
            pos = start + i
            if pos < 100000:
                backbone[pos] = base

    full_seq = "".join(backbone)
    ctg1_seq = full_seq[0:40000]
    ctg2_seq = full_seq[50000:90000]

    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        f.write(">ctg1\n")
        for i in range(0, len(ctg1_seq), 60):
            f.write(ctg1_seq[i:i + 60] + "\n")
        f.write(">ctg2\n")
        for i in range(0, len(ctg2_seq), 60):
            f.write(ctg2_seq[i:i + 60] + "\n")
        path = f.name

    yield path
    os.unlink(path)


@pytest.fixture
def tiny_block_tree():
    """A minimal Block Tree JSON for Stage 2 testing."""
    import json
    tree = {
        "node_id": "root",
        "level": "root",
        "ref_chr": None,
        "ref_start": None,
        "ref_end": None,
        "query_contigs": [],
        "orientation": "+",
        "synteny_score": 0.0,
        "confidence": 0.0,
        "flags": [],
        "gene_pair_count": 0,
        "children": [{
            "node_id": "SG_0",
            "level": "subgenome",
            "ref_chr": None,
            "ref_start": None,
            "ref_end": None,
            "query_contigs": ["ctg1", "ctg2"],
            "orientation": "+",
            "synteny_score": 200.0,
            "confidence": 1.0,
            "flags": [],
            "gene_pair_count": 20,
            "children": [{
                "node_id": "SG_0_Chr1",
                "level": "chromosome",
                "ref_chr": "Chr1",
                "ref_start": 5000,
                "ref_end": 95000,
                "query_contigs": ["ctg1", "ctg2"],
                "orientation": "+",
                "synteny_score": 180.0,
                "confidence": 0.9,
                "flags": [],
                "gene_pair_count": 16,
                "children": [{
                    "node_id": "SG_0_Chr1_lc0",
                    "level": "local",
                    "ref_chr": "Chr1",
                    "ref_start": 5000,
                    "ref_end": 45000,
                    "query_contigs": ["ctg1", "ctg2"],
                    "orientation": "+",
                    "synteny_score": 150.0,
                    "confidence": 0.85,
                    "flags": [],
                    "gene_pair_count": 12,
                    "children": [],
                }],
            }],
        }],
    }
    with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as f:
        json.dump(tree, f)
        path = f.name
    yield path
    os.unlink(path)
