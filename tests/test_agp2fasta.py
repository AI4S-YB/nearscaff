"""Tests for AGP -> scaffold FASTA conversion (nearscaff.agp2fasta)."""
from nearscaff.agp2fasta import agp_to_fasta, collect_placed_components


def _write(path, text):
    with open(path, 'w') as f:
        f.write(text)
    return str(path)


def _read_fasta(path):
    """Return [(full_header, sequence), ...] in file order."""
    records = []
    header = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(parts)))
                header = line.rstrip('\n')
                parts = []
            else:
                parts.append(line.strip())
    if header is not None:
        records.append((header, ''.join(parts)))
    return records


def test_agp_to_fasta_basic_and_unplaced(tmp_path):
    agp = _write(tmp_path / "out.agp",
                 "nearscaff_0001\t1\t8\t1\tW\tctg1\t1\t8\t+\n"
                 "nearscaff_0001\t9\t12\t2\tN\t4\tscaffold\tyes\tna\n"
                 "nearscaff_0001\t13\t16\t3\tW\tctg2\t1\t4\t-\n")
    # ctg2 is placed in reverse orientation: revcomp(GATC) = GATC -> use
    # an asymmetric sequence to make the check meaningful
    query = _write(tmp_path / "query.fa",
                   ">ctg1\nACGTACGT\n"
                   ">ctg2\nAAAC\n"
                   ">ctg3 unplaced contig description\nTTTTGG\n")

    out = tmp_path / "scaffolds.fa"
    n_scaf, n_unplaced = agp_to_fasta(agp, query, str(out))

    assert n_scaf == 1
    assert n_unplaced == 1

    records = _read_fasta(out)
    assert len(records) == 2

    # Scaffold: ctg1 (+) + 4 Ns + revcomp(AAAC) = GTTT
    assert records[0] == (">nearscaff_0001", "ACGTACGT" + "N" * 4 + "GTTT")

    # Unplaced contig appended with its full original header
    assert records[1] == (">ctg3 unplaced contig description", "TTTTGG")


def test_collect_placed_components(tmp_path):
    agp = _write(tmp_path / "out.agp",
                 "# comment\n"
                 "s1\t1\t10\t1\tW\tctgA\t1\t10\t+\n"
                 "s1\t11\t20\t2\tU\t10\tscaffold\tyes\tna\n"
                 "s2\t1\t5\t1\tW\tctgB\t1\t5\t-\n")
    assert collect_placed_components(agp) == {"ctgA", "ctgB"}


def test_agp_to_fasta_all_placed(tmp_path):
    agp = _write(tmp_path / "out.agp",
                 "s1\t1\t6\t1\tW\tctg1\t1\t6\t+\n")
    query = _write(tmp_path / "query.fa", ">ctg1\nAACCGG\n")
    out = tmp_path / "scaffolds.fa"
    n_scaf, n_unplaced = agp_to_fasta(agp, query, str(out))
    assert (n_scaf, n_unplaced) == (1, 0)
    assert _read_fasta(out) == [(">s1", "AACCGG")]
