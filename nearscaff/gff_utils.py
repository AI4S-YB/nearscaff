"""GFF3 parsing for CDS extraction and translation."""
from collections import defaultdict

# Standard genetic code
CODON_TABLE = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C','TGT':'C','TGA':'*','TGG':'W',
}

COMPLEMENT = {'A':'T','T':'A','C':'G','G':'C','N':'N'}


def _parse_gff_attrs(attr_str: str) -> dict:
    """Parse GFF3 attribute column (key=value;...) to dict."""
    result = {}
    for part in attr_str.split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            result[k.strip()] = v.strip()
    return result


def _reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq))


def extract_cds(gff3_text: str, fasta_text: str) -> dict:
    """Extract CDS sequences from GFF3 + FASTA.

    Args:
        gff3_text: GFF3 format string
        fasta_text: FASTA format string

    Returns:
        {mrna_id: concatenated_cds_sequence}
    """
    # Parse GFF3: collect CDS features by mRNA parent
    cds_by_mrna = defaultdict(list)
    for line in gff3_text.strip().split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) != 9:
            continue
        if fields[2] != 'CDS':
            continue
        attrs = _parse_gff_attrs(fields[8])
        parent = attrs.get('Parent', '')
        if parent:
            cds_by_mrna[parent].append({
                'chr': fields[0],
                'start': int(fields[3]) - 1,  # 0-based
                'end': int(fields[4]),         # exclusive
                'strand': fields[6],
                'phase': int(fields[7]),
            })

    # Parse FASTA
    seqs = {}
    current_id = None
    current_seq = []
    for line in fasta_text.strip().split('\n'):
        if line.startswith('>'):
            if current_id:
                seqs[current_id] = ''.join(current_seq).upper()
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_id:
        seqs[current_id] = ''.join(current_seq).upper()

    # Extract CDS for each mRNA
    result = {}
    for mrna_id, cds_list in cds_by_mrna.items():
        chr_name = cds_list[0]['chr']
        seq = seqs.get(chr_name, '')
        if not seq:
            continue
        cds_list.sort(key=lambda c: c['start'])
        cds_parts = []
        for cds in cds_list:
            cds_parts.append(seq[cds['start']:cds['end']])
        full_cds = ''.join(cds_parts)
        if cds_list[0]['strand'] == '-':
            full_cds = _reverse_complement(full_cds)
        result[mrna_id] = full_cds
    return result


def translate_cds(cds_seq: str) -> str:
    """Translate a CDS nucleotide sequence to protein."""
    aa = []
    for i in range(0, len(cds_seq) - 2, 3):
        codon = cds_seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)
