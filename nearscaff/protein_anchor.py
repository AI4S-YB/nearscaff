"""miniprot wrapper for protein-to-genome alignment."""
import subprocess, os
from nearscaff.gff_utils import extract_cds, translate_cds
from nearscaff.types import GeneAnchor


def extract_proteins_from_gff(gff_path: str, fasta_path: str,
                               output_path: str) -> dict:
    """Extract CDS from GFF3, translate, write proteins FASTA.

    Returns: {mrna_id: protein_sequence}
    """
    with open(gff_path) as gf, open(fasta_path) as ff:
        gff_text = gf.read()
        fasta_text = ff.read()

    cds_seqs = extract_cds(gff_text, fasta_text)
    proteins = {}
    with open(output_path, 'w') as out:
        for mrna_id, cds in sorted(cds_seqs.items()):
            aa = translate_cds(cds)
            # Strip trailing stop codon for miniprot
            if aa.endswith('*'):
                aa = aa[:-1]
            if len(aa) >= 30:
                proteins[mrna_id] = aa
                out.write(f">{mrna_id}\n")
                for i in range(0, len(aa), 60):
                    out.write(aa[i:i+60] + "\n")
    return proteins


def build_protein_index(ref_fasta: str, index_path: str, threads: int = 4):
    """Build miniprot index for reference genome."""
    cmd = ["miniprot", "-t", str(threads), "-d", index_path, ref_fasta]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"miniprot index failed: {result.stderr}")


def run_protein_map(index_path: str, protein_fasta: str,
                     threads: int = 4, splice_model: int = 1) -> str:
    """Run miniprot mapping, return PAF output as string."""
    cmd = ["miniprot", "-t", str(threads), f"-j{splice_model}",
           index_path, protein_fasta]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"miniprot map failed: {result.stderr}")
    return result.stdout


def parse_protein_paf(paf_text: str, ref_genes: dict = None) -> list:
    """Parse miniprot PAF output to GeneAnchor list.

    PAF fields: qname qlen qstart qend strand tname tlen tstart tend
                matches blocklen mapq [tags...]

    Args:
        paf_text: PAF-format output from miniprot
        ref_genes: optional {mrna_id: gene_id} mapping

    Returns:
        List of GeneAnchor objects
    """
    if ref_genes is None:
        ref_genes = {}
    anchors = []
    for line in paf_text.strip().split('\n'):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 12:
            continue

        qname = fields[0]
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        strand = fields[4]
        tname = fields[5]
        tlen = int(fields[6])
        tstart = int(fields[7])
        tend = int(fields[8])
        matches = int(fields[9])
        blocklen = int(fields[10])
        mapq = int(fields[11])

        # Parse optional tags
        tags = {}
        for f in fields[12:]:
            parts = f.split(':')
            if len(parts) >= 3:
                tag_name = parts[0]
                tag_type = parts[1]
                tag_value = ':'.join(parts[2:])
                if tag_type == 'i':
                    tags[tag_name] = int(tag_value)
                elif tag_type == 'f':
                    tags[tag_name] = float(tag_value)
                elif tag_type == 'Z' or tag_type == 'z':
                    tags[tag_name] = tag_value

        ref_gene = ref_genes.get(qname, qname)
        identity = tags.get('np', 0) / max(matches, 1)

        anchors.append(GeneAnchor(
            query_contig=tname,
            query_gene=qname,
            ref_chr=tname,
            ref_gene=ref_gene,
            q_start=qstart,
            q_end=qend,
            r_start=tstart,
            r_end=tend,
            strand=strand,
            score=float(tags.get('AS', 0)),
            identity=identity,
            n_exons=1,
            n_frameshifts=tags.get('fs', 0),
            n_stop_codons=tags.get('st', 0),
        ))
    return anchors
