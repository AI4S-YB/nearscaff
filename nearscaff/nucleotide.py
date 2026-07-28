"""Constrained nucleotide alignment using minimap2 subprocess.

Unlike RagTag's full genome-to-genome alignment, minimap2 is run
only on regions specified by Block Tree LocalClusters.
"""
import subprocess, tempfile, os, logging

logger = logging.getLogger("nearscaff")


def run_constrained_alignment(ref_fasta: str, query_contigs_fasta: str,
                               region_chr: str, region_start: int,
                               region_end: int, margin: int = 50000,
                               preset: str = "asm5",
                               threads: int = 4) -> str:
    """Run minimap2 on a constrained region of the reference.

    Args:
        ref_fasta: path to reference genome FASTA
        query_contigs_fasta: path to query contigs FASTA (subset)
        region_chr: reference chromosome name
        region_start, region_end: region bounds (0-based)
        margin: padding added to both sides
        preset: minimap2 -x preset
        threads: number of threads

    Returns:
        PAF text output from minimap2
    """
    region_start = max(0, region_start - margin)
    region_end = region_end + margin
    region_str = f"{region_chr}:{region_start}-{region_end}"

    with tempfile.NamedTemporaryFile(suffix='.fa', mode='w', delete=False) as tf:
        ref_region_path = tf.name

    try:
        result = subprocess.run(
            ["samtools", "faidx", ref_fasta, region_str, "-o", ref_region_path],
            capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"samtools faidx failed: {result.stderr}")

        cmd = ["minimap2", "-t", str(threads), "-x", preset,
               "-c", ref_region_path, query_contigs_fasta]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {result.stderr}")
        return result.stdout
    finally:
        if os.path.exists(ref_region_path):
            os.unlink(ref_region_path)


def parse_nucleotide_paf(paf_text: str) -> list[dict]:
    """Parse minimap2 PAF output into lightweight dicts.

    Each dict has keys: query, qlen, qstart, qend, strand,
    ref_chr, rlen, rstart, rend, nmatch, hitlen, mapq.
    """
    entries = []
    for line in paf_text.strip().split('\n'):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 12:
            continue
        entries.append({
            'query': fields[0], 'qlen': int(fields[1]),
            'qstart': int(fields[2]), 'qend': int(fields[3]),
            'strand': fields[4],
            'ref_chr': fields[5], 'rlen': int(fields[6]),
            'rstart': int(fields[7]), 'rend': int(fields[8]),
            'nmatch': int(fields[9]), 'hitlen': int(fields[10]),
            'mapq': int(fields[11]),
        })
    return entries


def _extract_contigs(fasta_path: str, contig_ids: list[str], output_path: str):
    """Extract specific contigs from a FASTA file using samtools faidx."""
    contig_set = set(contig_ids)
    with open(fasta_path) as f_in, open(output_path, 'w') as f_out:
        write = False
        for line in f_in:
            if line.startswith('>'):
                name = line[1:].split()[0]
                write = name in contig_set
            if write:
                f_out.write(line)
