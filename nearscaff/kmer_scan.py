"""Jellyfish k-mer counting helpers for the phasing modules.

Vendored here so kmer_phasing/homeolog stay independent of the
scaffolding pipeline.  jellyfish count+dump outputs are cached next to
the input FASTA as ``<prefix>_<k>.fa``.
"""

import logging
import os
import subprocess

logger = logging.getLogger("nearscaff.kmer_scan")


def run_jellyfish_count(fasta_path: str, k: int, lower_count: int = 3,
                        hash_size: str = "100M", threads: int = 2) -> str:
    """Run jellyfish count+dump on a FASTA file, return path to dump output."""
    prefix = fasta_path.rsplit('.', 1)[0]
    jf_path = f"{prefix}_{k}.jf"
    dump_path = f"{prefix}_{k}.fa"

    if os.path.exists(dump_path):
        return dump_path

    cmd = ["jellyfish", "count", "-t", str(threads), "-m", str(k),
           "-s", hash_size, "--canonical", fasta_path, "-o", jf_path]
    subprocess.run(cmd, capture_output=True, text=True)

    cmd = ["jellyfish", "dump", "-c", "-o", dump_path, jf_path,
           "-L", str(lower_count)]
    subprocess.run(cmd, capture_output=True, text=True)

    if os.path.exists(jf_path):
        os.unlink(jf_path)
    return dump_path if os.path.exists(dump_path) else ""


def _contig_kmer_profile(seq: str, k: int, kmer_to_idx: dict,
                         n_features: int) -> list | None:
    """Count feature k-mers in a contig sequence.  Returns None if too short."""
    if len(seq) < k:
        return None
    counts = [0] * n_features
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        idx = kmer_to_idx.get(kmer)
        if idx is not None:
            counts[idx] += 1
    return counts


def scan_contig_profiles(
    contig_fasta: str, k: int, kmer_to_idx: dict, n_features: int
) -> tuple[list[str], list[list[int]]]:
    """Scan contigs and build k-mer count vectors.  Returns (names, profiles)."""
    contig_names = []
    profiles = []
    current_name = None
    current_seq = []

    with open(contig_fasta) as f:
        for line in f:
            if line.startswith('>'):
                if current_name and current_seq:
                    seq = ''.join(current_seq).upper()
                    vec = _contig_kmer_profile(seq, k, kmer_to_idx,
                                               n_features)
                    if vec is not None:
                        contig_names.append(current_name)
                        profiles.append(vec)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_name and current_seq:
            seq = ''.join(current_seq).upper()
            vec = _contig_kmer_profile(seq, k, kmer_to_idx, n_features)
            if vec is not None:
                contig_names.append(current_name)
                profiles.append(vec)

    return contig_names, profiles
