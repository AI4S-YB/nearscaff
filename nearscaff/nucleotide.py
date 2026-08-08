"""Constrained nucleotide alignment using minimap2 subprocess.

Unlike RagTag's full genome-to-genome alignment, minimap2 is run
only on regions specified by Block Tree LocalClusters.
"""
import subprocess, tempfile, os, logging, shutil

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


def _faidx_available(fasta_path: str) -> bool:
    """True if samtools is installed and <fasta>.fai exists."""
    return shutil.which("samtools") is not None and os.path.exists(fasta_path + ".fai")


def ensure_query_faid(query_fasta: str) -> bool:
    """Ensure a samtools .fai index exists for the query FASTA.

    Returns True if faidx is usable (index present or just created), else False.
    """
    if shutil.which("samtools") is None:
        return False
    if os.path.exists(query_fasta + ".fai"):
        return True
    result = subprocess.run(["samtools", "faidx", query_fasta],
                            capture_output=True, text=True)
    return result.returncode == 0


def _extract_contigs(fasta_path: str, contig_ids: list[str], output_path: str):
    """Extract specific contigs from a FASTA.

    Prefers samtools faidx (O(1) per contig); falls back to a linear scan when
    samtools/.fai is unavailable or faidx fails for any reason.
    """
    contig_set = set(contig_ids)
    if _faidx_available(fasta_path):
        try:
            result = subprocess.run(
                ["samtools", "faidx", fasta_path] + list(contig_ids),
                capture_output=True, text=True)
            if result.returncode == 0 and result.stdout:
                with open(output_path, "w") as f_out:
                    f_out.write(result.stdout)
                return
            logger.debug("samtools faidx extraction failed (%s); linear scan",
                         result.stderr.strip())
        except OSError as e:
            # E.g. Argument list too long (E2BIG) when extracting a very large
            # number of contigs in one invocation — fall back to streaming scan.
            logger.debug("samtools faidx extraction raised %s; linear scan", e)
    # Fallback: linear scan
    with open(fasta_path) as f_in, open(output_path, "w") as f_out:
        write = False
        for line in f_in:
            if line.startswith(">"):
                name = line[1:].split()[0]
                write = name in contig_set
            if write:
                f_out.write(line)


_ALIGN_CACHE_HEADER = ("contig", "chr", "r_start", "r_end", "strand",
                       "mapq", "hitlen", "identity")


def write_align_cache(path: str, cache: dict) -> None:
    """Write the per-contig best-alignment cache to a TSV.

    cache: {contig: {chr, r_start, r_end, strand, mapq, hitlen, identity}}
    """
    with open(path, "w") as f:
        f.write("\t".join(_ALIGN_CACHE_HEADER) + "\n")
        for contig in sorted(cache):
            e = cache[contig]
            f.write(f"{contig}\t{e['chr']}\t{e['r_start']}\t{e['r_end']}\t"
                    f"{e['strand']}\t{e['mapq']}\t{e['hitlen']}\t{e['identity']:.4f}\n")


def read_align_cache(path: str) -> dict:
    """Read the per-contig best-alignment cache; {} if missing/unreadable."""
    cache: dict[str, dict] = {}
    if not os.path.exists(path):
        return cache
    try:
        with open(path) as f:
            f.readline()  # header
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                cache[parts[0]] = {
                    "chr": parts[1],
                    "r_start": int(parts[2]),
                    "r_end": int(parts[3]),
                    "strand": parts[4],
                    "mapq": int(parts[5]),
                    "hitlen": int(parts[6]),
                    "identity": float(parts[7]),
                }
    except (ValueError, OSError):
        logger.warning("Align cache read failed at %s; ignoring", path)
        return {}
    return cache


def index_path_for(ref_fasta: str, preset: str, index_dir: str) -> str:
    """Stable, preset-specific path for a reusable minimap2 index."""
    os.makedirs(index_dir, exist_ok=True)
    return os.path.join(index_dir, f"{os.path.basename(ref_fasta)}.{preset}.mmi")


def build_ref_index(ref_fasta: str, preset: str, index_dir: str,
                    threads: int = 4) -> str:
    """Build (once) and return a reusable minimap2 reference index.

    Idempotent: reuses an existing preset-specific index. Raises RuntimeError
    if minimap2 index construction fails.
    """
    idx_path = index_path_for(ref_fasta, preset, index_dir)
    if os.path.exists(idx_path):
        # Reuse only if the index is not older than the reference FASTA;
        # a stale index from a swapped reference file would silently mis-align.
        try:
            if os.path.getmtime(idx_path) >= os.path.getmtime(ref_fasta):
                logger.info("  Reusing reference index %s", idx_path)
                return idx_path
            logger.info("  Reference index stale (ref newer); rebuilding")
        except OSError:
            pass  # mtime check failed — rebuild to be safe
    logger.info("  Building reference index (preset=%s) ...", preset)
    cmd = ["minimap2", "-t", str(threads), "-x", preset,
           "-d", idx_path, ref_fasta]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 index build failed: {result.stderr.strip()}")
    return idx_path


def _align_full_ref_cmd(ref_index: str, query_fa: str, preset: str,
                        secondary: int = 5, with_cigar: bool = False,
                        threads: int = 4) -> list[str]:
    """Build the minimap2 command used by align_to_full_reference.

    secondary > 0 -> -N (retain up to N secondary alignments, for extension
    recovery); secondary == 0 -> --secondary=no (refine only needs the best
    alignment, so skip the wasteful secondary computation).
    """
    cmd = ["minimap2", "-t", str(threads), "-x", preset]
    if with_cigar:
        cmd.append("-c")
    if secondary and secondary > 0:
        cmd += ["-N", str(secondary)]
    elif secondary == 0:
        cmd.append("--secondary=no")
    cmd += [ref_index, query_fa]
    return cmd


def align_to_full_reference(ref_index: str, query_fa: str, preset: str,
                            secondary: int = 5, with_cigar: bool = False,
                            threads: int = 4) -> str:
    """Align query contigs to a prebuilt reference index; return PAF text.

    secondary > 0: minimap2 -N (retain up to N secondary alignments).
    secondary == 0: disable secondary alignments (--secondary=no).
    with_cigar: include -c (CIGAR rescue) — slower; use only for precise passes.
    """
    cmd = _align_full_ref_cmd(ref_index, query_fa, preset,
                              secondary=secondary, with_cigar=with_cigar,
                              threads=threads)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed: {result.stderr.strip()}")
    return result.stdout
