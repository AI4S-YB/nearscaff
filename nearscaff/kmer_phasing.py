"""SubPhaser-style k-mer subgenome phasing.

STANDALONE MODULE: not part of the default nearscaff pipeline (Signal A is
skipped unless ``--kmer-phasing`` is given).  Run it directly via
``nearscaff kmer-phase`` — see docs/kmer_phasing.md for usage, expected
input quality (chromosome-scale or long-contig assemblies), and
limitations.

Internalizes the SubPhaser method (``Jellyfish.py::_filter_kmer`` and
``Cluster.py`` clustering/bootstrap) instead of wrapping it:

  1. Genome-wide jellyfish count (canonical) -> repetitive k-mer filter
     (total copies in [min_freq, max_freq], no top-N truncation).
  2. Per-contig count matrix, normalized by contig length.
  3. Fold-change filter: within each homology group, keep k-mers whose
     max/submax normalized frequency ratio >= min_fold in >= group_ratio
     of groups.  SubPhaser uses configured homoeologous chromosome groups;
     here Signal B homology groups (contigs clustered by reference
     coordinates) serve as the group proxy.  Without usable groups the
     filter degrades to the repetitive-k-mer filter only.
  4. Per-column Z-score -> KMeans(k) directly (no PCA).  When homology
     groups of exactly 2 members are available and k = 2, clustering is
     instead HG-aware: k-mer direction (which homeolog is enriched) must
     be consistent across groups — pair-orientations are estimated
     spectrally from the fold-direction matrix, globally directional
     k-mers are kept, and contigs are assigned by their summed signed
     Z-score (see ``hg_orientation_clustering``).  Plain KMeans splits
     outlier chromosomes instead of subgenomes on this signal.
  5. Bootstrap: resample k-mer columns with replacement, re-cluster,
     align clusters via Hungarian matching; per-contig confidence =
     fraction of replicates agreeing with the main clustering.
  6. Cluster IDs are arbitrary, so they are normalized to sg_0..sg_{k-1}
     by maximizing agreement with *seed_labels* (e.g. Signal C protein
     labels, permutation-enumerated for k <= 6); without seeds, clusters
     are named by total contig length (deterministic).
  7. Post-hoc top1-vs-top2 t-test per k-mer defines SG-specific k-mers
     (written to ``kmer_specific.tsv`` when *output_dir* is given).

sklearn/scipy are dev-only dependencies: they are imported lazily and
their absence degrades the signal to empty with a warning.
"""

import itertools
import logging
import os
import re
from collections import Counter, defaultdict

from nearscaff.config import PhasingConfig

logger = logging.getLogger("nearscaff.kmer_phasing")

_EPS = 1e-20


# ============================================================================
# Candidate k-mer selection (SubPhaser _filter_kmer, repetitive filter)
# ============================================================================

def _load_candidate_kmers(dump_path: str, min_freq: int,
                          max_freq: int) -> list[str]:
    """Load k-mers whose genome-wide copy number is in [min_freq, max_freq].

    This is SubPhaser's repetitive-k-mer filter: subgenome-differential
    signal lives in repeats (TEs), not in single-copy sequence.  No top-N
    truncation — the fold filter does the discriminative selection.
    """
    kmers = []
    with open(dump_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                count = int(parts[1])
                if min_freq <= count <= max_freq:
                    kmers.append(parts[0])
    return kmers


def _read_fasta_lengths(fasta_path: str) -> dict[str, int]:
    """Read sequence lengths from a FASTA file."""
    lengths = {}
    name = None
    total = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0]
                total = 0
            else:
                total += len(line.strip())
    if name is not None:
        lengths[name] = total
    return lengths


# ============================================================================
# Fold-change filter (SubPhaser _filter_kmer, differential filter)
# ============================================================================

def fold_filter_kmers(
    X_norm: "object",  # np.ndarray (n_contigs, n_kmers), length-normalized
    contig_names: list[str],
    homology_groups: dict[str, set[str]],
    min_fold: float = 2.0,
    group_ratio: float = 0.5,
    min_total_freq: int = 0,
    max_total_freq: int = 1_000_000_000,
    n_features: int | None = None,
    raw_col_sums: "object | None" = None,  # np.ndarray of raw per-kmer totals
) -> "object":  # np.ndarray of selected column indices
    """SubPhaser-style fold-change k-mer filtering.

    For each k-mer (column of X_norm):
      1. Optional raw abundance filter on [min_total_freq, max_total_freq].
      2. For each homology group (>= 2 members), compute
         fold = max / submax of normalized frequencies across members.
         Zero frequencies participate (a k-mer present in one homeologous
         copy and absent from the other is maximally differential); the
         denominator is floored at epsilon.
      3. Keep the k-mer if fold >= min_fold in >= group_ratio of groups.

    Returns indices of selected columns, ranked by the fraction of groups
    passed.  *n_features* caps the output (None = no truncation).
    """
    import numpy as np

    n_contigs, n_kmers = X_norm.shape
    name_to_idx = {name: i for i, name in enumerate(contig_names)}

    group_members: list[list[int]] = []
    for members in homology_groups.values():
        idxs = [name_to_idx[m] for m in members if m in name_to_idx]
        if len(idxs) >= 2:
            group_members.append(idxs)

    if len(group_members) < 2:
        logger.warning("  Fewer than 2 usable homology groups — "
                       "fold filter keeps all candidates")
        return np.arange(n_kmers, dtype=np.intp)

    n_groups = len(group_members)
    min_groups_pass = max(1, int(n_groups * group_ratio))

    selected = []
    scores = []
    for ki in range(n_kmers):
        if raw_col_sums is not None:
            total_freq = raw_col_sums[ki]
            if total_freq < min_total_freq or total_freq > max_total_freq:
                continue

        passes = 0
        for members in group_members:
            freqs = sorted((X_norm[mi, ki] for mi in members), reverse=True)
            if freqs[0] / (freqs[1] + _EPS) >= min_fold:
                passes += 1

        if passes >= min_groups_pass:
            selected.append(ki)
            scores.append(passes / n_groups)

    if n_features is not None and len(selected) > n_features:
        order = np.argsort(scores)[-n_features:]
        return np.array(selected, dtype=np.intp)[np.sort(order)]
    return np.array(selected, dtype=np.intp)


# ============================================================================
# Clustering + bootstrap (SubPhaser Cluster.py)
# ============================================================================

def zscore_columns(X: "object") -> "object":
    """Z-normalize each column (k-mer) so every k-mer weighs equally."""
    import numpy as np

    means = X.mean(axis=0)
    stds = X.std(axis=0)
    stds[stds == 0] = 1
    return (X - means) / stds


def _align_to_reference_labels(ref_labels: "object", labels: "object",
                               n_clusters: int) -> "object":
    """Map bootstrap cluster IDs onto reference cluster IDs.

    Uses Hungarian assignment on the cluster contingency table (falls back
    to greedy matching when scipy is unavailable).
    """
    import numpy as np

    contingency = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for r, l in zip(ref_labels, labels):
        contingency[r, l] += 1
    try:
        from scipy.optimize import linear_sum_assignment
        rows, cols = linear_sum_assignment(-contingency)
        mapping = {c: r for r, c in zip(rows, cols)}
    except ImportError:
        mapping = {}
        used = set()
        for r in range(n_clusters):
            for c in np.argsort(-contingency[r]):
                if c not in used:
                    mapping[int(c)] = r
                    used.add(int(c))
                    break
    return np.array([mapping.get(int(l), -1) for l in labels])


def cluster_with_bootstrap(X: "object", n_clusters: int,
                           n_replicates: int = 100,
                           random_state: int = 42
                           ) -> tuple["object", "object"]:
    """KMeans on Z-scored matrix + bootstrap confidence per sample.

    Each replicate resamples k-mer columns with replacement and
    re-clusters; confidence = fraction of replicates where the sample
    lands in the (Hungarian-aligned) same cluster as the main run.

    Returns (labels, confidences) as np.ndarray of int / float.
    """
    import numpy as np
    from sklearn.cluster import KMeans

    n_samples, n_kmers = X.shape
    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=random_state)
    labels = km.fit_predict(X)

    if n_replicates <= 0 or n_kmers < 2:
        return labels, np.ones(n_samples)

    rng = np.random.RandomState(random_state)
    support = np.zeros(n_samples)
    for _ in range(n_replicates):
        cols = rng.choice(n_kmers, size=n_kmers, replace=True)
        km_bs = KMeans(n_clusters=n_clusters, n_init=5,
                       random_state=rng.randint(0, 10000))
        bs_labels = km_bs.fit_predict(X[:, cols])
        aligned = _align_to_reference_labels(labels, bs_labels, n_clusters)
        support += (aligned == labels)

    return labels, support / n_replicates


# ============================================================================
# HG-aware orientation clustering (2 subgenomes, homeologous pairs)
# ============================================================================
#
# Plain KMeans on the Z-scored matrix degenerates on chromosome-scale
# inputs: most fold-selected k-mers are differential in only one homology
# group (pair-local repeats), so the dominant axis separates outlier
# chromosomes, not subgenomes (measured on a mixed fragmented two-species assembly:
# KMeans splits [14, 2], accuracy 0.625).  The subgenome signal lives in
# k-mers whose within-pair direction is consistent ACROSS pairs.  Since
# pair member order is arbitrary, the global orientation is unknown — it
# is estimated from the data:
#
#   1. Direction matrix D (k-mers x pairs): d = sign(x_a - x_b) where the
#      within-pair fold >= min_fold, else 0.
#   2. Orientation o in {+-1}^pairs: leading eigenvector of C = D_sub^T
#      D_sub over k-mers directional in >= 75% of pairs (pair-local noise
#      is incoherent across pairs and cancels; globally directional
#      k-mers dominate C), then a few EM refinement passes.
#   3. Keep globally directional k-mers (|D o| >= 0.6 * n_directional)
#      with direction m = sign(D o); contig score = sum of m_k * z_ck.
#   4. Bootstrap resamples k-mer columns and re-runs the whole estimate;
#      confidence = fraction of replicates agreeing with the main sign.


def _hg_pairs(homology_groups: dict[str, set[str]],
              contig_names: list[str]) -> list[tuple[int, int]]:
    """Index pairs from 2-member homology groups present in contig_names."""
    name_to_idx = {name: i for i, name in enumerate(contig_names)}
    pairs = []
    for members in homology_groups.values():
        idxs = [name_to_idx[m] for m in members if m in name_to_idx]
        if len(idxs) == 2 and len(members) == 2:
            pairs.append((idxs[0], idxs[1]))
    return pairs


def _direction_matrix(X: "object", pairs: list[tuple[int, int]],
                      min_fold: float) -> "object":
    """D (k-mers x pairs): sign(x_a - x_b) where within-pair fold passes."""
    import numpy as np

    A = X[[a for a, _ in pairs]].T
    B = X[[b for _, b in pairs]].T
    fold = (np.maximum(A, B) + _EPS) / (np.minimum(A, B) + _EPS)
    return np.where(fold >= min_fold, np.sign(A - B), 0.0)


def _estimate_orientation(D: "object") -> "object | None":
    """Pair orientations o in {+-1}^P from the direction matrix.

    Spectral init on k-mers directional in >= 75% of pairs, then EM
    refinement.  Returns None when too few k-mers are usable.
    """
    import numpy as np

    n_kmers, n_pairs = D.shape
    nnz = (D != 0).sum(axis=1)
    # Adaptive subset: k-mers directional in >= 75% of pairs is the ideal;
    # with many sparse pairs (small/fragmented units) almost no k-mer
    # reaches that, so fall back to the top-5% most directional k-mers.
    thr = max(2, int(np.ceil(0.75 * n_pairs)))
    if (nnz >= thr).sum() < 100:
        thr = max(2, int(np.percentile(nnz, 95)))
    sub = D[nnz >= thr]
    if len(sub) < 100:
        sub = D[nnz >= nnz.max()]
        if len(sub) < 100:
            return None
    C = sub.T @ sub
    _w, V = np.linalg.eigh(C)
    o = np.sign(V[:, -1])
    o[o == 0] = 1
    for _ in range(50):
        m = np.sign(D @ o)
        o_new = np.sign(m @ D)
        o_new[o_new == 0] = o[o_new == 0]
        if (o_new == o).all():
            break
        o = o_new
    return o


def _orientation_scores(X_z: "object", D: "object",
                        o: "object") -> tuple["object", "object"]:
    """(scores, kept k-mer mask) given pair orientations."""
    import numpy as np

    nnz = (D != 0).sum(axis=1)
    agree = D @ o
    m = np.sign(agree)
    keep = np.abs(agree) >= 0.6 * np.maximum(nnz, 1)
    if keep.sum() < 10:
        return None, keep
    return X_z[:, keep] @ m[keep], keep


def hg_orientation_clustering(
    X: "object",  # np.ndarray (n_contigs, n_kmers), length-normalized
    contig_names: list[str],
    homology_groups: dict[str, set[str]],
    min_fold: float = 2.0,
    n_bootstrap: int = 100,
    random_state: int = 42,
) -> "tuple[object, object] | None":
    """HG-aware 2-group clustering via cross-pair direction consistency.

    Returns (labels, confidences) like ``cluster_with_bootstrap``, or
    None when the homology groups are unusable (< 2 two-member pairs, or
    too few directional k-mers) so the caller can fall back to KMeans.
    """
    import numpy as np

    pairs = _hg_pairs(homology_groups, contig_names)
    if len(pairs) < 2:
        return None

    D = _direction_matrix(X, pairs, min_fold)
    o = _estimate_orientation(D)
    if o is None:
        logger.info("  HG orientation: too few directional k-mers — "
                    "falling back to KMeans")
        return None

    X_z = zscore_columns(X)
    scores, keep = _orientation_scores(X_z, D, o)
    if scores is None:
        logger.info("  HG orientation: too few consistent k-mers — "
                    "falling back to KMeans")
        return None
    labels = (scores < 0).astype(np.int64)
    logger.info("  HG orientation: %d pairs, %d/%d directional k-mers, "
                "group sizes %s", len(pairs), int(keep.sum()), D.shape[0],
                sorted(np.bincount(labels, minlength=2).tolist(),
                       reverse=True))

    if n_bootstrap <= 0 or X.shape[1] < 2:
        return labels, np.ones(len(contig_names))

    rng = np.random.RandomState(random_state)
    support = np.zeros(len(contig_names))
    n_kmers = X.shape[1]
    for _ in range(n_bootstrap):
        cols = rng.choice(n_kmers, size=n_kmers, replace=True)
        D_bs = D[cols]
        o_bs = _estimate_orientation(D_bs)
        if o_bs is None:
            continue
        s_bs, _keep_bs = _orientation_scores(X_z[:, cols], D_bs, o_bs)
        if s_bs is None:
            continue
        flip = np.sign(np.dot(s_bs, scores)) or 1.0
        support += (np.sign(s_bs * flip) == np.sign(scores))
    return labels, support / n_bootstrap


def enforce_hg_constraints(labels: "object", confidence: "object",
                           contig_names: list[str],
                           homology_groups: dict[str, set[str]],
                           ) -> tuple["object", "object"]:
    """Force the two members of each homology pair into different groups.

    When both members land on the same side, the one with the lower
    confidence is flipped to the other group and its confidence is capped
    at 0.5 (i.e. flagged unreliable, unphased at the default threshold).
    """
    import numpy as np

    labels = labels.copy()
    confidence = confidence.copy()
    n_flipped = 0
    for a, b in _hg_pairs(homology_groups, contig_names):
        if labels[a] == labels[b]:
            weak = a if confidence[a] <= confidence[b] else b
            labels[weak] = 1 - labels[weak]
            confidence[weak] = min(confidence[weak], 0.5)
            n_flipped += 1
    if n_flipped:
        logger.info("  HG constraint: %d pair member(s) forcibly split "
                    "(confidence capped at 0.5)", n_flipped)
    return labels, confidence


# ============================================================================
# Label normalization (cluster ID -> stable sg_N names)
# ============================================================================


def _canon_sg(label: str) -> str:
    """Canonicalize a subgenome label to 0-indexed ``sg_N`` form."""
    label = str(label).strip()
    m = re.match(r'^[Ss][Gg]_(\d+)$', label)
    if m:
        return f"sg_{int(m.group(1))}"
    m = re.match(r'^SG0*(\d+)$', label)
    if m:
        return f"sg_{int(m.group(1)) - 1}"
    m = re.match(r'^[Ss][Gg](\d+)$', label)
    if m:
        return f"sg_{int(m.group(1))}"
    return label


def best_label_permutation(source: dict[str, str],
                           reference: dict[str, str]) -> dict[str, str]:
    """Find the source->reference label mapping maximizing agreement.

    Enumerates all permutations when both label sets have <= 6 members,
    otherwise falls back to greedy matching.  Only shared keys vote.
    """
    shared = [key for key in source if key in reference]
    if not shared:
        return {}

    src_labels = sorted({source[k] for k in shared})
    ref_labels = sorted({reference[k] for k in shared})

    if len(src_labels) == len(ref_labels) and len(src_labels) <= 6:
        best_map: dict[str, str] = {}
        best_score = -1
        for perm in itertools.permutations(ref_labels):
            mapping = dict(zip(src_labels, perm))
            score = sum(1 for k in shared
                        if mapping[source[k]] == reference[k])
            if score > best_score:
                best_score, best_map = score, mapping
        return best_map

    # Greedy fallback for larger / unequal label sets
    overlap = defaultdict(Counter)
    for key in shared:
        overlap[reference[key]][source[key]] += 1
    mapping = {}
    used = set()
    for ref_label in sorted(overlap, key=lambda r: -sum(overlap[r].values())):
        for src_label, count in overlap[ref_label].most_common():
            if src_label not in used:
                mapping[src_label] = ref_label
                used.add(src_label)
                break
    return mapping


def normalize_cluster_labels(labels: "object", contig_names: list[str],
                             contig_lengths: dict[str, int],
                             seed_labels: dict[str, str] | None,
                             n_subgenomes: int) -> list[str]:
    """Map arbitrary KMeans cluster IDs to stable sg_0..sg_{k-1} names.

    With *seed_labels* (e.g. Signal C protein labels) the mapping is the
    agreement-maximizing permutation; without seeds, clusters are ranked
    by total contig length (descending) for determinism.
    """
    cluster_ids = sorted(set(int(l) for l in labels))
    cluster_of = {name: int(labels[i]) for i, name in enumerate(contig_names)}

    if seed_labels:
        seeds = {ctg: _canon_sg(sg) for ctg, sg in seed_labels.items()
                 if ctg in cluster_of}
        clusters = {ctg: f"sg_{cluster_of[ctg]}" for ctg in seeds}
        if len(seeds) >= 10:
            mapping = best_label_permutation(clusters, seeds)
            n_agree = sum(1 for ctg in seeds
                          if mapping.get(clusters[ctg]) == seeds[ctg])
            logger.info("  Cluster naming: %d/%d seed labels agree after "
                        "permutation alignment", n_agree, len(seeds))
            return [mapping.get(f"sg_{cluster_of[name]}",
                                f"sg_{cluster_of[name]}")
                    for name in contig_names]
        logger.warning("  Too few seed labels (%d < 10) — naming clusters "
                       "by size instead", len(seeds))

    # No usable seeds: rank clusters by total contig length, descending
    cluster_len = Counter()
    for name in contig_names:
        cluster_len[cluster_of[name]] += contig_lengths.get(name, 0)
    order = sorted(cluster_ids, key=lambda c: (-cluster_len[c], c))
    mapping = {c: f"sg_{rank}" for rank, c in enumerate(order)}
    return [mapping[cluster_of[name]] for name in contig_names]


# ============================================================================
# SG-specific k-mers (SubPhaser Cluster.output_kmers)
# ============================================================================

def find_specific_kmers(X: "object", kmer_names: list[str],
                        labels: list[str], pval: float = 0.05
                        ) -> list[tuple[str, str, float, list[float]]]:
    """Per k-mer, t-test of the top-1 vs top-2 subgenome group means.

    Returns [(kmer, sg_label, pvalue, group_means)] for k-mers with
    p < *pval* — the SG-specific k-mer set (SubPhaser step 6).
    """
    import numpy as np
    from scipy import stats

    sg_labels = sorted(set(labels))
    groups = {sg: [i for i, l in enumerate(labels) if l == sg]
              for sg in sg_labels}

    specific = []
    for ki, kmer in enumerate(kmer_names):
        grouped = [X[groups[sg], ki] for sg in sg_labels]
        means = [float(g.mean()) if len(g) else 0.0 for g in grouped]
        order = sorted(range(len(sg_labels)), key=lambda i: -means[i])
        if len(order) < 2:
            continue
        top1, top2 = grouped[order[0]], grouped[order[1]]
        if means[order[0]] <= 0:
            continue
        test = stats.ttest_ind(top1, top2, equal_var=False)
        if test.pvalue < pval:
            specific.append((kmer, sg_labels[order[0]],
                             float(test.pvalue), means))
    return specific


def _write_specific_kmers(specific: list[tuple[str, str, float, list]],
                          out_path: str):
    """Write the SG-specific k-mer table (TSV)."""
    with open(out_path, 'w') as f:
        f.write("kmer\tsubgenome\tp_value\tgroup_means\n")
        for kmer, sg, p, means in specific:
            f.write(f"{kmer}\t{sg}\t{p:.3g}\t"
                    f"{','.join(f'{m:.6g}' for m in means)}\n")


# ============================================================================
# Main entry point
# ============================================================================

def run_subgenome_kmer_phasing(
    contig_fasta: str,
    *,
    k: int = 15,
    n_subgenomes: int = 2,
    ncpu: int = 4,
    homology_groups: dict[str, set[str]] | None = None,
    contig_lengths: dict[str, int] | None = None,
    seed_labels: dict[str, str] | None = None,
    cfg: PhasingConfig | None = None,
    output_dir: str | None = None,
) -> dict[str, tuple[str, float]]:
    """SubPhaser-style unsupervised k-mer subgenome phasing.

    Returns ``{contig: (sg_label, confidence_0_to_1)}``; contigs below
    ``cfg.bootstrap_threshold`` confidence are omitted (unphased).
    """
    cfg = cfg or PhasingConfig()

    try:
        import numpy as np
        from sklearn.cluster import KMeans  # noqa: F401
    except ImportError:
        logger.warning("  sklearn/numpy not installed — k-mer phasing "
                       "unavailable (dev dependency)")
        return {}

    # Lazy import: subgenome imports this module from run_kmer_phasing
    from nearscaff.kmer_scan import run_jellyfish_count as _run_jellyfish_count, scan_contig_profiles as _scan_contig_profiles

    logger.info("  SubPhaser-style k-mer phasing (k=%d, min_freq=%d, "
                "min_fold=%.1f)", k, cfg.min_freq, cfg.min_fold)

    # Step 1: genome-wide k-mer count
    dump_path = _run_jellyfish_count(contig_fasta, k=k, lower_count=3,
                                     threads=ncpu)
    if not dump_path:
        logger.warning("  jellyfish failed — k-mer signal absent")
        return {}

    # Step 2: repetitive-k-mer candidates (no top-N truncation)
    candidates = _load_candidate_kmers(dump_path, cfg.min_freq, cfg.max_freq)
    if len(candidates) < 10:
        logger.warning("  Too few candidate k-mers (%d) — k-mer signal absent",
                       len(candidates))
        return {}
    logger.info("  %d candidate k-mers (copies in [%d, %d])",
                len(candidates), cfg.min_freq, cfg.max_freq)

    # Step 3: per-contig count matrix, length-normalized
    kmer_to_idx = {km: i for i, km in enumerate(candidates)}
    contig_names, raw_profiles = _scan_contig_profiles(
        contig_fasta, k, kmer_to_idx, len(candidates))

    if contig_lengths is None:
        contig_lengths = _read_fasta_lengths(contig_fasta)

    keep = [i for i, name in enumerate(contig_names)
            if contig_lengths.get(name, 0) >= cfg.min_contig_len]
    n_skipped = len(contig_names) - len(keep)
    if n_skipped:
        logger.info("  %d contigs below min_contig_len=%d skipped",
                    n_skipped, cfg.min_contig_len)
    contig_names = [contig_names[i] for i in keep]
    if len(contig_names) < n_subgenomes:
        logger.warning("  %d contigs < %d subgenomes — k-mer signal absent",
                       len(contig_names), n_subgenomes)
        return {}

    X_raw = np.array([raw_profiles[i] for i in keep], dtype=np.float64)
    lengths = np.array([max(contig_lengths.get(name, k), k)
                        for name in contig_names], dtype=np.float64)
    X_norm = X_raw / lengths[:, np.newaxis]
    logger.info("  %d contigs profiled", len(contig_names))

    # Step 4: fold-change filter (degrades to repetitive-only without groups)
    hg_contigs: set[str] = set()
    for members in (homology_groups or {}).values():
        hg_contigs.update(members)
    hg_coverage = len(hg_contigs & set(contig_names)) / max(len(contig_names), 1)

    if homology_groups and len(homology_groups) >= 2 and hg_coverage >= 0.1:
        logger.info("  Fold-change filtering with %d homology groups "
                    "(%.1f%% contig coverage) ...",
                    len(homology_groups), hg_coverage * 100)
        selected = fold_filter_kmers(
            X_norm, contig_names, homology_groups,
            min_fold=cfg.min_fold, group_ratio=cfg.group_ratio)
        if len(selected) < 10:
            logger.warning("  Fold filter kept only %d k-mers — using all "
                           "%d repetitive candidates", len(selected),
                           X_norm.shape[1])
            selected = np.arange(X_norm.shape[1], dtype=np.intp)
    else:
        if homology_groups:
            logger.info("  Homology group coverage too low (%.1f%%) — "
                        "repetitive k-mers only", hg_coverage * 100)
        else:
            logger.info("  No homology groups — repetitive k-mers only")
        selected = np.arange(X_norm.shape[1], dtype=np.intp)

    logger.info("  %d k-mers selected for clustering", len(selected))
    selected_kmers = [candidates[i] for i in selected]

    # Step 5: per-column Z-score, then KMeans directly (no PCA)
    X_z = zscore_columns(X_norm[:, selected])

    # Step 6: clustering + bootstrap confidence.  With 2-member homology
    # groups and k = 2, use HG-aware orientation clustering (cross-pair
    # direction consistency); plain KMeans degenerates to splitting
    # outlier chromosomes on this signal.  Fall back to KMeans when the
    # groups are unusable.
    labels = confidence = None
    if n_subgenomes == 2 and homology_groups:
        hg_res = hg_orientation_clustering(
            X_norm[:, selected], contig_names, homology_groups,
            min_fold=cfg.min_fold, n_bootstrap=cfg.bootstrap_reps)
        if hg_res is not None:
            labels, confidence = hg_res
            labels, confidence = enforce_hg_constraints(
                labels, confidence, contig_names, homology_groups)
    if labels is None:
        labels, confidence = cluster_with_bootstrap(
            X_z, n_subgenomes, n_replicates=cfg.bootstrap_reps)

    # Step 6b: EM LLR refinement on copy-band k-mer composition.
    # For closely related subgenomes the repetitive spectrum has not
    # diverged and the KMeans above degenerates; k-mer LLR composition
    # still separates (see module notes).  KMeans labels seed the EM; when
    # the KMeans init is collapsed, protein seed labels are used instead
    # if they cover both groups well.  A collapsed init without usable
    # seeds makes the EM skip itself (returns {} -> KMeans kept).
    if cfg.em_iterations > 0:
        init = {name: int(labels[i]) for i, name in enumerate(contig_names)}
        min_needed = max(50, int(0.15 * len(init)))
        cluster_sizes = Counter(init.values())
        collapsed = any(cluster_sizes.get(g, 0) < min_needed
                        for g in range(n_subgenomes))
        init_source = "KMeans clusters"
        if collapsed and seed_labels:
            seed_init = _seed_labels_as_int(seed_labels, contig_names)
            if seed_init is not None:
                init = seed_init
                init_source = "protein seed labels"
        logger.info("  EM init: %s (group sizes %s)", init_source,
                    sorted(Counter(init.values()).values(), reverse=True))
        em_result = em_marker_phasing(
            contig_fasta, k=cfg.em_k, n_subgenomes=n_subgenomes,
            init_labels=init, contig_lengths=contig_lengths,
            ncpu=ncpu, cfg=cfg)
        if em_result:
            return em_result
        logger.warning("  EM refinement unavailable — keeping KMeans result")

    # Step 7: normalize cluster IDs to stable sg_N names
    sg_names = normalize_cluster_labels(labels, contig_names, contig_lengths,
                                        seed_labels, n_subgenomes)

    # Step 8: SG-specific k-mers (t-test, SubPhaser output_kmers)
    try:
        specific = find_specific_kmers(X_norm[:, selected], selected_kmers,
                                       sg_names, pval=cfg.specific_pval)
        logger.info("  %d SG-specific k-mers (p < %g)",
                    len(specific), cfg.specific_pval)
        if output_dir:
            out_path = os.path.join(output_dir, "kmer_specific.tsv")
            _write_specific_kmers(specific, out_path)
            logger.info("  SG-specific k-mers written to %s", out_path)
    except ImportError:
        logger.warning("  scipy not installed — SG-specific k-mer "
                       "test skipped")

    result = {}
    for i, name in enumerate(contig_names):
        if confidence[i] >= cfg.bootstrap_threshold:
            result[name] = (sg_names[i], float(confidence[i]))
    n_unphased = len(contig_names) - len(result)
    logger.info("  %d contigs phased (confidence >= %.2f), %d unphased",
                len(result), cfg.bootstrap_threshold, n_unphased)
    return result


# ============================================================================
# EM LLR phasing (log-likelihood-ratio k-mer scoring)
# ============================================================================
#
# For closely related subgenomes the repetitive-k-mer spectrum has not
# diverged yet, so the SubPhaser frequency signal collapses (measured on
# On closely related two-species mixes: clustering
# degenerates to one cluster, ~60% accuracy).  But fixed differences make
# copy-band k-mers group-informative, and every contig carries hundreds of
# them.  EM alternates:
#
#   M-step: from the current contig assignments, accumulate per-group k-mer
#           copy counts (np.bincount over one cached scan); each k-mer gets
#           an LLR weight  w = log((c0+2)/size0) - log((c1+2)/size1)
#           (frequencies normalized by group total hits, smoothed by ~2
#           copies), zeroed when total copies < em_min_total.
#   E-step: score each contig by the summed weight of its k-mer hits and
#           reassign by sign.
#
# The fasta is scanned ONCE into (hit_contig, hit_kmer) arrays; all passes
# iterate on them in memory.  Design findings (synthetic benchmarks +
# oracle experiments on the real mixes):
#   * LLR scoring strictly beats the earlier binary "marker" scheme
#     (k-mer is a group marker if >= em_marker_frac of its size-normalized
#     copies sit in one group): 94.8% vs 88.6% oracle accuracy.
#   * A single M->E pass is the sweet spot; further iterations degrade good
#     initializations and degenerate on bad ones.
#   * The full copy band [3, 200000] without subsampling beats the old
#     [3, 50] band (94.8% vs 92.3% oracle).
#   * A collapsed initialization (one group nearly empty) makes EM skip
#     entirely — running it can only produce garbage and wastes ~20 min.


def _em_base_tables():
    """2-bit base encoding table (255 = invalid) + complement table."""
    import numpy as np
    code = np.full(256, 255, dtype=np.uint8)
    code[ord('A')] = 0
    code[ord('C')] = 1
    code[ord('G')] = 2
    code[ord('T')] = 3
    comp = np.array([3, 2, 1, 0], dtype=np.uint8)
    return code, comp


def _em_canonical_codes(seq: bytes, k: int, stride: int):
    """int64 canonical 2-bit k-mer codes at every stride-th position.

    Windows containing non-ACGT bases are coded as -1.
    """
    import numpy as np
    base_code, comp = _em_base_tables()
    b = base_code[np.frombuffer(seq, dtype=np.uint8)]
    n = len(b)
    if n < k:
        return np.empty(0, dtype=np.int64)
    pos = np.arange(0, n - k + 1, stride, dtype=np.int64)
    bad = (b == 255).astype(np.int64)
    cs = np.concatenate(([0], np.cumsum(bad)))
    valid = (cs[pos + k] - cs[pos]) == 0
    b = np.minimum(b, 3)  # keep indexing safe; masked via `valid` below
    fwd = np.zeros(len(pos), dtype=np.int64)
    rc = np.zeros(len(pos), dtype=np.int64)
    for o in range(k):
        fwd = (fwd << 2) | b[pos + o]
    for o in range(k):
        rc |= comp[b[pos + o]].astype(np.int64) << np.int64(2 * o)
    code = np.minimum(fwd, rc)
    code[~valid] = -1
    return code


def _em_load_candidates(dump_path: str, min_copy: int, max_copy: int,
                        max_candidates: int, seed: int = 42):
    """Jellyfish dump -> (codes, counts) for the copy band.

    codes: sorted unique int64 array of 2-bit k-mer codes; counts: int32
    total copies per code (canonical-mate duplicates merged).  When the
    band exceeds *max_candidates*, a deterministic seeded subsample is
    taken.
    """
    import numpy as np

    enc = {"A": 0, "C": 1, "G": 2, "T": 3}
    codes = []
    counts = []
    with open(dump_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            c = int(parts[1])
            if not (min_copy <= c <= max_copy):
                continue
            v = 0
            for ch in parts[0]:
                v = (v << 2) | enc[ch]
            codes.append(v)
            counts.append(c)
    if not codes:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int32)
    codes = np.array(codes, dtype=np.int64)
    counts = np.array(counts, dtype=np.int64)
    if max_candidates and len(codes) > max_candidates:
        rng = np.random.RandomState(seed)
        keep = np.sort(rng.choice(len(codes), max_candidates, replace=False))
        codes, counts = codes[keep], counts[keep]
        logger.info("  EM: subsampled candidates to %d (seed=%d)",
                    len(codes), seed)
    order = np.argsort(codes, kind="stable")
    codes, counts = codes[order], counts[order]
    uniq, starts = np.unique(codes, return_index=True)
    merged = np.add.reduceat(counts, starts)
    return uniq, merged.astype(np.int32)


def _em_scan(contig_fasta: str, k: int, cand, stride: int):
    """Single vectorized pass over the contigs.

    Returns (hit_contig, hit_kmer, names, lengths): parallel int32 arrays
    mapping each candidate k-mer occurrence to its contig index and
    candidate index, plus the contig names and lengths in scan order.
    """
    import numpy as np

    n_cand = len(cand)
    names: list[str] = []
    lengths: list[int] = []
    hit_c: list = []
    hit_j: list = []

    def flush(nm, seq):
        ci = len(names)
        names.append(nm)
        lengths.append(len(seq))
        code = _em_canonical_codes(seq.encode(), k, stride)
        code = code[code >= 0]
        if not len(code):
            return
        at = np.clip(np.searchsorted(cand, code), 0, n_cand - 1)
        m = cand[at] == code
        j = at[m]
        hit_c.append(np.full(len(j), ci, dtype=np.int32))
        hit_j.append(j.astype(np.int32))

    name = None
    parts = []
    with open(contig_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    flush(name, "".join(parts).upper())
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        flush(name, "".join(parts).upper())

    hc = np.concatenate(hit_c) if hit_c else np.empty(0, dtype=np.int32)
    hj = np.concatenate(hit_j) if hit_j else np.empty(0, dtype=np.int32)
    return hc, hj, names, np.array(lengths, dtype=np.int64)


def _seed_labels_as_int(seed_labels: dict[str, str],
                        contig_names: list[str],
                        min_per_group: int = 50) -> dict[str, int] | None:
    """Map sg-style seed labels to 0/1 ints for *contig_names*.

    Returns None unless exactly 2 seed groups each cover >= min_per_group
    of the contigs.
    """
    names = set(contig_names)
    pairs = [(c, _canon_sg(sg)) for c, sg in seed_labels.items()
             if c in names]
    uniq = sorted({sg for _c, sg in pairs})
    if len(uniq) != 2:
        return None
    lmap = {sg: i for i, sg in enumerate(uniq)}
    out = {c: lmap[sg] for c, sg in pairs}
    sizes = Counter(out.values())
    if any(sizes.get(g, 0) < min_per_group for g in range(2)):
        return None
    return out


def em_marker_phasing(contig_fasta: str, *, k: int = 21, n_subgenomes: int = 2,
                      init_labels: dict[str, int] | None = None,
                      contig_lengths: dict[str, int] | None = None,
                      ncpu: int = 4, cfg: PhasingConfig | None = None,
                      seed: int = 42) -> dict[str, tuple[str, float]]:
    """EM LLR phasing on copy-band k-mer composition.  See module notes.

    *init_labels* maps contig -> integer cluster (e.g. from the KMeans
    stage).  EM requires a usable init: without labels, or with a
    nearly-empty group (< 50 contigs on either side), it is skipped
    entirely and {} is returned so the caller keeps its own (KMeans)
    result.  After refinement, results agreeing with the init on < 60%
    of shared contigs are discarded as well (a bad init converges to
    noise, and this catches it).

    Returns {contig: (sg_label, confidence)}; contigs with fewer than
    ``cfg.em_min_markers`` informative (nonzero-weight) k-mer hits are
    omitted (unphased).  Confidence = |score| / sum(|hit weights|), 0..1.
    """
    import time

    import numpy as np

    cfg = cfg or PhasingConfig()
    from nearscaff.kmer_scan import run_jellyfish_count as _run_jellyfish_count

    if n_subgenomes != 2:
        logger.warning("  EM: LLR scheme supports 2 groups only (got %d) — "
                       "skipping EM", n_subgenomes)
        return {}

    # Init guard: a nearly-empty group makes the M-step LLR weights
    # statistically meaningless (and costs a full scan + jellyfish run),
    # so skip before doing any work.  The LLR is size-normalized, so an
    # imbalanced but well-populated init (e.g. 52 vs 354) is fine — only
    # the absolute minority size matters.
    if not init_labels:
        logger.warning("  EM: no initialization labels — skipping EM")
        return {}
    group_sizes = [sum(1 for g in init_labels.values() if int(g) == grp)
                   for grp in range(n_subgenomes)]
    if min(group_sizes) < cfg.em_min_init:
        logger.warning("  EM: initialization collapsed (group sizes %s, "
                       "need >= %d each) — skipping EM, keeping init result",
                       group_sizes, cfg.em_min_init)
        return {}

    logger.info("  EM LLR phasing (k=%d, copy band [%d, %d])",
                k, cfg.em_min_copy, cfg.em_max_copy)
    dump_path = _run_jellyfish_count(contig_fasta, k=k,
                                     lower_count=cfg.em_min_copy,
                                     threads=ncpu)
    if not dump_path:
        logger.warning("  jellyfish failed — EM phasing unavailable")
        return {}
    cand, _band_counts = _em_load_candidates(
        dump_path, cfg.em_min_copy, cfg.em_max_copy,
        cfg.em_max_candidates, seed)
    if len(cand) < 100:
        logger.warning("  EM: too few candidate k-mers (%d)", len(cand))
        return {}
    logger.info("  EM: %d candidate k-mers", len(cand))

    # ONE scan; all EM passes iterate on the cached hit arrays in memory.
    t0 = time.time()
    hit_contig, hit_kmer, names, lengths = _em_scan(
        contig_fasta, k, cand, cfg.em_scan_stride)
    n_contigs = len(names)
    n_cand = len(cand)
    logger.info("  EM: scan cached %d candidate hits over %d contigs "
                "in %.0fs", len(hit_contig), n_contigs, time.time() - t0)
    name_to_idx = {nm: i for i, nm in enumerate(names)}

    assign = np.full(n_contigs, -1, dtype=np.int64)
    for nm, g in init_labels.items():
        i = name_to_idx.get(nm)
        if i is not None and 0 <= int(g) < n_subgenomes:
            assign[i] = int(g)

    def m_step():
        """Per-group k-mer copy counts -> LLR weight per candidate."""
        gc = np.zeros((2, n_cand), dtype=np.int64)
        for g in (0, 1):
            sel = assign[hit_contig] == g
            if sel.any():
                gc[g] = np.bincount(hit_kmer[sel], minlength=n_cand)
        # Frequencies normalized by group total hits, smoothed by ~2
        # copies: w = log((c0+2)/size0) - log((c1+2)/size1).
        totals = gc.sum(axis=0)
        sizes = np.maximum(gc.sum(axis=1).astype(np.float64), 1.0)
        w = (np.log((gc[0] + 2.0) / sizes[0])
             - np.log((gc[1] + 2.0) / sizes[1]))
        w[totals < cfg.em_min_total] = 0.0
        return w

    def e_step(w):
        """Summed LLR per contig + informative-hit stats."""
        w_hits = w[hit_kmer]
        s = np.bincount(hit_contig, weights=w_hits, minlength=n_contigs)
        abs_w = np.bincount(hit_contig, weights=np.abs(w_hits),
                            minlength=n_contigs)
        n_inf = np.bincount(hit_contig, weights=w_hits != 0,
                            minlength=n_contigs)
        return s, abs_w, n_inf

    weights = None
    for it in range(max(cfg.em_iterations, 1)):
        weights = m_step()
        s, _abs_w, n_inf = e_step(weights)
        new_assign = np.where(n_inf >= cfg.em_min_markers,
                              (s < 0).astype(np.int64), -1)
        changed = int(((new_assign != assign) & (new_assign >= 0)).sum())
        logger.info("  EM pass %d: LLR k-mers = [%d, %d], %d contigs "
                    "reassigned", it + 1, int((weights > 0).sum()),
                    int((weights < 0).sum()), changed)
        assign = new_assign

    # Final scoring of ALL contigs with the last weights
    s, abs_w, n_inf = e_step(weights)
    pred = np.where(n_inf >= cfg.em_min_markers,
                    (s < 0).astype(np.int64), -1)
    conf = np.abs(s) / np.maximum(abs_w, 1e-12)

    # Post-EM sanity: agreement with the initialization.  EM from a bad
    # init (~50-60% correct) converges to noise and typically reassigns
    # most labeled contigs; a good init keeps high agreement.  Discarding
    # the result here lets the caller keep its pre-EM labels instead of
    # confident garbage (measured: 55% init -> 48% agreement -> noise).
    init_arr = np.full(n_contigs, -1, dtype=np.int64)
    for nm, g in init_labels.items():
        i = name_to_idx.get(nm)
        if i is not None and 0 <= int(g) < n_subgenomes:
            init_arr[i] = int(g)
    both = (init_arr >= 0) & (pred >= 0)
    if both.sum() >= 50:
        agree = float((pred[both] == init_arr[both]).mean())
        if agree < 0.6:
            logger.warning("  EM: only %.0f%% agreement with init — init "
                           "was unreliable, discarding EM result",
                           agree * 100)
            return {}

    # Deterministic sg naming: group with more total contig length -> sg_0
    if contig_lengths is None:
        contig_lengths = {nm: int(lengths[i]) for i, nm in enumerate(names)}
    group_len = Counter()
    for i, g in enumerate(pred):
        if g >= 0:
            group_len[int(g)] += contig_lengths.get(names[i], 1000)
    order = sorted(range(n_subgenomes),
                   key=lambda g: (-group_len.get(g, 0), g))
    sg_of_group = {g: f"sg_{rank}" for rank, g in enumerate(order)}

    result = {}
    for i, nm in enumerate(names):
        g = int(pred[i])
        if g < 0:
            continue
        result[nm] = (sg_of_group[g], float(min(conf[i], 1.0)))
    logger.info("  EM: %d/%d contigs phased (>= %d informative hits)",
                len(result), n_contigs, cfg.em_min_markers)
    return result


# ============================================================================
# Block-guided phasing (Allo4D-style collinear pairs -> orientation -> EM)
# ============================================================================

def run_block_guided_phasing(
    query_fa: str,
    ref_genome: str,
    ref_pep: str,
    work_dir: str,
    *,
    ncpu: int = 4,
    cfg: PhasingConfig | None = None,
) -> dict[str, tuple[str, float]]:
    """Diploid-anchored homeolog-block guided subgenome phasing.

    Chain (validated on a mixed 3rd-gen two-species assembly: 0.94 contig /
    0.998 length-weighted accuracy):

      1. miniprot-annotate the diploid reference and the query with the
         reference proteins (gene-level).
      2. jcvi ortholog anchors -> 1v2 homeolog block pairs
         (Allo4D-style, robust reimplementation in nearscaff.homeolog).
      3. Block-level HG orientation clustering -> labels for the (few
         dozen) long block-host contigs.
      4. Those labels init the EM LLR refinement, whose pooled k-mer
         weights then score EVERY contig — including short ones that
         carry no signal of their own.

    Returns {contig: (sg_label, confidence)}.  When any stage cannot
    produce usable labels (too few pairs, orientation or EM guard
    tripped), returns the best available partial labels, possibly {}.
    """
    import numpy as np

    cfg = cfg or PhasingConfig()
    from nearscaff import homeolog
    from nearscaff.kmer_scan import run_jellyfish_count as _run_jellyfish_count, scan_contig_profiles as _scan_contig_profiles

    os.makedirs(work_dir, exist_ok=True)
    logger.info("Block-guided phasing (annotate -> pairs -> orientation "
                "-> EM)")

    # ---- 1. annotations ----
    ref_pep_gl = homeolog.gene_level_pep(
        ref_pep, os.path.join(work_dir, "ref.gene.pep"))
    _g2, _p2, sp2_cds, sp2_bed = homeolog.annotate_genome(
        ref_genome, ref_pep_gl, work_dir, "sp2", threads=ncpu, prefix="DP")
    _g4, _p4, sp4_cds, sp4_bed = homeolog.annotate_genome(
        query_fa, ref_pep_gl, work_dir, "sp4", threads=ncpu, prefix="MP")

    # ---- 2. collinear homeolog pairs (jcvi) ----
    clusters = homeolog.find_homeolog_pairs(
        sp4_cds, sp2_cds, sp4_bed, sp2_bed, work_dir, threads=ncpu)
    blocks = homeolog.cluster_side_blocks(clusters, sp4_bed)
    if len(blocks) < 10:
        logger.warning("  Too few homeolog blocks (%d) — block-guided "
                       "phasing unavailable", len(blocks))
        return {}

    # ---- 3. contig k-mer profiles ----
    dump_path = _run_jellyfish_count(query_fa, k=cfg.k, lower_count=3,
                                     threads=ncpu)
    if not dump_path:
        logger.warning("  jellyfish failed — block-guided phasing absent")
        return {}
    candidates = _load_candidate_kmers(dump_path, cfg.min_freq, cfg.max_freq)
    kmer_to_idx = {km: i for i, km in enumerate(candidates)}
    contig_names, raw = _scan_contig_profiles(
        query_fa, cfg.k, kmer_to_idx, len(candidates))
    contig_lengths = _read_fasta_lengths(query_fa)
    keep = [i for i, n in enumerate(contig_names)
            if contig_lengths.get(n, 0) >= cfg.min_contig_len]
    contig_names = [contig_names[i] for i in keep]
    Xc = np.array([raw[i] for i in keep], dtype=np.float64)
    lens = np.array([contig_lengths[n] for n in contig_names],
                    dtype=np.float64)
    Xc = Xc / lens[:, np.newaxis]
    cidx = {n: i for i, n in enumerate(contig_names)}

    # ---- 4. block units + orientation ----
    unit_names, rows, pair_list, host = [], [], [], {}
    for cid, sides in sorted(blocks.items()):
        unit_pair = []
        for side in ("A", "B"):
            scaf, _s, _e = sides[side]
            if scaf not in cidx:
                continue
            uname = f"{cid}:{side}"
            unit_names.append(uname)
            rows.append(Xc[cidx[scaf]])
            host[uname] = scaf
            unit_pair.append(uname)
        if len(unit_pair) == 2:
            pair_list.append(tuple(unit_pair))
    if len(pair_list) < 5:
        logger.warning("  Too few block pairs on long contigs (%d)",
                       len(pair_list))
        return {}
    X = np.array(rows)
    pairs = {f"p{i}": set(p) for i, p in enumerate(pair_list)}

    hg_res = hg_orientation_clustering(X, unit_names, pairs,
                                       min_fold=cfg.min_fold,
                                       n_bootstrap=cfg.bootstrap_reps)
    if hg_res is None:
        logger.warning("  Block orientation failed — block-guided "
                       "phasing absent")
        return {}
    labels, _conf = hg_res

    init = {}
    for uname, g in zip(unit_names, labels):
        init[host[uname]] = int(g)
    sizes = Counter(init.values())
    logger.info("  Block orientation: %d block-host contigs init %s",
                len(init), dict(sizes))

    # Honesty gate: the orientation's two groups should cover comparable
    # amounts of sequence.  A heavily skewed split means the orientation
    # collapsed onto one subgenome — report block-host labels at low
    # confidence instead of letting EM amplify garbage genome-wide.
    group_len = Counter()
    for nm, g in init.items():
        group_len[g] += contig_lengths.get(nm, 0)
    if group_len and min(group_len.values()) < 0.25 * max(group_len.values()):
        logger.warning("  Block orientation groups heavily skewed "
                       "(%s bp) — likely unreliable; reporting block-host "
                       "labels at confidence 0.5", dict(group_len))
        order = sorted(range(2), key=lambda g: (-group_len.get(g, 0), g))
        sg_of = {g: f"sg_{r}" for r, g in enumerate(order)}
        return {nm: (sg_of[g], 0.5) for nm, g in init.items()}

    # ---- 5. EM LLR refinement over all contigs ----
    # The em_min_init guard is tuned for contig-scale KMeans inits (many
    # small contigs).  Block-guided inits are few but curated and long —
    # 30 Mb-scale contigs per group give plenty of LLR evidence — so the
    # floor is relaxed here (still >= 10 per side as a sanity bound).
    import dataclasses

    em_cfg = cfg
    if min(sizes.values()) >= 10 and min(sizes.values()) < cfg.em_min_init:
        em_cfg = dataclasses.replace(
            cfg, em_min_init=max(10, min(sizes.values())))
    result = em_marker_phasing(query_fa, k=cfg.em_k, n_subgenomes=2,
                               init_labels=init, contig_lengths=contig_lengths,
                               ncpu=ncpu, cfg=em_cfg)
    if result:
        return result

    # EM skipped/discarded: fall back to block-host labels only
    logger.warning("  EM unavailable — reporting block-host labels only")
    group_len = Counter()
    for nm, g in init.items():
        group_len[g] += contig_lengths.get(nm, 1000)
    order = sorted(range(2), key=lambda g: (-group_len.get(g, 0), g))
    sg_of = {g: f"sg_{r}" for r, g in enumerate(order)}
    return {nm: (sg_of[g], 0.5) for nm, g in init.items()}
