"""Unit tests for kmer_phasing.py — SubPhaser-style Signal A phasing."""

import os
import random
import shutil
import tempfile

import pytest

from nearscaff.config import PhasingConfig
from nearscaff.kmer_phasing import (
    _load_candidate_kmers,
    best_label_permutation,
    fold_filter_kmers,
    normalize_cluster_labels,
    zscore_columns,
)

np = pytest.importorskip("numpy", reason="k-mer matrix tests need numpy")


# ===========================================================================
# Fold-change filter
# ===========================================================================

def _make_fold_matrix():
    """8 contigs (4 homology groups of 2) x 4 k-mers.

    kmer0: differential in all 4 groups (zero in second member)
    kmer1: uniform across all contigs
    kmer2: differential in only 1 group
    kmer3: differential in exactly 2 groups
    """
    names = [f"c{i}" for i in range(8)]
    groups = {
        "g1": {"c0", "c1"},
        "g2": {"c2", "c3"},
        "g3": {"c4", "c5"},
        "g4": {"c6", "c7"},
    }
    X = np.zeros((8, 4))
    X[:, 1] = 1.0  # uniform
    # kmer0: high in first member of each group, zero in second
    for i in (0, 2, 4, 6):
        X[i, 0] = 2.0
    # kmer2: differential only in g1
    X[0, 2] = 2.0
    # kmer3: differential in g1 and g2
    X[0, 3] = 2.0
    X[2, 3] = 2.0
    return X, names, groups


class TestFoldFilter:
    def test_selects_differential_kmers(self):
        X, names, groups = _make_fold_matrix()
        selected = fold_filter_kmers(X, names, groups, min_fold=2.0,
                                     group_ratio=0.5)
        assert set(selected) == {0, 3}

    def test_zero_frequency_counts_as_differential(self):
        """A k-mer present in one homeolog and absent from the other must
        pass the fold filter (epsilon denominator), not be dropped."""
        X, names, groups = _make_fold_matrix()
        selected = fold_filter_kmers(X, names, groups, min_fold=100.0,
                                     group_ratio=0.5)
        assert 0 in set(selected)

    def test_raw_abundance_filter(self):
        X, names, groups = _make_fold_matrix()
        raw = np.array([100.0, 1000.0, 100.0, 100.0])
        selected = fold_filter_kmers(X, names, groups, min_fold=2.0,
                                     group_ratio=0.5,
                                     raw_col_sums=raw, min_total_freq=200)
        assert set(selected) == set()

    def test_n_features_cap(self):
        X, names, groups = _make_fold_matrix()
        selected = fold_filter_kmers(X, names, groups, min_fold=2.0,
                                     group_ratio=0.5, n_features=1)
        assert len(selected) == 1

    def test_too_few_groups_keeps_all(self):
        X, names, groups = _make_fold_matrix()
        selected = fold_filter_kmers(X, names, {"g1": groups["g1"]},
                                     min_fold=2.0, group_ratio=0.5)
        assert len(selected) == X.shape[1]


# ===========================================================================
# Z-score + KMeans + bootstrap (synthetic two clusters)
# ===========================================================================

class TestClusterBootstrap:
    def test_two_clusters_recovered(self):
        pytest.importorskip("sklearn")
        from nearscaff.kmer_phasing import cluster_with_bootstrap

        rng = np.random.RandomState(7)
        n_per, n_kmers = 12, 40
        X = np.zeros((2 * n_per, n_kmers))
        # cluster A high on first half of k-mers, cluster B on second half
        X[:n_per, :n_kmers // 2] = 2.0
        X[n_per:, n_kmers // 2:] = 2.0
        X += rng.normal(0, 0.05, X.shape)

        X_z = zscore_columns(X)
        labels, confidence = cluster_with_bootstrap(
            X_z, 2, n_replicates=20, random_state=42)

        # cluster IDs are arbitrary: check agreement with truth up to flip
        truth = np.array([0] * n_per + [1] * n_per)
        agree = max(np.mean(labels == truth), np.mean(labels != truth))
        assert agree >= 0.95
        assert confidence.min() >= 0.8

    def test_zscore_equalizes_columns(self):
        X = np.array([[1.0, 100.0], [3.0, 200.0], [5.0, 300.0]])
        Z = zscore_columns(X)
        assert np.allclose(Z.mean(axis=0), 0)
        assert np.allclose(Z.std(axis=0), 1)


# ===========================================================================
# Permutation alignment + cluster naming
# ===========================================================================

class TestBestLabelPermutation:
    def test_flipped_two_labels(self):
        truth = {f"c{i}": ("sg_0" if i < 5 else "sg_1") for i in range(10)}
        flipped = {c: ("sg_1" if t == "sg_0" else "sg_0")
                   for c, t in truth.items()}
        mapping = best_label_permutation(flipped, truth)
        assert mapping == {"sg_1": "sg_0", "sg_0": "sg_1"}

    def test_three_labels_rotation(self):
        truth = {f"c{i}": f"sg_{i % 3}" for i in range(30)}
        rot = {"sg_0": "sg_2", "sg_1": "sg_0", "sg_2": "sg_1"}
        source = {c: rot[t] for c, t in truth.items()}
        mapping = best_label_permutation(source, truth)
        assert all(mapping[rot[f"sg_{i}"]] == f"sg_{i}" for i in range(3))

    def test_no_shared_keys(self):
        assert best_label_permutation({"a": "sg_0"}, {"b": "sg_0"}) == {}

    def test_unequal_label_sets_greedy(self):
        truth = {f"c{i}": ("sg_0" if i < 8 else "sg_1") for i in range(10)}
        source = dict(truth)
        source["c0"] = "sg_5"  # extra label, one contig
        mapping = best_label_permutation(source, truth)
        assert mapping.get("sg_1") == "sg_1"


class TestNormalizeClusterLabels:
    def test_seed_labels_flip_recovered(self):
        names = [f"c{i}" for i in range(10)]
        lengths = {n: 1000 for n in names}
        labels = np.array([0] * 5 + [1] * 5)
        # seeds say the opposite of raw cluster IDs
        seeds = {f"c{i}": ("sg_1" if i < 5 else "sg_0") for i in range(10)}
        sg = normalize_cluster_labels(labels, names, lengths, seeds, 2)
        assert sg == ["sg_1"] * 5 + ["sg_0"] * 5

    def test_no_seeds_length_descending(self):
        names = ["a1", "a2", "b1", "b2"]
        lengths = {"a1": 100, "a2": 100, "b1": 5000, "b2": 5000}
        labels = np.array([0, 0, 1, 1])
        sg = normalize_cluster_labels(labels, names, lengths, None, 2)
        assert sg == ["sg_1", "sg_1", "sg_0", "sg_0"]

    def test_too_few_seeds_falls_back_to_length(self):
        names = [f"c{i}" for i in range(10)]
        lengths = {n: 1000 for n in names}
        labels = np.array([0] * 5 + [1] * 5)
        seeds = {"c0": "sg_1", "c1": "sg_1"}  # < 10 seeds
        sg = normalize_cluster_labels(labels, names, lengths, seeds, 2)
        assert sg == ["sg_0"] * 5 + ["sg_1"] * 5


# ===========================================================================
# HG-aware orientation clustering
# ===========================================================================

def _make_hg_matrix(n_directional=300, n_noise=200, seed=13):
    """8 contigs (4 pairs) x k-mers.

    Directional k-mers: enriched in the same-side member of every pair
    (member 0 of each pair = group 0 truth).  Noise k-mers: differential
    in one pair only, random side.  Returns (X, names, groups, truth).
    """
    rng = np.random.RandomState(seed)
    names = [f"c{i}" for i in range(8)]
    groups = {f"g{p}": {f"c{2 * p}", f"c{2 * p + 1}"} for p in range(4)}
    pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]
    n_kmers = n_directional + n_noise
    X = rng.uniform(0.9, 1.1, (8, n_kmers))  # ~uniform background
    for j in range(n_directional):
        for a, b in pairs:
            X[a, j] = 3.0
            X[b, j] = 0.3
    for j in range(n_directional, n_kmers):
        a, b = pairs[rng.randint(4)]
        if rng.rand() < 0.5:
            a, b = b, a
        X[a, j] = 3.0
        X[b, j] = 0.3
    truth = np.array([0, 1, 0, 1, 0, 1, 0, 1])
    return X, names, groups, truth


class TestHgOrientationClustering:
    def test_recovers_groups_and_splits_pairs(self):
        from nearscaff.kmer_phasing import hg_orientation_clustering

        X, names, groups, truth = _make_hg_matrix()
        res = hg_orientation_clustering(X, names, groups, min_fold=2.0,
                                        n_bootstrap=20)
        assert res is not None
        labels, confidence = res
        agree = max(np.mean(labels == truth), np.mean(labels != truth))
        assert agree == 1.0
        assert confidence.min() >= 0.8

    def test_pair_local_noise_does_not_dominate(self):
        """Noise k-mers outnumbering directional ones must not flip the
        global orientation (the KMeans failure mode)."""
        from nearscaff.kmer_phasing import hg_orientation_clustering

        X, names, groups, truth = _make_hg_matrix(n_directional=300,
                                                  n_noise=2000)
        res = hg_orientation_clustering(X, names, groups, min_fold=2.0,
                                        n_bootstrap=10)
        assert res is not None
        labels, _conf = res
        agree = max(np.mean(labels == truth), np.mean(labels != truth))
        assert agree == 1.0

    def test_too_few_pairs_returns_none(self):
        from nearscaff.kmer_phasing import hg_orientation_clustering

        X, names, groups, _truth = _make_hg_matrix()
        res = hg_orientation_clustering(X, names, {"g1": groups["g1"]},
                                        min_fold=2.0, n_bootstrap=0)
        assert res is None

    def test_too_few_directional_kmers_returns_none(self):
        from nearscaff.kmer_phasing import hg_orientation_clustering

        X, names, groups, _truth = _make_hg_matrix(n_directional=50,
                                                   n_noise=0)
        res = hg_orientation_clustering(X, names, groups, min_fold=2.0,
                                        n_bootstrap=0)
        assert res is None


class TestEnforceHgConstraints:
    def test_same_side_pair_split_at_weaker_member(self):
        from nearscaff.kmer_phasing import enforce_hg_constraints

        names = [f"c{i}" for i in range(4)]
        groups = {"g1": {"c0", "c1"}, "g2": {"c2", "c3"}}
        labels = np.array([0, 0, 0, 1])
        conf = np.array([0.9, 0.8, 0.95, 0.9])
        new_labels, new_conf = enforce_hg_constraints(labels, conf, names,
                                                      groups)
        assert new_labels[0] != new_labels[1]
        assert new_labels[0] == 0            # stronger member keeps label
        assert new_labels[1] == 1
        assert new_conf[1] == 0.5            # flipped member flagged
        assert new_labels[2] == 0 and new_labels[3] == 1  # untouched

    def test_already_split_pairs_untouched(self):
        from nearscaff.kmer_phasing import enforce_hg_constraints

        names = ["c0", "c1"]
        groups = {"g1": {"c0", "c1"}}
        labels = np.array([1, 0])
        conf = np.array([0.9, 0.85])
        new_labels, new_conf = enforce_hg_constraints(labels, conf, names,
                                                      groups)
        assert list(new_labels) == [1, 0]
        assert list(new_conf) == [0.9, 0.85]


# ===========================================================================
# SG-specific k-mers (t-test)
# ===========================================================================

class TestSpecificKmers:
    def test_top1_vs_top2(self):
        pytest.importorskip("scipy")
        from nearscaff.kmer_phasing import find_specific_kmers

        rng = np.random.RandomState(3)
        # 10 contigs per group; kmer0 specific to sg_0, kmer1 not specific
        X = np.zeros((20, 2))
        X[:10, 0] = 2.0 + rng.normal(0, 0.01, 10)
        X[10:, 0] = 0.1 + rng.normal(0, 0.01, 10)
        X[:, 1] = 1.0 + rng.normal(0, 0.01, 20)
        labels = ["sg_0"] * 10 + ["sg_1"] * 10

        specific = find_specific_kmers(X, ["kmer0", "kmer1"], labels,
                                       pval=0.05)
        by_kmer = {k: sg for k, sg, _p, _m in specific}
        assert by_kmer.get("kmer0") == "sg_0"
        assert "kmer1" not in by_kmer


# ===========================================================================
# Candidate k-mer loading (repetitive filter, no top-N truncation)
# ===========================================================================

class TestLoadCandidateKmers:
    def test_freq_window(self):
        with tempfile.NamedTemporaryFile(suffix=".fa", mode="w",
                                         delete=False) as f:
            f.write("AAAA 50\nCCCC 200\nGGGG 5000\nTTTT 99999999999\n")
            path = f.name
        try:
            kmers = _load_candidate_kmers(path, min_freq=200, max_freq=10**9)
            assert kmers == ["CCCC", "GGGG"]
        finally:
            os.unlink(path)


# ===========================================================================
# Integration: full phasing on a synthetic two-subgenome FASTA
# (needs jellyfish + sklearn)
# ===========================================================================

def _has_tool(name):
    return shutil.which(name) is not None


HAS_JELLYFISH = _has_tool("jellyfish")


def _synthetic_subgenome_fasta(path, n_contigs=6, contig_len=2000, seed=11):
    """Two subgenomes with group-specific 40-mer repeats."""
    rng = random.Random(seed)
    bases = "ACGT"
    repeat_a = "".join(rng.choice(bases) for _ in range(40))
    repeat_b = "".join(rng.choice(bases) for _ in range(40))
    truth = {}
    with open(path, "w") as f:
        for sg, repeat in (("sgA", repeat_a), ("sgB", repeat_b)):
            for i in range(n_contigs):
                name = f"{sg}_ctg{i}"
                bg = [rng.choice(bases) for _ in range(contig_len // 2)]
                seq = "".join(bg) + repeat * (contig_len // 2 // 40)
                f.write(f">{name}\n")
                for j in range(0, len(seq), 60):
                    f.write(seq[j:j + 60] + "\n")
                truth[name] = sg
    return truth


@pytest.mark.skipif(not HAS_JELLYFISH, reason="jellyfish not installed")
class TestRunSubgenomeKmerPhasing:
    def test_end_to_end_two_subgenomes(self):
        pytest.importorskip("sklearn")
        from nearscaff.kmer_phasing import run_subgenome_kmer_phasing

        with tempfile.TemporaryDirectory() as d:
            fasta = os.path.join(d, "query.fa")
            truth = _synthetic_subgenome_fasta(fasta)
            cfg = PhasingConfig(k=15, min_freq=50, min_contig_len=500,
                                bootstrap_reps=10, bootstrap_threshold=0.5)
            result = run_subgenome_kmer_phasing(
                fasta, k=15, n_subgenomes=2, ncpu=2, cfg=cfg, output_dir=d)

            assert len(result) >= len(truth) * 0.8
            pred = {c: sg for c, (sg, _conf) in result.items()}
            mapping = best_label_permutation(pred, truth)
            correct = sum(1 for c in pred if mapping.get(pred[c]) == truth[c])
            assert correct / len(pred) >= 0.9
            assert os.path.exists(os.path.join(d, "kmer_specific.tsv"))


@pytest.mark.skipif(not HAS_JELLYFISH, reason="jellyfish not installed")
class TestEmMarkerPhasing:
    @staticmethod
    def _marker_fasta(path, k=21, n_markers=100, n_shared=60,
                      n_contigs=8, positions=120, seed=5):
        """Two groups; group markers at medium copy, shared k-mers in both."""
        rng = random.Random(seed)
        bases = "ACGT"

        def rand_kmer():
            return "".join(rng.choice(bases) for _ in range(k))

        markers = {"A": set(), "B": set()}
        while len(markers["A"]) < n_markers:
            markers["A"].add(rand_kmer())
        while True:
            km = rand_kmer()
            if km not in markers["A"]:
                markers["B"].add(km)
            if len(markers["B"]) == n_markers:
                break
        shared = set()
        while len(shared) < n_shared:
            km = rand_kmer()
            if km not in markers["A"] and km not in markers["B"]:
                shared.add(km)
        mk_a, mk_b, sh = sorted(markers["A"]), sorted(markers["B"]), sorted(shared)

        truth = {}
        with open(path, "w") as f:
            for group, own in (("A", mk_a), ("B", mk_b)):
                for i in range(n_contigs):
                    name = f"sg{group}_ctg{i}"
                    parts = []
                    for _ in range(positions):
                        if rng.random() < 0.85:
                            parts.append(rng.choice(own))
                        else:
                            parts.append(rng.choice(sh))
                        parts.append("".join(rng.choice(bases)
                                             for _ in range(rng.randint(0, 9))))
                    f.write(f">{name}\n{''.join(parts)}\n")
                    truth[name] = f"sg{group}"
        return truth

    def test_em_recovers_groups_noisy_init(self):
        """Single-pass LLR EM refines a noisy (70% correct) initialization,
        as produced by the degenerate KMeans stage on closely related
        subgenomes."""
        from nearscaff.kmer_phasing import em_marker_phasing, best_label_permutation

        with tempfile.TemporaryDirectory() as d:
            fasta = os.path.join(d, "query.fa")
            truth = self._marker_fasta(fasta, n_contigs=60)
            rng = random.Random(3)
            init = {}
            for name, sg in truth.items():
                g = 0 if sg == "sgA" else 1
                if rng.random() < 0.3:
                    g = 1 - g
                init[name] = g
            cfg = PhasingConfig(em_k=21, em_min_copy=3,
                                em_max_copy=200000,
                                em_min_markers=10, em_iterations=1,
                                em_scan_stride=1)
            result = em_marker_phasing(fasta, k=21, n_subgenomes=2,
                                       init_labels=init, ncpu=2, cfg=cfg)

            assert len(result) >= len(truth) * 0.8
            pred = {c: sg for c, (sg, _c) in result.items()}
            mapping = best_label_permutation(pred, truth)
            correct = sum(1 for c in pred if mapping.get(pred[c]) == truth[c])
            assert correct / len(pred) >= 0.9
            # confidence must be a proper fraction in [0, 1]
            assert all(0.0 <= c <= 1.0 for _sg, c in result.values())

    def test_em_skips_collapsed_init(self):
        """A collapsed init (one group nearly empty) must skip EM entirely
        and return {} so the caller keeps its own (KMeans) result."""
        from nearscaff.kmer_phasing import em_marker_phasing

        with tempfile.TemporaryDirectory() as d:
            fasta = os.path.join(d, "query.fa")
            truth = self._marker_fasta(fasta, n_contigs=60)
            init = {}
            for i, (name, sg) in enumerate(sorted(truth.items())):
                init[name] = 1 if i == 0 else 0  # 1 vs 119: collapsed
            cfg = PhasingConfig(em_k=21, em_min_markers=10,
                                em_iterations=1, em_scan_stride=1)
            result = em_marker_phasing(fasta, k=21, n_subgenomes=2,
                                       init_labels=init, ncpu=2, cfg=cfg)
            assert result == {}

    def test_em_skips_without_init(self):
        """No init labels at all -> EM skipped ({}), no random restarts."""
        from nearscaff.kmer_phasing import em_marker_phasing

        with tempfile.TemporaryDirectory() as d:
            fasta = os.path.join(d, "query.fa")
            self._marker_fasta(fasta, n_contigs=60)
            cfg = PhasingConfig(em_k=21, em_min_markers=10,
                                em_iterations=1, em_scan_stride=1)
            result = em_marker_phasing(fasta, k=21, n_subgenomes=2,
                                       init_labels=None, ncpu=2, cfg=cfg)
            assert result == {}
