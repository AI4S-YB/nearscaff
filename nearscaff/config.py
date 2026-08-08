"""Configuration dataclasses and divergence presets."""
from dataclasses import dataclass, field


@dataclass
class ProteinConfig:
    min_orf_len: int = 30
    max_intron: str = "auto"
    splice_model: int = 1
    frame_penalty: int = 23
    max_secondary: int = 30


@dataclass
class SyntenyConfig:
    cscore_threshold: float = 0.7
    gene_dist_x: str = "auto"
    gene_dist_y: str = "auto"
    min_cluster_size: int = 4
    max_iterations: int = 100


@dataclass
class BlockTreeConfig:
    interleave_overlap: float = 0.5
    max_inter_block_gap: int = 500000
    min_synteny_score: float = 50.0
    orientation_penalty: float = 0.2
    interleaving_penalty: float = 0.3
    nucleotide_weight: float = 0.0     # weight for nuc identity diff (disabled: too noisy)
    protein_weight: float = 0.0        # weight for protein identity diff (disabled: too noisy)
    min_divergence: float = 0.001      # minimum divergence for k-mer Jaccard-based edges


@dataclass
class NucleotideConfig:
    aligner: str = "minimap2"
    preset: str = "asm5"
    region_margin: int = 50000
    boundary_window: int = 2000
    nucleotide_passes: list[str] = field(default_factory=lambda: ["asm5", "asm20"])
    secondary_alignments: int = 5          # minimap2 -N：挽回 scaffold 附近的次级比对
    reuse_ref_index: bool = True           # 构建/复用 minimap2 参考索引
    index_dir: str | None = None           # 索引与缓存目录；None 落到 <output>/intermediate/


@dataclass
class ScaffoldConfig:
    gap_min: int = 0
    gap_max: int = 500000
    unknown_gap_size: int = 100
    candidate_scaling: float = 0.75
    best_buddy_scale: bool = True


@dataclass
class PhasingConfig:
    """Configuration for k-mer subgenome phasing (kmer-phase)."""
    enabled: bool = False           # standalone module only
    k: int = 15
    min_freq: int = 200             # min genome-wide total copies for candidate k-mers
    max_freq: int = 1_000_000_000   # max genome-wide total copies
    min_fold: float = 2.0           # max/submax fold within homology groups
    group_ratio: float = 0.5        # fraction of groups where fold must hold
    min_contig_len: int = 10000     # contigs shorter than this are not clustered
    bootstrap_reps: int = 100
    bootstrap_threshold: float = 0.7  # contigs below this confidence stay unphased
    specific_pval: float = 0.05     # t-test p-value for SG-specific k-mers
    semi_supervised: bool = False   # opt-in legacy RF path (off by default)
    # EM LLR refinement (for closely related subgenomes whose
    # repetitive-k-mer spectra are not yet diverged).  Signal lives in
    # copy-band k-mers: per-group copy frequencies define a
    # log-likelihood-ratio weight per k-mer (M-step), and each contig is
    # scored by the summed LLR of its k-mer hits (E-step).  A single M->E
    # pass is best — further iterations degrade good initializations.
    # em_iterations=0 disables.
    em_k: int = 21                  # longer k: more fixed-difference markers
    em_min_copy: int = 3            # candidate k-mer copy band (genome-wide)
    em_max_copy: int = 200000
    em_marker_frac: float = 0.75    # (legacy binary-marker threshold, unused
                                    # by the LLR scheme)
    em_min_markers: int = 20        # min informative (nonzero-weight) k-mer
                                    # hits for a contig to be phased
    em_min_total: int = 3           # k-mer total copy floor for a nonzero
                                    # LLR weight
    em_iterations: int = 1          # number of M->E passes
    em_scan_stride: int = 2         # sample every Nth k-mer position in scans
    em_max_candidates: int = 30_000_000  # subsample cap for candidates
    em_min_init: int = 50           # min contigs per group in the EM init


@dataclass
class NearscaffConfig:
    protein: ProteinConfig = field(default_factory=ProteinConfig)
    synteny: SyntenyConfig = field(default_factory=SyntenyConfig)
    blocktree: BlockTreeConfig = field(default_factory=BlockTreeConfig)
    nucleotide: NucleotideConfig = field(default_factory=NucleotideConfig)
    scaffold: ScaffoldConfig = field(default_factory=ScaffoldConfig)
    phasing: PhasingConfig = field(default_factory=PhasingConfig)
    threads: int = 4
    output_dir: str = "nearscaff_output"
    keep_intermediate: bool = False


DIVERGENCE_PRESETS = {
    "close": {"miniprot_j": 1, "minimap2_preset": "asm5"},
    "near":  {"miniprot_j": 1, "minimap2_preset": "asm10"},
    "mid":   {"miniprot_j": 2, "minimap2_preset": "asm15"},
    "far":   {"miniprot_j": 2, "minimap2_preset": "asm20"},
}
