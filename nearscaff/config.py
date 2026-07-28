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
    nucleotide_passes: list[str] = field(default_factory=lambda: ["asm5", "asm10", "asm20"])


@dataclass
class ScaffoldConfig:
    gap_min: int = 0
    gap_max: int = 500000
    unknown_gap_size: int = 100
    candidate_scaling: float = 0.75
    best_buddy_scale: bool = True


@dataclass
class NearscaffConfig:
    protein: ProteinConfig = field(default_factory=ProteinConfig)
    synteny: SyntenyConfig = field(default_factory=SyntenyConfig)
    blocktree: BlockTreeConfig = field(default_factory=BlockTreeConfig)
    nucleotide: NucleotideConfig = field(default_factory=NucleotideConfig)
    scaffold: ScaffoldConfig = field(default_factory=ScaffoldConfig)
    threads: int = 4
    output_dir: str = "nearscaff_output"
    keep_intermediate: bool = False


DIVERGENCE_PRESETS = {
    "close": {"miniprot_j": 1, "minimap2_preset": "asm5"},
    "near":  {"miniprot_j": 1, "minimap2_preset": "asm10"},
    "mid":   {"miniprot_j": 2, "minimap2_preset": "asm15"},
    "far":   {"miniprot_j": 2, "minimap2_preset": "asm20"},
}
