"""Core data types for the nearscaff pipeline."""
from dataclasses import dataclass, field
from typing import Optional
from enum import Enum


class EdgeType(str, Enum):
    PROTEIN_SYNTENY = "protein_synteny"
    NUCLEOTIDE_CHAIN = "nucleotide_chain"
    BLOCK_TREE_ADJACENCY = "block_tree_adjacency"
    HIC_LINKAGE = "hic_linkage"


EDGE_TYPE_WEIGHTS = {
    EdgeType.PROTEIN_SYNTENY: 1.0,
    EdgeType.NUCLEOTIDE_CHAIN: 0.85,
    EdgeType.BLOCK_TREE_ADJACENCY: 0.5,
    EdgeType.HIC_LINKAGE: 0.3,
}


@dataclass
class GeneAnchor:
    """Single gene-level anchor from protein alignment."""
    query_contig: str
    query_gene: str
    ref_chr: str
    ref_gene: str
    q_start: int
    q_end: int
    r_start: int
    r_end: int
    strand: str
    score: float
    identity: float
    n_exons: int
    n_frameshifts: int
    n_stop_codons: int

    @property
    def cscore(self) -> float:
        return getattr(self, '_cscore', 0.0)


@dataclass
class SyntenyBlock:
    """A collinear block of gene anchors."""
    block_id: str
    ref_chr: str
    ref_start: int
    ref_end: int
    query_contigs: set[str]
    gene_pairs: list  # [(query_gene, ref_gene, score), ...]
    orientation: str
    score: float
    anchor_count: int
    nucleotide_identity: float = 0.0  # mean nmatch/hitlen from minimap2 PAF
    protein_identity: float = 0.0     # mean miniprot identity from GeneAnchor


@dataclass
class BlockTreeNode:
    """Hierarchical node in the Block Tree."""
    node_id: str
    level: str  # 'subgenome' | 'chromosome' | 'local'
    ref_chr: Optional[str] = None
    ref_start: Optional[int] = None
    ref_end: Optional[int] = None
    query_contigs: list[str] = field(default_factory=list)
    orientation: str = '+'
    gene_pairs: list = field(default_factory=list)
    synteny_score: float = 0.0
    parent: Optional['BlockTreeNode'] = None
    children: list['BlockTreeNode'] = field(default_factory=list)
    confidence: float = 0.0
    flags: set[str] = field(default_factory=set)
    sg_label: Optional[str] = None  # subgenome label: "SG01", "SG02", etc.

    def add_child(self, child: 'BlockTreeNode'):
        child.parent = self
        self.children.append(child)

    def iter_level(self, level: str):
        """Iterate all nodes at a given level in this subtree."""
        if self.level == level:
            yield self
        for child in self.children:
            yield from child.iter_level(level)

    def to_dict(self) -> dict:
        """Serialize to dict for JSON output."""
        return {
            'node_id': self.node_id,
            'level': self.level,
            'ref_chr': self.ref_chr,
            'ref_start': self.ref_start,
            'ref_end': self.ref_end,
            'query_contigs': self.query_contigs,
            'orientation': self.orientation,
            'synteny_score': self.synteny_score,
            'confidence': self.confidence,
            'flags': list(self.flags),
            'sg_label': self.sg_label,
            'gene_pair_count': len(self.gene_pairs),
            'children': [c.to_dict() for c in self.children],
        }
