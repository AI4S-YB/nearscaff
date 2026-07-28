"""Block Tree Builder — the core innovation of nearly_scaffold.

Partitions synteny blocks into a 3-layer hierarchy:
  Layer 1: Subgenome partition (INTERLEAVED detection + k-coloring)
  Layer 2: Chromosome assignment (majority rule)
  Layer 3: LocalCluster construction (adjacent same-orientation merge)

The INTERLEAVED signal is THE key discriminator for subgenome splitting
in polyploid genomes: two blocks that overlap on the same reference
chromosome but carry disjoint query contigs indicate homeologous regions
from different subgenomes.
"""

import json
import logging
import random
from collections import defaultdict

from nearscaff.types import SyntenyBlock, BlockTreeNode

logger = logging.getLogger("nearscaff.blocktree")


# ---------------------------------------------------------------------------
# Step 1: INTERLEAVED Detection
# ---------------------------------------------------------------------------

def detect_interleaved(
    b1: SyntenyBlock,
    b2: SyntenyBlock,
    overlap_threshold: float = 0.5,
) -> str | None:
    """Return 'INTERLEAVED' if *b1* and *b2* are interleaved, else None.

    Two synteny blocks are INTERLEAVED when:

    * They map to the same reference chromosome.
    * Their reference intervals overlap by more than *overlap_threshold*
      (overlap length / length of the shorter block).
    * Their query contig sets are disjoint (no shared contigs).
    * Their orientations are consistent (both '+' or both '-').

    Returns
    -------
    str or None
        ``"INTERLEAVED"`` when the criteria are met, otherwise ``None``.
    """
    # 1. Same reference chromosome
    if b1.ref_chr != b2.ref_chr:
        return None

    # 2. Handle zero- or negative-length blocks gracefully
    len1 = b1.ref_end - b1.ref_start
    len2 = b2.ref_end - b2.ref_start
    if len1 <= 0 or len2 <= 0:
        return None

    # 3. Reference interval overlap
    overlap_start = max(b1.ref_start, b2.ref_start)
    overlap_end = min(b1.ref_end, b2.ref_end)
    overlap_len = max(0, overlap_end - overlap_start)
    shorter_len = min(len1, len2)
    if overlap_len / shorter_len < overlap_threshold:
        return None

    # 4. Disjoint query contigs
    if b1.query_contigs & b2.query_contigs:
        return None

    # 5. Consistent orientation
    if b1.orientation != b2.orientation:
        return None

    return "INTERLEAVED"


# ---------------------------------------------------------------------------
# Step 2: Conflict Graph
# ---------------------------------------------------------------------------

def build_conflict_graph(
    blocks: list[SyntenyBlock],
    overlap_threshold: float = 0.5,
) -> dict:
    """Build a conflict graph from pairwise INTERLEAVED detection.

    Parameters
    ----------
    blocks : list[SyntenyBlock]
        All synteny blocks discovered during Stage 1.
    overlap_threshold : float
        Minimum overlap ratio for INTERLEAVED classification.

    Returns
    -------
    dict
        Nested dict ``{chr_name: {block_id: {conflicting_block_ids}}}``.
    """
    graph: dict = defaultdict(lambda: defaultdict(set))

    n = len(blocks)
    for i in range(n):
        bi = blocks[i]
        for j in range(i + 1, n):
            bj = blocks[j]
            if detect_interleaved(bi, bj, overlap_threshold) == "INTERLEAVED":
                chr_name = bi.ref_chr
                graph[chr_name][bi.block_id].add(bj.block_id)
                graph[chr_name][bj.block_id].add(bi.block_id)

    # Convert inner defaultdict(set) to plain dict of sets
    return {chr_name: dict(adj) for chr_name, adj in graph.items()}


# ---------------------------------------------------------------------------
# Step 2b: Weighted Divergence Graph (replaces binary INTERLEAVED detection)
# ---------------------------------------------------------------------------

def build_homology_regions_from_blocks(
    blocks: list[SyntenyBlock],
) -> dict[str, list[list[int]]]:
    """Group blocks by reference chromosome into homology regions.

    A homology region is a set of blocks whose reference intervals overlap.
    Overlap is transitive: if A overlaps B and B overlaps C, then A, B, C
    form one region.

    Returns ``{chr_name: [[block_index, ...], ...]}`` — each inner list
    holds indices into the *blocks* list.  Single-block regions are dropped.
    """
    if not blocks:
        return {}

    # Group block indices by chromosome
    chr_indices: dict[str, list[int]] = defaultdict(list)
    for i, b in enumerate(blocks):
        chr_indices[b.ref_chr].append(i)

    regions: dict[str, list[list[int]]] = {}
    for chr_name, idx_list in chr_indices.items():
        # Sort by ref_start
        idx_list.sort(key=lambda i: blocks[i].ref_start)

        chr_regions: list[list[int]] = []
        current: list[int] = [idx_list[0]]
        current_end = blocks[idx_list[0]].ref_end

        for idx in idx_list[1:]:
            b = blocks[idx]
            if b.ref_start < current_end:
                # Overlaps current region — join
                current.append(idx)
                current_end = max(current_end, b.ref_end)
            else:
                # No overlap — flush current, start new
                if len(current) >= 2:
                    chr_regions.append(current)
                current = [idx]
                current_end = b.ref_end
        if len(current) >= 2:
            chr_regions.append(current)

        if chr_regions:
            regions[chr_name] = chr_regions

    total_regions = sum(len(v) for v in regions.values())
    logger.info("  %d homology regions across %d chromosomes",
                total_regions, len(regions))
    return regions


def build_divergence_graph(
    blocks: list[SyntenyBlock],
    homology_regions: dict[str, list[list[int]]],
    nucleotide_weight: float = 0.5,
    protein_weight: float = 0.5,
    min_divergence: float = 0.05,
    query_fasta: str = "",
    kmer_jaccard_k: int = 21,
    kmer_jaccard_weight: float = 1.0,
    kmer_jaccard_samples: int = 10000,
) -> dict[str, dict[str, dict[str, float]]]:
    """Build weighted divergence graph from homology regions.

    For each pair of blocks in the same homology region, the pairwise
    divergence is computed from three signals:
      - nucleotide identity difference (from minimap2 PAF)
      - protein identity difference (from miniprot anchors)
      - k-mer Jaccard distance (contig-to-contig, within the region)

    Edges with weight >= *min_divergence* are kept.

    Returns ``{chr_name: {block_id: {neighbor_id: divergence}}}``.
    """
    graph: dict[str, dict[str, dict[str, float]]] = defaultdict(
        lambda: defaultdict(dict))

    # Pre-compute block k-mer sets if query_fasta is provided
    block_kmer_sets: dict[str, set[str]] = {}
    if query_fasta:
        block_kmer_sets = _build_block_kmer_sets(
            blocks, query_fasta, k=kmer_jaccard_k,
            n_samples=kmer_jaccard_samples, stride=10)

    for chr_name, chr_regions in homology_regions.items():
        for region in chr_regions:
            n = len(region)
            for pi in range(n):
                for pj in range(pi + 1, n):
                    bi = blocks[region[pi]]
                    bj = blocks[region[pj]]

                    # Same orientation required
                    if bi.orientation != bj.orientation:
                        continue

                    # Compute divergence from identity differences
                    nuc_div = abs(bi.nucleotide_identity -
                                  bj.nucleotide_identity)
                    prot_div = abs(bi.protein_identity -
                                   bj.protein_identity)

                    # K-mer Jaccard distance (contig-level sequence comparison)
                    kmer_div = 0.0
                    has_kmer = False
                    if block_kmer_sets:
                        kset_i = block_kmer_sets.get(bi.block_id, set())
                        kset_j = block_kmer_sets.get(bj.block_id, set())
                        if kset_i and kset_j:
                            has_kmer = True
                            intersection = len(kset_i & kset_j)
                            union = len(kset_i | kset_j)
                            if union > 0:
                                kmer_div = 1.0 - intersection / union

                    # Fallback: if both blocks have no identity data,
                    # use binary INTERLEAVED check with default threshold
                    if bi.nucleotide_identity == 0.0 and \
                       bj.nucleotide_identity == 0.0 and \
                       bi.protein_identity == 0.0 and \
                       not has_kmer:
                        result = detect_interleaved(bi, bj, 0.5)
                        if result == "INTERLEAVED":
                            divergence = 1.0
                        else:
                            continue
                    else:
                        # Normalize weights so they sum to 1
                        w_sum = nucleotide_weight + protein_weight
                        if has_kmer:
                            w_sum += kmer_jaccard_weight
                            divergence = (
                                nucleotide_weight * nuc_div +
                                protein_weight * prot_div +
                                kmer_jaccard_weight * kmer_div
                            ) / w_sum
                        else:
                            divergence = (
                                nucleotide_weight * nuc_div +
                                protein_weight * prot_div
                            ) / max(w_sum, 0.001)

                    if divergence >= min_divergence:
                        graph[chr_name][bi.block_id][bj.block_id] = divergence
                        graph[chr_name][bj.block_id][bi.block_id] = divergence

    # Convert defaultdict to plain dict
    return {chr: dict(adj) for chr, adj in graph.items()}


def _build_block_kmer_sets(
    blocks: list[SyntenyBlock],
    query_fasta: str,
    k: int = 21,
    n_samples: int = 10000,
    stride: int = 100,
) -> dict[str, set[str]]:
    """Build a sampled k-mer set for each block from its contig sequences.

    Uses reservoir-style sampling: every k-mer has an equal chance of being
    kept, with a hard cap of *n_samples* per block to bound memory.
    """
    random.seed(42)

    contig_to_block: dict[str, str] = {}
    for b in blocks:
        for ctg in b.query_contigs:
            contig_to_block[ctg] = b.block_id

    # Reserve space: list of k-mers per block for reservoir sampling
    block_samples: dict[str, list[str]] = {b.block_id: [] for b in blocks}
    block_seen: dict[str, int] = {b.block_id: 0 for b in blocks}

    current_name = None
    current_seq = []
    with open(query_fasta) as f:
        for line in f:
            if line.startswith('>'):
                if current_name and current_seq:
                    bid = contig_to_block.get(current_name)
                    if bid:
                        seq = ''.join(current_seq).upper()
                        _sample_kmers_reservoir(
                            seq, k, block_samples[bid], block_seen,
                            bid, n_samples, stride=stride)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_name and current_seq:
            bid = contig_to_block.get(current_name)
            if bid:
                seq = ''.join(current_seq).upper()
                _sample_kmers_reservoir(
                    seq, k, block_samples[bid], block_seen,
                    bid, n_samples, stride=stride)

    result = {bid: set(samples) for bid, samples in block_samples.items()}
    sizes = {bid: len(s) for bid, s in result.items()}
    logger.info("  Block k-mer sets: %d blocks, sizes range %d–%d",
                len(result), min(sizes.values()) if sizes else 0,
                max(sizes.values()) if sizes else 0)
    return result


def _sample_kmers_reservoir(
    seq: str, k: int, samples: list[str], seen: dict,
    block_id: str, cap: int, stride: int = 100,
):
    """Strided reservoir-sample k-mers from *seq* into *samples*, capped at *cap*.

    Only every *stride*-th k-mer is considered, reducing processing 100×.
    """
    if len(seq) < k:
        return
    for i in range(0, len(seq) - k + 1, stride):
        kmer = seq[i:i + k]
        rev = _revcomp(kmer)
        canonical = kmer if kmer < rev else rev
        n = seen[block_id]
        seen[block_id] = n + 1
        if len(samples) < cap:
            samples.append(canonical)
        else:
            j = random.randint(0, n)
            if j < cap:
                samples[j] = canonical


def _revcomp(seq: str) -> str:
    """Reverse complement."""
    cmap = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(cmap)[::-1]


def _greedy_k_coloring_weighted(
    blocks: list[SyntenyBlock],
    adj: dict[str, dict[str, float]],
    k: int,
) -> dict[str, int]:
    """Greedy k-coloring with weighted conflict minimization.

    Like :func:`_greedy_k_coloring` but sums divergence weights instead
    of counting binary conflicts.  Length-balance tie-breaking is preserved.
    """
    color_of: dict[str, int] = {}
    color_total_len: dict[int, int] = {c: 0 for c in range(k)}

    for b in blocks:
        block_len = b.ref_end - b.ref_start
        conflict_score: dict[int, float] = {c: 0.0 for c in range(k)}
        for nbr, weight in adj.get(b.block_id, {}).items():
            if nbr in color_of:
                conflict_score[color_of[nbr]] += weight

        min_score = min(conflict_score.values())
        # Among colors with minimum conflict, pick the smallest so far
        best_color = min(
            (c for c in range(k) if conflict_score[c] == min_score),
            key=lambda c: color_total_len[c],
        )
        color_of[b.block_id] = best_color
        color_total_len[best_color] += block_len

    return color_of


# ---------------------------------------------------------------------------
# Internal helpers for colouring
# ---------------------------------------------------------------------------

def _find_connected_components(adj: dict) -> list[list]:
    """Find connected components in an undirected graph.

    Parameters
    ----------
    adj : dict
        Adjacency dict ``{node_id: {neighbor_ids}}``.

    Returns
    -------
    list[list]
        Each sub-list is one connected component (list of node ids).
    """
    visited: set = set()
    components: list[list] = []

    for node in adj:
        if node in visited:
            continue
        comp: list = []
        stack = [node]
        while stack:
            v = stack.pop()
            if v not in visited:
                visited.add(v)
                comp.append(v)
                stack.extend(adj.get(v, set()) - visited)
        components.append(comp)

    return components


def _greedy_k_coloring(
    blocks: list[SyntenyBlock],
    adj: dict,
    k: int,
) -> dict:
    """Greedy k-coloring — assign each block to the colour with fewest conflicts.

    Tie-breaking is by smallest *total length* currently assigned to the
    colour, so blocks are balanced across subgenomes rather than all piling
    into colour 0.
    """
    color_of: dict = {}
    color_total_len: dict[int, int] = {c: 0 for c in range(k)}

    for b in blocks:
        block_len = b.ref_end - b.ref_start
        conflicts_per_color = {c: 0 for c in range(k)}
        for nbr in adj.get(b.block_id, set()):
            if nbr in color_of:
                conflicts_per_color[color_of[nbr]] += 1

        min_conflicts = min(conflicts_per_color.values())
        # Pick among least-conflict colours: prefer smallest current size
        best_color = min(
            (c for c in range(k) if conflicts_per_color[c] == min_conflicts),
            key=lambda c: color_total_len[c],
        )
        color_of[b.block_id] = best_color
        color_total_len[best_color] += block_len

    return color_of


def _majority_orientation(blocks: list[SyntenyBlock]) -> str:
    """Return the majority orientation ('+' or '-') across *blocks*."""
    if not blocks:
        return "+"
    plus = sum(1 for b in blocks if b.orientation == "+")
    minus = len(blocks) - plus
    return "+" if plus >= minus else "-"


# ---------------------------------------------------------------------------
# Step 3: Subgenome Assignment
# ---------------------------------------------------------------------------

def assign_subgenomes(
    blocks: list[SyntenyBlock],
    conflict_graph: dict,
    ploidy_hint: int | None = None,
    weighted: bool = False,
) -> list[BlockTreeNode]:
    """Assign every synteny block to a subgenome via k-coloring.

    Two modes:

    **Binary** (*weighted*=False, default): connected components of the
    conflict graph determine *k*, and :func:`_greedy_k_coloring` does
    binary conflict counting.

    **Weighted** (*weighted*=True): the conflict graph carries float
    divergence weights.  *k* is taken from *ploidy_hint* (or 2 by
    default), and :func:`_greedy_k_coloring_weighted` sums weights.

    Parameters
    ----------
    blocks : list[SyntenyBlock]
        All blocks.
    conflict_graph : dict
        ``{chr_name: {block_id: ...}}`` — structure depends on *weighted*.
    ploidy_hint : int or None
        If given, *k* is at least this value.
    weighted : bool
        If True, *conflict_graph* is a weighted adjacency.

    Returns
    -------
    list[BlockTreeNode]
        One ``BlockTreeNode`` per assigned subgenome (level ``"subgenome"``).
    """
    if not blocks:
        return []

    # Group blocks by reference chromosome
    chr_blocks: dict = defaultdict(list)
    for b in blocks:
        chr_blocks[b.ref_chr].append(b)

    # Determine global k
    if weighted:
        # Count actual conflict edges; if graph is empty, don't force k=2
        total_edges = sum(
            len(neighbors)
            for ch in conflict_graph.values()
            for neighbors in ch.values()
        )
        if total_edges == 0 and ploidy_hint is None:
            k = 1
        else:
            k = max(2, ploidy_hint or 2)
    else:
        max_components = 0
        for chr_name, chr_adj in conflict_graph.items():
            comps = _find_connected_components(chr_adj)
            max_components = max(max_components, len(comps))
        k = max_components
        if k < 1:
            k = 1
        if ploidy_hint is not None:
            k = max(k, ploidy_hint)

    # Colour every chromosome's blocks
    color_of: dict = {}
    for chr_name, chr_list in chr_blocks.items():
        chr_adj = conflict_graph.get(chr_name, {})
        if weighted:
            chr_colors = _greedy_k_coloring_weighted(chr_list, chr_adj, k)
        else:
            chr_colors = _greedy_k_coloring(chr_list, chr_adj, k)
        color_of.update(chr_colors)

    # Build subgenome nodes — only for colours that actually have blocks
    sg_nodes: list[BlockTreeNode] = []
    for i in range(k):
        sg_blocks = [b for b in blocks if color_of.get(b.block_id, 0) == i]
        if not sg_blocks:
            continue

        sg_node = BlockTreeNode(
            node_id=f"sg_{i}",
            level="subgenome",
            query_contigs=sorted(set().union(*(b.query_contigs for b in sg_blocks))),
            gene_pairs=[gp for b in sg_blocks for gp in b.gene_pairs],
            synteny_score=(
                sum(b.score for b in sg_blocks) / len(sg_blocks)
                if sg_blocks
                else 0.0
            ),
            orientation=_majority_orientation(sg_blocks),
        )
        # Stash block ids for downstream layer construction
        sg_node._block_ids = [b.block_id for b in sg_blocks]  # type: ignore[attr-defined]
        sg_nodes.append(sg_node)

    return sg_nodes


# ---------------------------------------------------------------------------
# Step 4: Chromosome Assignment
# ---------------------------------------------------------------------------

def build_chromosome_layer(
    sg_node: BlockTreeNode,
    blocks: list[SyntenyBlock],
) -> list[BlockTreeNode]:
    """Populate *sg_node.children* with chromosome-level nodes.

    Within a subgenome, blocks are grouped by ``ref_chr``.  Each group
    becomes one chromosome node (level ``"chromosome"``) and is attached
    as a child of *sg_node*.

    Parameters
    ----------
    sg_node : BlockTreeNode
        A subgenome-level node with ``_block_ids`` stashed.
    blocks : list[SyntenyBlock]
        All blocks (will be filtered to those belonging to *sg_node*).

    Returns
    -------
    list[BlockTreeNode]
        The newly-created chromosome nodes.
    """
    sg_block_ids = set(getattr(sg_node, "_block_ids", []))
    sg_blocks = [b for b in blocks if b.block_id in sg_block_ids]

    # Group by reference chromosome
    chr_groups: dict = defaultdict(list)
    for b in sg_blocks:
        chr_groups[b.ref_chr].append(b)

    chr_nodes: list[BlockTreeNode] = []
    for chr_name, chr_list in chr_groups.items():
        chr_node = BlockTreeNode(
            node_id=f"{sg_node.node_id}_{chr_name}",
            level="chromosome",
            ref_chr=chr_name,
            ref_start=min(b.ref_start for b in chr_list),
            ref_end=max(b.ref_end for b in chr_list),
            query_contigs=sorted(set().union(*(b.query_contigs for b in chr_list))),
            gene_pairs=[gp for b in chr_list for gp in b.gene_pairs],
            synteny_score=(
                sum(b.score for b in chr_list) / len(chr_list)
                if chr_list
                else 0.0
            ),
            orientation=_majority_orientation(chr_list),
        )
        chr_node._block_ids = [b.block_id for b in chr_list]  # type: ignore[attr-defined]
        sg_node.add_child(chr_node)
        chr_nodes.append(chr_node)

    return chr_nodes


# ---------------------------------------------------------------------------
# Step 5: LocalCluster Construction
# ---------------------------------------------------------------------------

def build_local_clusters(
    chr_node: BlockTreeNode,
    blocks: list[SyntenyBlock],
    max_inter_block_gap: int = 500_000,
) -> list[BlockTreeNode]:
    """Populate *chr_node.children* with local-cluster nodes.

    Blocks belonging to *chr_node* are sorted by ``ref_start``.
    Adjacent blocks are merged into a single ``LocalCluster`` when:

    * The gap (*next.ref_start* - *prev.ref_end*) is < *max_inter_block_gap*.
    * They share the same orientation.

    Parameters
    ----------
    chr_node : BlockTreeNode
        A chromosome-level node with ``_block_ids`` stashed.
    blocks : list[SyntenyBlock]
        All blocks (will be filtered to those belonging to *chr_node*).
    max_inter_block_gap : int
        Maximum gap (bp) tolerated between adjacent blocks in a cluster.

    Returns
    -------
    list[BlockTreeNode]
        The newly-created local-cluster nodes.
    """
    chr_block_ids = set(getattr(chr_node, "_block_ids", []))
    chr_blocks = sorted(
        [b for b in blocks if b.block_id in chr_block_ids],
        key=lambda b: b.ref_start,
    )

    if not chr_blocks:
        return []

    # Merge adjacent blocks
    clusters: list[list[SyntenyBlock]] = []
    current_cluster = [chr_blocks[0]]

    for b in chr_blocks[1:]:
        prev = current_cluster[-1]
        gap = b.ref_start - prev.ref_end
        if gap < max_inter_block_gap and b.orientation == prev.orientation:
            current_cluster.append(b)
        else:
            clusters.append(current_cluster)
            current_cluster = [b]
    clusters.append(current_cluster)

    local_nodes: list[BlockTreeNode] = []
    for i, cluster in enumerate(clusters):
        local_node = BlockTreeNode(
            node_id=f"{chr_node.node_id}_lc{i}",
            level="local",
            ref_chr=chr_node.ref_chr,
            ref_start=min(b.ref_start for b in cluster),
            ref_end=max(b.ref_end for b in cluster),
            query_contigs=sorted(
                set().union(*(b.query_contigs for b in cluster))
            ),
            gene_pairs=[gp for b in cluster for gp in b.gene_pairs],
            synteny_score=(
                sum(b.score for b in cluster) / len(cluster)
                if cluster
                else 0.0
            ),
            orientation=_majority_orientation(cluster),
        )
        chr_node.add_child(local_node)
        local_nodes.append(local_node)

    return local_nodes


# ---------------------------------------------------------------------------
# Step 6: Confidence Scoring
# ---------------------------------------------------------------------------

def _compute_confidence_scores(
    root: BlockTreeNode,
    blocks: list[SyntenyBlock],
) -> None:
    """Assign confidence scores to every node in the tree.

    Leaf (local) nodes are scored from their constituent blocks using
    synteny continuity, anchor density, and orientation discord.
    Parent nodes average their children's scores.
    """
    w1, w2, w3 = 0.4, 0.3, 0.2
    block_map = {b.block_id: b for b in blocks}

    # --- leaf: local clusters ---
    for node in root.iter_level("local"):
        bids = getattr(node, "_block_ids", [])
        node_blocks = [block_map[bid] for bid in bids if bid in block_map]
        if not node_blocks:
            node.confidence = 0.5
            continue

        span = (node.ref_end or 0) - (node.ref_start or 0)
        if span <= 0:
            node.confidence = 0.5
            continue

        # continuity: fraction of span covered by actual blocks
        total_block_len = sum(b.ref_end - b.ref_start for b in node_blocks)
        continuity = min(total_block_len / span, 1.0)

        # anchor density (normalised — 0.01 anchors/bp ≈ 1.0)
        total_anchors = sum(b.anchor_count for b in node_blocks)
        density = total_anchors / span
        norm_density = min(density / 0.01, 1.0)

        # orientation discord
        plus_cnt = sum(1 for b in node_blocks if b.orientation == "+")
        discord = 1.0 - max(plus_cnt, len(node_blocks) - plus_cnt) / len(node_blocks)

        conf = w1 * continuity + w2 * norm_density - w3 * discord
        node.confidence = max(0.0, min(1.0, conf))

    # --- chromosome: average of local children ---
    for node in root.iter_level("chromosome"):
        if node.children:
            node.confidence = sum(c.confidence for c in node.children) / len(node.children)
        else:
            node.confidence = 0.5

    # --- subgenome: average of chromosome children ---
    for node in root.iter_level("subgenome"):
        if node.children:
            node.confidence = sum(c.confidence for c in node.children) / len(node.children)
        else:
            node.confidence = 0.5

    # --- root: average of subgenome children ---
    if root.children:
        root.confidence = sum(c.confidence for c in root.children) / len(root.children)
    else:
        root.confidence = 0.5


# ---------------------------------------------------------------------------
# Step 7: Main Entry Point
# ---------------------------------------------------------------------------

def build_block_tree(
    blocks: list[SyntenyBlock],
    ref_chromosomes: list[str],
    ploidy_hint: int | None = None,
    overlap_threshold: float = 0.5,
    max_inter_block_gap: int = 500_000,
    nucleotide_weight: float = 0.5,
    protein_weight: float = 0.5,
    min_divergence: float = 0.05,
    query_fasta: str = "",
) -> BlockTreeNode:
    """Build the full 3-layer Block Tree from a flat list of synteny blocks.

    By default uses a **weighted divergence graph** that replaces binary
    INTERLEAVED detection with continuous nucleotide + protein identity
    difference scoring.  Falls back to binary INTERLEAVED only when no
    blocks carry identity data.

    Parameters
    ----------
    blocks : list[SyntenyBlock]
        Synteny blocks discovered by Stage 1.
    ref_chromosomes : list[str]
        Ordered list of reference chromosome names.
    ploidy_hint : int or None
        Minimum subgenome count (e.g. 4 for tetraploid).
    overlap_threshold : float
        Fallback: overlap ratio for binary INTERLEAVED detection.
    max_inter_block_gap : int
        Maximum gap (bp) in LocalCluster merging (default 500 kbp).
    nucleotide_weight : float
        Weight for nucleotide identity diff in divergence (default 0.5).
    protein_weight : float
        Weight for protein identity diff in divergence (default 0.5).
    min_divergence : float
        Minimum combined divergence to create a conflict edge (default 0.05).

    Returns
    -------
    BlockTreeNode
        Root node of the Block Tree (level ``"root"``).
    """
    root = BlockTreeNode(node_id="root", level="root")

    if not blocks:
        return root

    # Build the hierarchy
    if ploidy_hint == 1:
        # Haploid: skip conflict detection, single subgenome
        conflict_graph = {}
        weighted = False
    else:
        homology_regions = build_homology_regions_from_blocks(blocks)
        conflict_graph = build_divergence_graph(
            blocks, homology_regions,
            nucleotide_weight=nucleotide_weight,
            protein_weight=protein_weight,
            min_divergence=min_divergence,
            query_fasta=query_fasta,
        )
        weighted = True
    sg_nodes = assign_subgenomes(blocks, conflict_graph, ploidy_hint,
                                  weighted=weighted)

    for sg_node in sg_nodes:
        root.add_child(sg_node)
        sg_label = sg_node.node_id
        sg_node.sg_label = sg_label
        chr_nodes = build_chromosome_layer(sg_node, blocks)
        for chr_node in chr_nodes:
            chr_node.sg_label = sg_label
            build_local_clusters(chr_node, blocks, max_inter_block_gap)
            for lc_node in chr_node.children:
                lc_node.sg_label = sg_label

    # Compute confidence scores
    _compute_confidence_scores(root, blocks)

    return root


# ---------------------------------------------------------------------------
# JSON Serialization (Stage 1 → Stage 2 handoff)
# ---------------------------------------------------------------------------

def block_tree_to_json(root: BlockTreeNode) -> str:
    """Serialize a Block Tree to a JSON string.

    Parameters
    ----------
    root : BlockTreeNode
        Root node of the Block Tree.

    Returns
    -------
    str
        JSON-encoded string representation.
    """
    return json.dumps(root.to_dict(), indent=2, default=str)


def _node_from_dict(data: dict) -> BlockTreeNode:
    """Recursively reconstruct a :class:`BlockTreeNode` from a dictionary."""
    node = BlockTreeNode(
        node_id=data["node_id"],
        level=data["level"],
        ref_chr=data.get("ref_chr"),
        ref_start=data.get("ref_start"),
        ref_end=data.get("ref_end"),
        query_contigs=data.get("query_contigs", []),
        orientation=data.get("orientation", "+"),
        synteny_score=data.get("synteny_score", 0.0),
        confidence=data.get("confidence", 0.0),
        flags=set(data.get("flags", [])),
        sg_label=data.get("sg_label"),
        gene_pairs=[],
    )
    # Restore gene pair count as dummy entries (exact pairs not needed
    # for tree structure operations)
    gp_count = data.get("gene_pair_count", 0)
    if gp_count > 0:
        node.gene_pairs = [("", "", 0)] * gp_count

    for child_data in data.get("children", []):
        node.add_child(_node_from_dict(child_data))

    return node


def block_tree_from_json(json_str: str) -> BlockTreeNode:
    """Deserialize a JSON string back into a Block Tree.

    Parameters
    ----------
    json_str : str
        JSON string produced by :func:`block_tree_to_json`.

    Returns
    -------
    BlockTreeNode
        Reconstructed root node.
    """
    data = json.loads(json_str)
    return _node_from_dict(data)
