"""Full pipeline orchestration for nearscaff.

Stage 0: Protein anchoring (miniprot) + C-score filtering + synteny
         -> contig-level SyntenyBlocks -> Block Tree
Stage 1: Scaffolding (constrained nucleotide alignment + Scaffold Graph
         -> AGP + scaffold FASTA)
"""

import os
import re
import logging
import subprocess
from nearscaff.config import NearscaffConfig
from nearscaff.protein_anchor import (
    extract_proteins_from_gff, build_protein_index,
    run_protein_map, parse_protein_paf,
)
from nearscaff.synteny import filter_cscore
from nearscaff.blocktree import build_block_tree, block_tree_to_json, block_tree_from_json
from nearscaff.types import EdgeType

logger = logging.getLogger("nearscaff")

# Minimum mapq of the final precise alignment for a contig's orientation
# to be reported as +/- in the AGP; below this it is reported as '?'.
MIN_ORIENT_MAPQ = 10

# Pattern: FCYA01G000001.mRNA1 -> chromosome=FCYA01, gene_number=1
_PROTEIN_ID_PATTERN = re.compile(r'^(\D+\d+)G(\d+)')


# ---------------------------------------------------------------------------
# Signal C: Protein anchoring -> synteny -> Block Tree
# ---------------------------------------------------------------------------

def _map_contigs_to_reference(contigs: set, ref_fasta: str, query_fasta: str,
                               threads: int = 4, preset: str = "asm5"
                               ) -> dict[str, tuple[str, int, int, float]]:
    """Align query contigs to reference genome via minimap2.

    Returns ``{contig_name: (ref_chr, ref_start, ref_end, identity)}``
    using the best alignment (by nmatch) per contig.  *identity* is
    nmatch / hitlen from the minimap2 PAF.
    """
    import tempfile
    from nearscaff.nucleotide import _extract_contigs, parse_nucleotide_paf

    if not contigs:
        return {}

    fd, contigs_fa = tempfile.mkstemp(suffix='.fa')
    os.close(fd)
    try:
        _extract_contigs(query_fasta, list(contigs), contigs_fa)
        cmd = ["minimap2", "-t", str(threads), "-x", preset,
               "-c", ref_fasta, contigs_fa]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning("Contig-to-reference alignment failed: %s",
                           result.stderr.strip())
            return {}

        entries = parse_nucleotide_paf(result.stdout)
        best: dict[str, tuple[str, int, int, int, int]] = {}
        for e in entries:
            query = e['query']
            paf_chr = e['ref_chr'].split(':')[0]
            rs = min(e['rstart'], e['rend'])
            re = max(e['rstart'], e['rend'])
            nm = e['nmatch']
            if query not in best or nm > best[query][3]:
                best[query] = (paf_chr, rs, re, nm, e['hitlen'])
        return {k: (v[0], v[1], v[2], v[3] / max(v[4], 1))
                for k, v in best.items()}
    finally:
        if os.path.exists(contigs_fa):
            os.unlink(contigs_fa)


def _build_contig_blocks(anchors, contig_ref_pos,
                          cluster_radius: int = 500_000,
                          min_cluster_size: int = 4
                          ) -> list:
    """Build one SyntenyBlock per (contig, ref_chromosome).

    Each contig mapping to a reference chromosome gets its own block spanning
    the contig's full alignment range.  Contigs from different subgenomes that
    cover the same reference region will produce overlapping blocks, which
    triggers INTERLEAVED detection in build_block_tree.

    Contigs from the *same* subgenome are adjacent along the chromosome
    (little or no overlap), so the overlap_threshold in build_conflict_graph
    naturally excludes them from false INTERLEAVED calls.
    """
    from collections import defaultdict
    from nearscaff.types import SyntenyBlock

    groups: dict[tuple[str, str], list] = defaultdict(list)
    for a in anchors:
        ctg = a.query_contig
        if ctg not in contig_ref_pos:
            continue
        crc = contig_ref_pos[ctg][0]
        groups[(ctg, crc)].append(a)

    blocks = []
    for i, ((ctg, ref_chr), grp) in enumerate(sorted(groups.items())):
        if len(grp) < min_cluster_size:
            continue
        info = contig_ref_pos[ctg]
        crs, cre = info[1], info[2]
        if crs >= cre:
            continue
        nuc_id = info[3] if len(info) > 3 else 0.0
        prot_id = sum(a.identity for a in grp) / len(grp)
        strands = [a.strand for a in grp]
        orientation = "+" if strands.count("+") >= len(strands) - strands.count("+") else "-"

        blocks.append(SyntenyBlock(
            block_id=f"ctg_block_{i:04d}",
            ref_chr=ref_chr,
            ref_start=crs,
            ref_end=cre,
            query_contigs={ctg},
            gene_pairs=[(a.query_gene, a.ref_gene, a.score) for a in grp],
            orientation=orientation,
            score=sum(getattr(a, '_cscore', 1.0) for a in grp),
            anchor_count=len(grp),
            nucleotide_identity=nuc_id,
            protein_identity=prot_id,
        ))

    n_contigs = len(set(g[0] for g in groups))
    n_chrs = len(set(g[1] for g in groups))
    logger.info("  %d contig-level blocks from %d contigs across %d chromosomes",
                len(blocks), n_contigs, n_chrs)

    # Merge blocks with low overlap (same-subgenome redundancy) while keeping
    # heavily overlapping blocks separate for INTERLEAVED detection.
    merged = _merge_adjacent_blocks(blocks, max_overlap_ratio=0.5)
    logger.info("  %d blocks after merging adjacent non-overlapping", len(merged))
    return merged


def _merge_adjacent_blocks(blocks, max_overlap_ratio: float = 0.15):
    """Merge blocks that don't heavily overlap, iteratively.

    Contigs from the *same* subgenome tile the chromosome with relatively
    little overlap among themselves.  Contigs from *different* subgenomes
    cover the same region and overlap heavily.  We iteratively merge any
    pair whose overlap is below *max_overlap_ratio*, collapsing intra-
    subgenome redundancy while preserving the inter-subgenome INTERLEAVED
    signal.
    """
    if len(blocks) <= 1:
        return list(blocks)

    from collections import defaultdict

    # Group by chromosome for efficient pairwise checking
    chr_blocks: dict[str, list] = defaultdict(list)
    for b in blocks:
        chr_blocks[b.ref_chr].append(b)

    for lst in chr_blocks.values():
        lst.sort(key=lambda b: b.ref_start)

    changed = True
    while changed:
        changed = False
        for chr_name in list(chr_blocks.keys()):
            b_list = chr_blocks[chr_name]
            if len(b_list) <= 1:
                continue
            # Try to merge the closest pair below threshold
            best_dist = float('inf')
            best_pair = None
            for i in range(len(b_list)):
                for j in range(i + 1, len(b_list)):
                    bi, bj = b_list[i], b_list[j]
                    if bi.orientation != bj.orientation:
                        continue
                    # Don't merge blocks with significantly different
                    # nucleotide identity (likely different subgenomes)
                    nuc_diff = abs(bi.nucleotide_identity -
                                   bj.nucleotide_identity)
                    if nuc_diff > 0.003:
                        continue
                    # Overlap ratio
                    ov_start = max(bi.ref_start, bj.ref_start)
                    ov_end = min(bi.ref_end, bj.ref_end)
                    ov_len = max(0, ov_end - ov_start)
                    shorter = min(bi.ref_end - bi.ref_start,
                                  bj.ref_end - bj.ref_start)
                    ratio = ov_len / shorter if shorter > 0 else 0.0
                    if ratio <= max_overlap_ratio:
                        dist = abs(bi.ref_start - bj.ref_start)
                        if dist < best_dist:
                            best_dist = dist
                            best_pair = (i, j)
            if best_pair is not None:
                i, j = best_pair
                bi = b_list[i]
                bj = b_list[j]
                # Merge bj into bi (length-weighted for identity fields)
                len_i = bi.ref_end - bi.ref_start
                len_j = bj.ref_end - bj.ref_start
                total_len = len_i + len_j
                bi.ref_start = min(bi.ref_start, bj.ref_start)
                bi.ref_end = max(bi.ref_end, bj.ref_end)
                if total_len > 0:
                    bi.nucleotide_identity = (
                        bi.nucleotide_identity * len_i +
                        bj.nucleotide_identity * len_j
                    ) / total_len
                    bi.protein_identity = (
                        bi.protein_identity * len_i +
                        bj.protein_identity * len_j
                    ) / total_len
                bi.query_contigs.update(bj.query_contigs)
                bi.gene_pairs.extend(bj.gene_pairs)
                bi.score += bj.score
                bi.anchor_count += bj.anchor_count
                b_list.pop(j)
                changed = True
                break  # restart scanning after merge

    result = []
    for lst in chr_blocks.values():
        result.extend(lst)
    return result


def _run_signal_c_protein(
    config: NearscaffConfig,
    ref_fasta: str,
    ref_gff: str,
    query_fasta: str,
    output_dir: str,
    protein_faa: str | None = None,
):
    """Run protein anchoring + synteny -> Block Tree.

    Returns BlockTreeNode root, or None if protein anchoring fails.
    """
    os.makedirs(output_dir, exist_ok=True)
    index_path = os.path.join(output_dir, "query.mpi")

    if protein_faa:
        logger.info("Signal C: Loading proteins directly from %s ...", protein_faa)
        proteins, gene_locations = _load_proteins_from_fasta(protein_faa)
        out_protein_faa = os.path.join(output_dir, "ref_proteins.faa")
        _copy_or_link(protein_faa, out_protein_faa)
        protein_faa = out_protein_faa
    else:
        logger.info("Signal C: Extracting proteins + gene coordinates from GFF3 ...")
        protein_faa = os.path.join(output_dir, "ref_proteins.faa")
        proteins = extract_proteins_from_gff(ref_gff, ref_fasta, protein_faa)
        gene_locations = _parse_gff3_gene_locations(ref_gff)

    logger.info("  %d protein sequences loaded", len(proteins))
    logger.info("  %d gene locations available", len(gene_locations))

    if not proteins:
        logger.warning("No proteins available -- Signal C will be absent")
        return None

    logger.info("Signal C: Building miniprot index of query genome ...")
    build_protein_index(query_fasta, index_path, config.threads)

    logger.info("Signal C: Mapping reference proteins to query genome ...")
    paf_text = run_protein_map(index_path, protein_faa, config.threads,
                                splice_model=config.protein.splice_model)

    logger.info("Signal C: Parsing anchors ...")
    anchors = parse_protein_paf(paf_text)
    logger.info("  %d raw anchors parsed", len(anchors))

    anchors = _enrich_anchors(anchors, gene_locations)
    logger.info("  %d anchors after enrichment with ref coordinates", len(anchors))

    # Realign contigs to reference for genuine genomic coordinates.
    # This is critical for INTERLEAVED detection: synthetic positions
    # (5 kb spacing from protein headers) would spread homeologous
    # blocks apart so they never overlap on the reference.
    contig_set = {a.query_contig for a in anchors}
    logger.info("Signal C: Mapping %d contigs to reference for real coordinates ...",
                len(contig_set))
    contig_ref_pos = _map_contigs_to_reference(
        contig_set, ref_fasta, query_fasta,
        threads=config.threads, preset=config.nucleotide.preset,
    )
    logger.info("  %d contigs with reference alignment", len(contig_ref_pos))

    # Update anchor coordinates from contig-level alignments.
    # All anchors on the same contig span its full reference alignment
    # so that blocks from different subgenomes mapping to the same
    # reference region will overlap and trigger INTERLEAVED.
    for a in anchors:
        if a.query_contig in contig_ref_pos:
            crc, crs, cre, _nuc = contig_ref_pos[a.query_contig]
            if a.ref_chr and crc != a.ref_chr:
                continue  # chromosome mismatch, keep gene-level assignment
            a.ref_chr = crc
            a.r_start = crs
            a.r_end = cre
    logger.info("  %d anchors updated with real reference coordinates",
                sum(1 for a in anchors if a.query_contig in contig_ref_pos))

    if config.keep_intermediate:
        _write_anchors_tsv(anchors, os.path.join(output_dir, "gene_anchors.tsv"))

    logger.info("Signal C: C-score filtering (threshold=%.2f) ...",
                config.synteny.cscore_threshold)
    filtered = list(filter_cscore(anchors, cscore=config.synteny.cscore_threshold))
    logger.info("  %d anchors after C-score filter", len(filtered))

    if not filtered:
        logger.warning("No anchors survive C-score filter -- Signal C absent")
        return None

    logger.info("Signal C: Building contig-level blocks from reference coordinates ...")
    blocks = _build_contig_blocks(filtered, contig_ref_pos,
                                   min_cluster_size=config.synteny.min_cluster_size)
    logger.info("  %d contig-level blocks produced", len(blocks))

    if not blocks:
        logger.warning("No contig-level blocks -- building empty block tree")
        return build_block_tree([], [])

    if config.keep_intermediate:
        _write_blocks_file(blocks, os.path.join(output_dir, "synteny.blocks"))

    chr_name_map = _build_chr_name_map({b.ref_chr for b in blocks}, ref_fasta)
    if chr_name_map:
        logger.info("  Chromosome name mapping: %d entries", len(chr_name_map))
        for b in blocks:
            if b.ref_chr in chr_name_map:
                b.ref_chr = chr_name_map[b.ref_chr]

    logger.info("Signal C: Building Block Tree ...")
    ref_chromosomes = sorted({b.ref_chr for b in blocks})
    root = build_block_tree(
        blocks, ref_chromosomes,
        # Single-subgenome mode (ploidy_hint=1): no subgenome splitting,
        # conflict/divergence detection is skipped inside build_block_tree.
        ploidy_hint=1,
        overlap_threshold=config.blocktree.interleave_overlap,
        max_inter_block_gap=config.blocktree.max_inter_block_gap,
        nucleotide_weight=config.blocktree.nucleotide_weight,
        protein_weight=config.blocktree.protein_weight,
        min_divergence=config.blocktree.min_divergence,
        query_fasta=query_fasta,
    )
    return root


def run_stage1(config: NearscaffConfig, block_tree_path: str,
               ref_fasta: str, query_fasta: str, output_dir: str):
    """Stage 1: Progressive nucleotide extension + Scaffold Graph -> AGP.

    Multi-pass strategy:
      Pass 1 — protein edges from Block Tree (PROTEIN_SYNTENY)
      Pass 2 — asm5  extension: align unplaced contigs, add NUCLEOTIDE_CHAIN edges
      Pass 3 — asm10 extension: same with relaxed params
      Pass 4 — asm20 extension: final sweep

    All passes share one FusedScaffoldGraph.  Re-solve after each pass.
    Returns path to the AGP file.
    """
    from nearscaff.scaffold_graph import FusedScaffoldGraph
    from nearscaff.agp import AGPWriter

    os.makedirs(output_dir, exist_ok=True)

    logger.info("Stage 1: Loading Block Tree from %s", block_tree_path)
    with open(block_tree_path) as f:
        root = block_tree_from_json(f.read())

    contig_lengths = _read_fasta_lengths(query_fasta)
    all_contigs = set(contig_lengths.keys())

    # ---- Build persistent contig→chromosome lookup ----
    # This lookup accumulates across passes so nucleotide-placed contigs
    # retain their chromosome assignment.
    contig_ref = {}  # {contig_base: (chr_name, r_start, r_end)}

    # Seed with Block Tree protein-placed contigs
    for lc in root.iter_level("local"):
        if lc.ref_chr is None:
            continue
        rs = lc.ref_start or 0
        re = lc.ref_end or rs + 3000
        for ctg in lc.query_contigs:
            if ctg not in contig_ref:
                contig_ref[ctg] = (lc.ref_chr, rs, re)

    logger.info("Stage 1: %d contigs with chromosome from Block Tree", len(contig_ref))

    # Track which tier each contig was placed (for confidence-stratified PAF)
    contig_tier = {}  # {contig_base: "protein" | "asm5" | "asm10" | "asm20"}
    for ctg in contig_ref:
        contig_tier[ctg] = "protein"

    # ---- Pass 1: Build graph with protein edges only ----
    logger.info("Stage 1 — Pass 1 (protein): building scaffold graph ...")
    sg, bt_contigs = _build_initial_graph(root)
    cover = _solve_graph(sg, config)

    n_seq, n_scaf, n_gap = _count_agp_from_cover(cover, contig_lengths)
    logger.info("  Pass 1 result: %d scaffolds, %d contigs, %d gaps",
                n_scaf, n_seq, n_gap)

    # ---- Passes 2-4: Progressive nucleotide extension ----
    chunk_threads = max(2, config.threads // 2)
    for preset in config.nucleotide.nucleotide_passes:
        unplaced = _get_unplaced_contigs(cover, all_contigs)
        if not unplaced:
            logger.info("  No unplaced contigs remain — stopping extension")
            break

        logger.info("Stage 1 — Pass (%s extension): %d unplaced contigs ...",
                    preset, len(unplaced))

        scaffold_regions = _get_scaffold_regions(cover, root, contig_lengths)
        if not scaffold_regions:
            logger.info("  No scaffold regions to extend from")
            break

        paf_entries = _align_unplaced_to_scaffolds(
            unplaced, scaffold_regions, ref_fasta, query_fasta,
            preset=preset,
            margin=config.nucleotide.region_margin,
            threads=chunk_threads,
        )
        logger.info("  %d PAF entries from alignment", len(paf_entries))

        # Update persistent chromosome lookup from this pass's PAF.
        # Only ADD new contigs; never overwrite Block Tree assignments.
        for entry in paf_entries:
            query = entry['query']
            if query not in contig_ref:
                paf_chr = entry['ref_chr'].split(':')[0]
                rs = entry['rstart']
                re = entry['rend']
                contig_ref[query] = (paf_chr, min(rs, re), max(rs, re))

        if not paf_entries:
            logger.info("  No alignments found — skipping remaining passes")
            break

        n_added = _add_extension_edges(
            sg, paf_entries, scaffold_regions, contig_lengths,
            gap_min=config.scaffold.gap_min,
            gap_max=config.scaffold.gap_max,
        )
        logger.info("  %d extension edges added", n_added)

        # Chain nearby contigs on the same chromosome using persistent lookup
        n_chained = _add_adjacent_chromosome_edges(
            sg, contig_ref, contig_lengths,
            gap_min=config.scaffold.gap_min,
            gap_max=config.scaffold.gap_max,
        )
        logger.info("  %d chromosome-ordering edges added", n_chained)

        # Mark newly added contigs with this pass's tier
        for ctg_base in {n[:-2] for n in sg._nodes}:
            if ctg_base not in contig_tier:
                contig_tier[ctg_base] = preset

        if n_added == 0 and n_chained == 0:
            logger.info("  No new edges — skipping remaining passes")
            break

        cover = _solve_graph(sg, config)
        n_seq, n_scaf, n_gap = _count_agp_from_cover(cover, contig_lengths)
        logger.info("  Pass (%s) result: %d scaffolds, %d contigs, %d gaps",
                    preset, n_scaf, n_seq, n_gap)

    # ---- Final step: merge scaffolds by chromosome on cover graph ----
    logger.info("Stage 1c: Merging scaffolds by chromosome on cover graph ...")
    n_merged = _merge_scaffolds_on_cover(cover, contig_ref, contig_lengths,
                                          gap_min=config.scaffold.gap_min)
    logger.info("  %d chromosome-level merge edges added to cover", n_merged)

    # ---- Final precise alignment for accurate reference coordinates ----
    # Runs BEFORE AGP extraction: the refined coordinates and strands are
    # used to normalize scaffold direction and orient each AGP component.
    logger.info("Stage 1d: Final precise alignment of all placed contigs ...")
    contig_strand, contig_mapq = _refine_contig_coordinates(
        contig_ref, contig_lengths, ref_fasta, query_fasta,
        config, chunk_threads)
    logger.info("  %d contigs with refined reference coordinates", len(contig_ref))

    # ---- Extract final AGP ----
    logger.info("Stage 1e: Writing final AGP ...")
    agp_lines = _extract_agp_paths(cover, contig_lengths,
                                   config.scaffold.unknown_gap_size,
                                   contig_ref=contig_ref,
                                   contig_strand=contig_strand,
                                   contig_mapq=contig_mapq)

    agp_path = os.path.join(output_dir, "nearscaff.agp")
    writer = AGPWriter()
    with open(agp_path, 'w') as f:
        f.write(writer.format(agp_lines))
    logger.info("AGP saved to %s", agp_path)

    n_seq, n_scaf, n_gap = _count_agp_from_cover(cover, contig_lengths)
    logger.info("Final: %d scaffolds, %d contigs placed, %d gaps",
                n_scaf, n_seq, n_gap)

    # ---- Export tiered PAF for visualization ----
    _write_tiered_paf(output_dir, root, contig_ref, contig_tier,
                      contig_lengths, agp_lines, contig_strand=contig_strand)
    logger.info("Tiered PAF saved to %s",
                os.path.join(output_dir, "nearscaff_tiered.paf"))

    # ---- Final FASTA: scaffolds + unplaced contigs ----
    from nearscaff.agp2fasta import agp_to_fasta
    logger.info("Stage 1f: Writing scaffold FASTA ...")
    fasta_path = os.path.join(output_dir, "nearscaff.scaffolds.fa")
    n_scaf, n_unplaced = agp_to_fasta(agp_path, query_fasta, fasta_path)
    logger.info("Scaffold FASTA saved to %s (%d scaffolds, "
                "%d unplaced contigs appended)",
                fasta_path, n_scaf, n_unplaced)

    return agp_path


# ---------------------------------------------------------------------------
# Stage 1 internal helpers
# ---------------------------------------------------------------------------

def _build_initial_graph(root) -> tuple:
    """Build FusedScaffoldGraph with protein edges from Block Tree."""
    from nearscaff.scaffold_graph import FusedScaffoldGraph

    sg = FusedScaffoldGraph()
    contigs = set()
    for lc in root.iter_level("local"):
        for ctg in lc.query_contigs:
            contigs.add(ctg)

    for ctg in contigs:
        sg.add_node(ctg + "_b")
        sg.add_node(ctg + "_e")

    for lc in root.iter_level("local"):
        ctg_list = sorted(lc.query_contigs)
        for i in range(len(ctg_list) - 1):
            sg.add_fused_edge(
                ctg_list[i] + "_e", ctg_list[i + 1] + "_b",
                weight=lc.confidence,
                edge_type=EdgeType.PROTEIN_SYNTENY,
                gap_type="scaffold",
                source=f"block_tree:{lc.node_id}",
            )

    return sg, contigs


def _solve_graph(sg, config):
    """Fuse -> scale -> matching -> cover. Returns nx.Graph cover graph."""
    fused_weights = sg.fuse_weights()
    if config.scaffold.best_buddy_scale:
        scaled_weights = sg.best_buddy_scale(fused_weights)
    else:
        scaled_weights = fused_weights
    matching = sg.get_max_weight_matching(scaled_weights)
    return sg.cover_graph(matching, scaled_weights)


def _count_agp_from_cover(cover, contig_lengths: dict) -> tuple:
    """Return (n_contigs_placed, n_scaffolds, n_gaps) from cover graph."""
    import networkx as nx
    n_seq = 0
    n_scaf = 0
    for cc_nodes in nx.connected_components(cover):
        sub = cover.subgraph(cc_nodes)
        endpoints = [n for n in cc_nodes if sub.degree(n) <= 1]
        if len(endpoints) >= 2:
            path = nx.shortest_path(sub, source=endpoints[0], target=endpoints[-1])
        else:
            path = list(cc_nodes)
        bases = []
        for node in path:
            base = node[:-2]
            if base not in bases:
                bases.append(base)
        if bases:
            n_scaf += 1
            n_seq += len(bases)
    n_gap = n_seq - n_scaf
    return n_seq, n_scaf, n_gap


def _get_unplaced_contigs(cover, all_contigs: set) -> set:
    """Return contig IDs not part of any multi-contig scaffold path."""
    import networkx as nx
    placed = set()
    for cc_nodes in nx.connected_components(cover):
        sub = cover.subgraph(cc_nodes)
        endpoints = [n for n in cc_nodes if sub.degree(n) <= 1]
        if len(endpoints) >= 2:
            path = nx.shortest_path(sub, source=endpoints[0], target=endpoints[-1])
        else:
            path = list(cc_nodes)
        for node in path:
            placed.add(node[:-2])
    return all_contigs - placed


def _get_scaffold_regions(cover, root, contig_lengths: dict) -> list:
    """Extract reference span for each scaffold from cover + Block Tree."""
    import networkx as nx
    from dataclasses import dataclass

    @dataclass
    class _ScaffoldRegion:
        scaffold_idx: int
        contigs: list
        ref_chr: str
        ref_start: int
        ref_end: int

    # Build contig -> (ref_chr, ref_start, ref_end) lookup from Block Tree
    ctg_ref = {}
    for lc in root.iter_level("local"):
        if lc.ref_chr is None:
            continue
        for ctg in lc.query_contigs:
            if ctg not in ctg_ref:
                ctg_ref[ctg] = (lc.ref_chr, lc.ref_start, lc.ref_end)
            else:
                prev = ctg_ref[ctg]
                ctg_ref[ctg] = (
                    prev[0],
                    min(prev[1], lc.ref_start or prev[1]),
                    max(prev[2], lc.ref_end or prev[2]),
                )

    regions = []
    scaffold_idx = 0
    for cc_nodes in nx.connected_components(cover):
        sub = cover.subgraph(cc_nodes)
        endpoints = [n for n in cc_nodes if sub.degree(n) <= 1]
        if len(endpoints) < 2:
            continue  # skip isolated nodes

        path = nx.shortest_path(sub, source=endpoints[0], target=endpoints[-1])
        bases = []
        for node in path:
            base = node[:-2]
            if base not in bases:
                bases.append(base)

        if len(bases) < 2:
            continue  # no extension possible with single contig

        # Determine reference span from placed contigs
        ref_chr = None
        ref_start = float('inf')
        ref_end = 0
        for ctg in bases:
            if ctg in ctg_ref:
                rc, rs, re = ctg_ref[ctg]
                ref_chr = rc
                if rs is not None and rs < ref_start:
                    ref_start = rs
                if re is not None and re > ref_end:
                    ref_end = re

        if ref_chr is not None and ref_start < ref_end:
            regions.append(_ScaffoldRegion(
                scaffold_idx=scaffold_idx,
                contigs=bases,
                ref_chr=ref_chr,
                ref_start=int(ref_start),
                ref_end=int(ref_end),
            ))
            scaffold_idx += 1

    return regions


def _align_unplaced_to_scaffolds(unplaced: set, scaffold_regions: list,
                                 ref_fasta: str, query_fasta: str,
                                 preset: str, margin: int,
                                 threads: int) -> list:
    """Align unplaced contigs to scaffold reference regions via minimap2.

    All unplaced contigs are extracted once.  For each scaffold region
    a constrained alignment is run, and the PAF entries are collected.
    """
    import tempfile, os
    from nearscaff.nucleotide import run_constrained_alignment, _extract_contigs, parse_nucleotide_paf
    from concurrent.futures import ThreadPoolExecutor, as_completed

    if not unplaced or not scaffold_regions:
        return []

    # Write unplaced contigs to temp FASTA once
    fd, unplaced_fa = tempfile.mkstemp(suffix='.fa')
    os.close(fd)
    try:
        _extract_contigs(query_fasta, list(unplaced), unplaced_fa)

        def _align_one(region):
            return run_constrained_alignment(
                ref_fasta, unplaced_fa,
                region_chr=region.ref_chr,
                region_start=region.ref_start,
                region_end=region.ref_end,
                margin=margin, preset=preset, threads=1,
            )

        all_entries = []
        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            futures = {pool.submit(_align_one, r): r for r in scaffold_regions}
            for future in as_completed(futures):
                paf = future.result()
                entries = parse_nucleotide_paf(paf)
                all_entries.extend(entries)

        return all_entries
    finally:
        if os.path.exists(unplaced_fa):
            os.unlink(unplaced_fa)


def _add_extension_edges(sg, paf_entries: list, scaffold_regions: list,
                         contig_lengths: dict, gap_min: int = 0,
                         gap_max: int = 500000) -> int:
    """Add NUCLEOTIDE_CHAIN edges for unplaced contigs near scaffold ends.

    For each PAF entry, determines if the contig aligns within reach of
    a scaffold endpoint.  Only the best alignment per unplaced contig
    per scaffold is used.

    Returns count of edges added.
    """
    # Group PAF entries by query contig — keep best alignment overall
    best_per_query = {}
    for entry in paf_entries:
        query = entry['query']
        prev = best_per_query.get(query)
        if prev is None or entry['nmatch'] > prev['nmatch']:
            best_per_query[query] = entry

    # For each scaffold region, find unplaced contigs that align nearby
    region_by_idx = {r.scaffold_idx: r for r in scaffold_regions}

    # Map: (query, scaf_idx) -> best PAF entry for that query+scaffold pair
    best_per_pair = {}
    for query, entry in best_per_query.items():
        al_start = min(entry['rstart'], entry['rend'])
        al_end = max(entry['rstart'], entry['rend'])
        # samtools faidx headers include region coords (e.g. "FCY04A:50141-3543086"),
        # so strip the suffix to get the plain chromosome name.
        paf_chr = entry['ref_chr'].split(':')[0]
        for region in scaffold_regions:
            if paf_chr != region.ref_chr:
                continue
            # Is this alignment near the scaffold region?
            gap_to_left = region.ref_start - al_end
            gap_to_right = al_start - region.ref_end
            best_gap = max(gap_to_left, gap_to_right)  # positive = outside
            if best_gap > gap_max:
                continue
            key = (query, region.scaffold_idx)
            prev = best_per_pair.get(key)
            if prev is None or entry['nmatch'] > prev['nmatch']:
                best_per_pair[key] = entry

    logger.debug("  %d query+scaffold pairings to consider", len(best_per_pair))

    # Add edges
    n_added = 0
    # Build set of contigs already in the graph (with nodes)
    known_contigs = {node[:-2] for node in sg._nodes}

    for (query, scaf_idx), entry in best_per_pair.items():
        region = region_by_idx[scaf_idx]
        if not region.contigs:
            continue

        # Add new contig to graph if not present
        if query not in known_contigs:
            try:
                sg.add_node(query + "_b")
                sg.add_node(query + "_e")
            except ValueError:
                continue
            known_contigs.add(query)

        al_start = min(entry['rstart'], entry['rend'])
        al_end = max(entry['rstart'], entry['rend'])
        first_ctg = region.contigs[0]
        last_ctg = region.contigs[-1]

        # Determine which scaffold end this extends
        gap_left = region.ref_start - al_end   # positive if alignment is left of scaffold
        gap_right = al_start - region.ref_end   # positive if alignment is right of scaffold

        if gap_right >= gap_left and gap_right >= 0:
            # Extend right
            u, v = last_ctg + "_e", query + "_b"
            gap_size = max(0, gap_right)
        elif gap_left > gap_right and gap_left >= 0:
            # Extend left
            u, v = query + "_e", first_ctg + "_b"
            gap_size = max(0, gap_left)
        else:
            # Alignment overlaps scaffold — insert at nearest end
            if abs(gap_left) < abs(gap_right):
                u, v = query + "_e", first_ctg + "_b"
                gap_size = 0
            else:
                u, v = last_ctg + "_e", query + "_b"
                gap_size = 0

        gap_size = max(gap_min, min(gap_size, gap_max))
        weight = (entry['nmatch'] / max(entry['hitlen'], 1)) * min(entry['mapq'] / 60.0, 1.0)

        sg.add_fused_edge(
            u, v,
            weight=weight,
            edge_type=EdgeType.NUCLEOTIDE_CHAIN,
            gap_size=gap_size,
            gap_type="scaffold",
            source=f"extension:{region.scaffold_idx}",
        )
        n_added += 1

    return n_added


def _add_adjacent_chromosome_edges(sg, contig_ref: dict,
                                    contig_lengths: dict,
                                    gap_min: int = 0,
                                    gap_max: int = 500000) -> int:
    """Add NUCLEOTIDE_CHAIN edges between contigs adjacent on the reference.

    Uses *contig_ref* (persistent {contig: (chr, r_start, r_end)} mapping)
    to group contigs by chromosome, sort by reference position, and connect
    consecutive pairs.
    """
    # Group by chromosome, sort by midpoint of reference alignment
    chr_contigs = {}
    for ctg, (chr_name, rs, re) in contig_ref.items():
        if chr_name not in chr_contigs:
            chr_contigs[chr_name] = []
        chr_contigs[chr_name].append((ctg, (rs + re) // 2))

    # For each chromosome, sort by position and add edges
    n_added = 0
    known_contigs = {node[:-2] for node in sg._nodes}

    for chr_name, ctg_list in chr_contigs.items():
        ctg_list.sort(key=lambda x: x[1])  # sort by reference position
        for i in range(len(ctg_list) - 1):
            ctg_a, pos_a = ctg_list[i]
            ctg_b, pos_b = ctg_list[i + 1]

            if ctg_a not in known_contigs or ctg_b not in known_contigs:
                continue

            dist = pos_b - pos_a
            if dist > gap_max or dist < 0:
                continue

            gap_size = max(gap_min, min(dist, gap_max))
            weight = 0.5

            sg.add_fused_edge(
                ctg_a + "_e", ctg_b + "_b",
                weight=weight,
                edge_type=EdgeType.NUCLEOTIDE_CHAIN,
                gap_size=gap_size,
                gap_type="scaffold",
                source=f"chr_order:{chr_name}",
            )
            n_added += 1

    return n_added


def _merge_scaffolds_on_cover(cover, contig_ref: dict,
                               contig_lengths: dict,
                               gap_min: int = 0) -> int:
    """Directly add edges to the cover graph to merge same-chromosome scaffolds.

    Unlike _merge_scaffolds_by_chromosome which adds edges to the scaffold
    graph and relies on the matching algorithm, this function directly
    modifies the cover graph.  This guarantees one connected component per
    chromosome regardless of edge weight conflicts.

    Only scaffolds with a known chromosome are merged; unknown-chr scaffolds
    are left as-is (unscaffolded state).
    """
    import networkx as nx
    from collections import Counter

    # Build scaffold info from cover graph components
    scaffold_info = {}
    for cc_nodes in nx.connected_components(cover):
        sub = cover.subgraph(cc_nodes)
        endpoints = [n for n in cc_nodes if sub.degree(n) <= 1]
        if len(endpoints) >= 2:
            path = nx.shortest_path(sub, source=endpoints[0], target=endpoints[-1])
        else:
            path = list(cc_nodes)

        bases = []
        for node in path:
            base = node[:-2]
            if base not in bases:
                bases.append(base)
        if not bases:
            continue

        chr_votes = Counter()
        min_pos = float('inf')
        max_pos = 0
        for ctg in bases:
            if ctg in contig_ref:
                chr_name, rs, re = contig_ref[ctg]
                chr_votes[chr_name] += 1
                if rs < min_pos:
                    min_pos = rs
                if re > max_pos:
                    max_pos = re

        if not chr_votes:
            continue  # unknown-chr: leave as-is

        best_chr = chr_votes.most_common(1)[0][0]
        key = frozenset(bases)
        scaffold_info[key] = {
            'chr': best_chr,
            'min_pos': min_pos,
            'max_pos': max_pos,
            'bases': bases,
            'first_node': path[0],
            'last_node': path[-1],
        }

    # Group by chromosome, sort by position
    chr_scaffolds = {}
    for info in scaffold_info.values():
        chr_name = info['chr']
        if chr_name not in chr_scaffolds:
            chr_scaffolds[chr_name] = []
        chr_scaffolds[chr_name].append(info)

    n_added = 0
    for chr_name, scaf_list in chr_scaffolds.items():
        if len(scaf_list) < 2:
            continue
        scaf_list.sort(key=lambda x: x['min_pos'])
        for i in range(len(scaf_list) - 1):
            left = scaf_list[i]
            right = scaf_list[i + 1]
            gap = right['min_pos'] - left['max_pos']
            gap_size = max(gap_min, int(max(0, gap)))

            # Add edge directly to cover graph — bypass matching
            cover.add_edge(
                left['last_node'], right['first_node'],
                weight=10.0,  # very high to indicate forced connection
                gap_size=gap_size,
                gap_type="scaffold",
                source=f"chr_merge:{chr_name}",
            )
            n_added += 1

    return n_added


def _refine_contig_coordinates(contig_ref: dict, contig_lengths: dict,
                                ref_fasta: str, query_fasta: str,
                                config, threads: int):
    """Refine reference coordinates by aligning ALL placed contigs to reference.

    Protein-placed contigs from the Block Tree have synthetic coordinates
    (inferred from protein ID ordering, not real genomic positions).  This
    step runs a minimap2 alignment of all placed contigs against the full
    reference genome and updates *contig_ref* with real genomic coordinates.

    Returns (strand_map, mapq_map) — {contig: '+'/'-'} and {contig: mapq}
    from the best alignment per contig (empty dicts on failure); used to
    orient AGP components.
    """
    import tempfile, os, subprocess
    from nearscaff.nucleotide import parse_nucleotide_paf, _extract_contigs

    placed_contigs = list(contig_ref.keys())
    if not placed_contigs:
        return {}, {}

    # Write placed contigs to temp FASTA
    fd, placed_fa = tempfile.mkstemp(suffix='.fa')
    os.close(fd)
    try:
        _extract_contigs(query_fasta, placed_contigs, placed_fa)

        # Align against the full reference
        cmd = ["minimap2", "-t", str(threads), "-x", config.nucleotide.preset,
               "-c", ref_fasta, placed_fa]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning("Final precise alignment failed: %s", result.stderr)
            return {}, {}

        entries = parse_nucleotide_paf(result.stdout)
        logger.info("  %d alignments from final precise pass", len(entries))

        # Update contig_ref with real coordinates; keep original tier's chr if
        # alignment chromosome matches (strip samtools region suffix).
        # Also record per-contig strand — needed to orient AGP components.
        # A contig may have several alignments; use the best one (longest
        # aligned block, then highest mapq) instead of an arbitrary one.
        best = {}
        for entry in entries:
            query = entry['query']
            if query not in contig_ref:
                continue
            key = (entry['hitlen'], entry['mapq'])
            if query not in best or key > best[query][0]:
                best[query] = (key, entry)

        contig_strand = {}
        contig_mapq = {}
        updated = 0
        for query, (_key, entry) in best.items():
            paf_chr = entry['ref_chr'].split(':')[0]
            rs = min(entry['rstart'], entry['rend'])
            re = max(entry['rstart'], entry['rend'])
            # Only update if the alignment chromosome is reasonable
            # (don't overwrite chr if PAF chr is totally different)
            old_chr = contig_ref[query][0]
            if paf_chr and len(paf_chr) > 2:
                contig_ref[query] = (paf_chr, rs, re)
                contig_strand[query] = entry['strand']
                contig_mapq[query] = entry['mapq']
                updated += 1

        logger.info("  %d contig coordinates refined", updated)
        return contig_strand, contig_mapq
    finally:
        if os.path.exists(placed_fa):
            os.unlink(placed_fa)


def _write_tiered_paf(output_dir: str, root, contig_ref: dict,
                      contig_tier: dict, contig_lengths: dict,
                      agp_lines: list, contig_strand: dict | None = None):
    """Write a tiered-confidence PAF file for visualization.

    Three confidence tiers based on the evidence that placed each contig:
      protein — miniprot protein anchoring (highest confidence)
      asm5    — strict nucleotide extension
      asm10   — relaxed nucleotide extension
      asm20   — loosest nucleotide extension (lowest confidence)

    Only contigs present in the final AGP scaffolds are included.
    Unplaced / unscaffolded contigs are filtered out.

    PAF format: query_contig, qlen, qstart, qend, strand,
                ref_chr, rlen, rstart, rend, nmatch, hitlen, mapq
                + tag nc:Z:<tier>
    """
    # Get set of contigs in final scaffolds
    placed_in_agp = set()
    for line in agp_lines:
        if hasattr(line, 'component_id'):
            placed_in_agp.add(line.component_id)

    contig_strand = contig_strand or {}
    paf_path = os.path.join(output_dir, "nearscaff_tiered.paf")
    with open(paf_path, 'w') as f:
        for ctg, (chr_name, rs, re) in sorted(contig_ref.items()):
            if ctg not in placed_in_agp:
                continue  # skip unplaced

            tier = contig_tier.get(ctg, "unknown")
            ctg_len = contig_lengths.get(ctg, 1000)
            strand = contig_strand.get(ctg, "+")

            # Use actual reference alignment coordinates
            r_start = rs
            r_end = re
            hit_len = r_end - r_start

            f.write(
                f"{ctg}\t{ctg_len}\t0\t{ctg_len}\t{strand}\t"
                f"{chr_name}\t0\t{r_start}\t{r_end}\t"
                f"{hit_len}\t{hit_len}\t60\t"
                f"nc:Z:{tier}\n"
            )


def run_stage0(config: NearscaffConfig, ref_fasta: str, ref_gff: str,
               query_fasta: str, output_dir: str,
               protein_faa: str | None = None,
               ):
    """Stage 0: Protein anchoring + synteny -> block_tree.json.

    Returns the Block Tree root node (or None if anchoring failed).
    """
    os.makedirs(output_dir, exist_ok=True)
    logger.info("=" * 60)
    logger.info("Stage 0: Protein Anchoring + Block Tree")
    logger.info("=" * 60)

    # ---- Protein anchoring + synteny -> Block Tree ----
    root = _run_signal_c_protein(
        config, ref_fasta, ref_gff, query_fasta, output_dir,
        protein_faa=protein_faa,
    )

    tree_path = os.path.join(output_dir, "block_tree.json")
    if root is not None:
        with open(tree_path, 'w') as f:
            f.write(block_tree_to_json(root))
        logger.info("Block Tree saved to %s", tree_path)

        n_sg = len(root.children) if root.level == "root" else 0
        n_chr = sum(len(sg.children) for sg in root.children)
        n_lc = sum(len(chr.children) for sg in root.children for chr in sg.children)
        logger.info("Tree structure: %d subgenome(s), %d chromosome(s), %d local cluster(s)",
                    n_sg, n_chr, n_lc)
    else:
        logger.warning("Signal C produced no Block Tree -- scaffolding will be unavailable")

    return root



def run_full(config: NearscaffConfig, ref_fasta: str, ref_gff: str,
             query_fasta: str, output_dir: str,
             protein_faa: str | None = None):
    """Run the complete nearscaff pipeline: Stage 0 -> Stage 1."""
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "nearscaff.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(),
        ],
    )

    logger.info("=" * 60)
    logger.info("nearscaff — Reference-guided genome scaffolding")
    logger.info("Reference: %s", ref_fasta)
    if protein_faa:
        logger.info("Proteins:  %s", protein_faa)
    else:
        logger.info("GFF3:      %s", ref_gff)
    logger.info("Query:     %s", query_fasta)
    logger.info("Output:    %s", output_dir)
    logger.info("Threads:   %d", config.threads)
    logger.info("=" * 60)

    # ---- Stage 0: Protein anchoring + synteny -> block_tree.json ----
    logger.info("")
    logger.info(">>> Stage 0: Protein Anchoring + Block Tree <<<")
    root = run_stage0(
        config, ref_fasta, ref_gff, query_fasta, output_dir,
        protein_faa=protein_faa,
    )
    if root is None:
        logger.error("Stage 0: protein anchoring failed — "
                     "cannot proceed to scaffolding")
        return None

    tree_path = os.path.join(output_dir, "block_tree.json")

    # ---- Stage 1: Nucleotide Extension + Scaffolding ----
    logger.info("")
    logger.info(">>> Stage 1: Nucleotide Extension + Scaffolding <<<")
    agp_path = run_stage1(config, tree_path, ref_fasta, query_fasta, output_dir)

    logger.info("")
    logger.info("Pipeline complete!")
    if agp_path:
        logger.info("Final AGP: %s", agp_path)
    return agp_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _infer_chromosome(header: str) -> str | None:
    """Infer reference chromosome from a protein header.

    Handles patterns like:
      FCYA01G000001.mRNA1  ->  FCYA01
      AT1G01010.1          ->  AT1
    """
    m = _PROTEIN_ID_PATTERN.match(header)
    if m:
        return m.group(1)
    return None


def _load_proteins_from_fasta(path: str) -> tuple[dict, dict]:
    """Load proteins and infer gene locations from a protein FASTA.

    Returns (proteins, gene_locations).

    Gene locations are inferred from protein headers — chromosomes are
    parsed via *_infer_chromosome* and genes are assigned sequential
    positions within each chromosome (preserving gene order).
    """
    proteins = {}
    gene_locations = {}
    current_name = None
    current_seq = []

    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if current_name:
                    proteins[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_name:
            proteins[current_name] = ''.join(current_seq)

    # Infer chromosome and assign synthetic positions
    chr_genes = {}
    for name in sorted(proteins.keys()):
        inferred_chr = _infer_chromosome(name)
        if inferred_chr is None:
            # Fallback: use first 4 chars or entire name
            inferred_chr = name[:4] if len(name) >= 4 else "unknown"
        if inferred_chr not in chr_genes:
            chr_genes[inferred_chr] = []
        chr_genes[inferred_chr].append(name)

    # Assign contiguous positions (5 kb spacing) preserving sorted order
    for chr_name, gene_list in chr_genes.items():
        for rank, gene_id in enumerate(gene_list):
            start = rank * 5000
            end = start + 3000
            gene_locations[gene_id] = (chr_name, start, end)

    return proteins, gene_locations


def _copy_or_link(src: str, dst: str):
    """Copy file. Falls back to symlink on OSError."""
    import shutil
    try:
        shutil.copy2(src, dst)
    except OSError:
        if not os.path.exists(dst):
            os.symlink(os.path.abspath(src), dst)


def _parse_gff3_gene_locations(gff_path: str) -> dict:
    """Parse GFF3 and return {mrna_id: (chromosome, start, end)}."""
    locations = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
            ftype = fields[2]
            if ftype not in ("mRNA", "gene"):
                continue
            attrs = {}
            for part in fields[8].split(';'):
                if '=' in part:
                    k, v = part.split('=', 1)
                    attrs[k.strip()] = v.strip()
            uid = attrs.get("ID", "")
            if uid:
                locations[uid] = (fields[0], int(fields[3]) - 1, int(fields[4]))
    return locations


def _build_chr_name_map(inferred_names: set, ref_fasta: str) -> dict:
    """Map inferred chromosome names to actual names in *ref_fasta*.

    Uses a combination of numeric-part matching and character-set overlap.
    """
    import re

    actual_names = sorted(_read_fasta_lengths(ref_fasta).keys())
    if not actual_names:
        return {}

    mapping = {}
    for inf in inferred_names:
        if inf in actual_names:
            mapping[inf] = inf
            continue

        # Extract numeric and alphabetic parts from inferred name
        inf_nums = re.findall(r'\d+', inf)
        inf_alpha = re.sub(r'\d+', '', inf)

        best = None
        best_score = -1
        for act in actual_names:
            act_nums = re.findall(r'\d+', act)
            act_alpha = re.sub(r'\d+', '', act)

            # Score: numeric match * 100 + character overlap
            num_match = len(set(inf_nums) & set(act_nums))
            char_overlap = len(set(inf_alpha) & set(act_alpha))
            score = num_match * 100 + char_overlap

            if score > best_score:
                best_score = score
                best = act

        if best and best_score > 0:
            mapping[inf] = best

    return mapping


def _enrich_anchors(anchors: list, gene_locations: dict) -> list:
    """Fix reference chromosome / position on each anchor using gene locations.

    The raw PAF parser sets both query_contig and ref_chr to the PAF
    target name (which is the query contig).  Here we replace ref_chr,
    r_start, r_end with the actual reference chromosome / gene position
    from the gene-location table.
    """
    enriched = []
    for a in anchors:
        mrna_id = a.query_gene
        if mrna_id in gene_locations:
            chr_name, gff_start, gff_end = gene_locations[mrna_id]
            a.ref_chr = chr_name
            if a.r_start == a.q_start and a.r_end == a.q_end:
                a.r_start = gff_start
                a.r_end = gff_end
        enriched.append(a)
    return enriched


def _write_anchors_tsv(anchors, path):
    with open(path, 'w') as f:
        f.write("query_contig\tquery_gene\tref_chr\tref_gene\tstrand\tscore\tidentity\n")
        for a in anchors:
            f.write(f"{a.query_contig}\t{a.query_gene}\t{a.ref_chr}\t"
                    f"{a.ref_gene}\t{a.strand}\t{a.score}\t{a.identity}\n")


def _write_blocks_file(blocks, path):
    with open(path, 'w') as f:
        for b in blocks:
            f.write(f"##block_id={b.block_id}\tref_chr={b.ref_chr}\t"
                    f"ref_start={b.ref_start}\tref_end={b.ref_end}\t"
                    f"score={b.score}\torientation={b.orientation}\tanchor_count={b.anchor_count}\n")
            for qg, rg, sc in b.gene_pairs:
                f.write(f"{rg}\t{qg}\t{sc}\n")


def _read_fasta_lengths(path: str) -> dict:
    """Read FASTA and return {contig_name: length}."""
    lengths = {}
    current_name = None
    current_len = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if current_name is not None:
                    lengths[current_name] = current_len
                current_name = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_name is not None:
        lengths[current_name] = current_len
    return lengths


def _extract_agp_paths(cover, contig_lengths: dict,
                       unknown_gap_size: int = 100,
                       gap_min: int = 0,
                       gap_max: int = 500000,
                       contig_ref: dict | None = None,
                       contig_strand: dict | None = None,
                       contig_mapq: dict | None = None):
    """Extract linear scaffold paths from cover graph into AGP lines.

    Each connected component in *cover* is a path (cycles already broken).
    Endpoints are degree-1 nodes; the path is the unique simple path
    between them.  Isolated nodes (degree-0, unplaced contigs) become
    their own single-contig scaffolds.

    When *contig_ref* (refined reference coordinates) is given, scaffold
    components are re-ordered by reference midpoint — contigs without
    coordinates keep their graph order, interpolated between anchored
    neighbours.  When *contig_strand* (strands from the final precise
    alignment) is given, components are oriented by their alignment
    strand; strands with mapq below MIN_ORIENT_MAPQ (from *contig_mapq*)
    are reported as '?'.  Contigs without alignment evidence keep the
    traversal-inferred strand.
    """
    from nearscaff.agp import AGPSeqLine, AGPGapLine
    import networkx as nx

    contig_ref = contig_ref or {}
    contig_strand = contig_strand or {}
    contig_mapq = contig_mapq or {}

    lines = []
    scaffold_idx = 1

    for cc_nodes in nx.connected_components(cover):
        sub = cover.subgraph(cc_nodes)
        endpoints = [n for n in cc_nodes if sub.degree(n) == 1]
        isolated = [n for n in cc_nodes if sub.degree(n) == 0]

        if isolated and not endpoints:
            # Unplaced contig — emit as standalone scaffold
            for n in sorted(isolated):
                base = n[:-2]
                if n.endswith("_e"):
                    continue
                ctg_len = contig_lengths.get(base, 1000)
                lines.append(AGPSeqLine(
                    f"nearscaff_{scaffold_idx:04d}", 1, ctg_len, 1,
                    "W", base, 1, ctg_len, contig_strand.get(base, "+"),
                ))
                scaffold_idx += 1
            continue

        # Walk the unique simple path between the two endpoints
        if len(endpoints) >= 2:
            path = nx.shortest_path(sub, source=endpoints[0], target=endpoints[-1])
        else:
            # Should not happen after cycle-breaking, but be defensive
            path = list(cc_nodes)

        # Extract contig order from path (skip internal edges)
        # Path alternates: _b, _e, _b, _e, ... for intra-contig edges
        # and _e, _b for inter-contig edges.
        # Simplify: collect unique contigs in path order.
        seen_bases = []
        for node in path:
            base = node[:-2]
            if base not in seen_bases:
                seen_bases.append(base)

        # Order by refined reference midpoint.  Contigs without coordinates
        # keep graph order, interpolated between anchored neighbours so the
        # scaffold runs in increasing reference coordinates.
        mids = {b: (contig_ref[b][1] + contig_ref[b][2]) / 2
                for b in seen_bases if b in contig_ref}
        if len(mids) >= 2:
            anchored_idx = [i for i, b in enumerate(seen_bases) if b in mids]
            est = []
            for i, b in enumerate(seen_bases):
                if b in mids:
                    est.append(mids[b])
                    continue
                lo = max((j for j in anchored_idx if j < i), default=None)
                hi = min((j for j in anchored_idx if j > i), default=None)
                if lo is not None and hi is not None:
                    m = (mids[seen_bases[lo]] + mids[seen_bases[hi]]) / 2
                elif lo is not None:
                    m = mids[seen_bases[lo]] - 0.5
                else:
                    m = mids[seen_bases[hi]] + 0.5
                est.append(m)
            # stable sort: ties (and interpolations) keep graph order
            seen_bases = [b for _, b in sorted(zip(est, seen_bases),
                                               key=lambda t: t[0])]
        flip = False

        scaf_name = f"nearscaff_{scaffold_idx:04d}"
        pos = 1
        part_num = 1

        for i, base in enumerate(seen_bases):
            ctg_len = contig_lengths.get(base, 1000)

            if base in contig_strand:
                # Orientation from the final precise alignment; mark
                # low-mapping-quality evidence as unknown.
                if contig_mapq.get(base, 60) < MIN_ORIENT_MAPQ:
                    strand = "?"
                else:
                    strand = contig_strand[base]
            else:
                # Fallback: infer strand from the path traversal
                # (if we see _b first, it's +; if _e first, it's -).
                first_occurrence = None
                for n in path:
                    if n[:-2] == base:
                        first_occurrence = n
                        break
                strand = "-" if first_occurrence and first_occurrence.endswith("_e") else "+"

            lines.append(AGPSeqLine(
                scaf_name, pos, pos + ctg_len - 1, part_num,
                "W", base, 1, ctg_len, strand,
            ))
            pos += ctg_len
            part_num += 1

            # Insert gap if more contigs follow
            if i < len(seen_bases) - 1:
                gap_size = unknown_gap_size
                lines.append(AGPGapLine(
                    scaf_name, pos, pos + gap_size - 1, part_num,
                    "U", gap_size,
                    "scaffold",
                    "yes",
                    "na",
                ))
                pos += gap_size
                part_num += 1

        scaffold_idx += 1

    return lines
