"""Fused Scaffold Graph — multi-modal edge fusion for scaffolding.

Extends the RagTag ScaffoldGraph concept with typed edges from
protein+synteny, nucleotide alignment, Block Tree topology, and Hi-C.
"""
from collections import defaultdict
import networkx as nx
from nearscaff.types import EdgeType, EDGE_TYPE_WEIGHTS


class FusedScaffoldGraph:
    """Graph for scaffolding with multi-type edge fusion.

    Nodes: contig_b / contig_e (as in RagTag).
    Each contig is split into begin(_b) and end(_e) half-nodes.
    Edges between half-nodes represent adjacencies with typed evidence.
    """

    def __init__(self):
        self.graph = nx.MultiGraph()
        self._nodes = set()

    def add_node(self, name: str):
        """Add a contig half-node. Name must end with _b or _e."""
        if not name.endswith("_b") and not name.endswith("_e"):
            raise ValueError(f"Node name must end with '_b' or '_e': {name}")
        self._nodes.add(name)
        self.graph.add_node(name)

    def add_fused_edge(self, u: str, v: str, weight: float,
                       edge_type: EdgeType,
                       gap_size: int = 100,
                       gap_type: str = "scaffold",
                       source: str = "unknown"):
        """Add a typed edge between two half-nodes."""
        self.graph.add_edge(u, v,
            weight=weight,
            edge_type=edge_type,
            gap_size=gap_size,
            gap_type=gap_type,
            source=source,
        )

    def fuse_weights(self) -> dict:
        """Fuse multi-type edges: linear weighted sum per node pair.

        Returns: {(u, v): total_weight}
        """
        pair_weights = defaultdict(float)
        for u, v, data in self.graph.edges(data=True):
            edge_type = data['edge_type']
            type_weight = EDGE_TYPE_WEIGHTS.get(edge_type, 0.5)
            fused = data['weight'] * type_weight
            key = (u, v) if u < v else (v, u)
            pair_weights[key] += fused
        return dict(pair_weights)

    def best_buddy_scale(self, fused_weights: dict) -> dict:
        """Apply SALSA2 best-buddy weight scaling.

        w_scaled(u,v) = w(u,v) / max(max_incident(u), max_incident(v))
        If current edge is the max for both, use second-best as denominator.
        Reference: Ghurye et al. (2019) PLOS Comp Bio
        """
        # Incident max per node
        max_incident = defaultdict(float)
        incident_weights = defaultdict(list)
        for (u, v), w in fused_weights.items():
            incident_weights[u].append(w)
            incident_weights[v].append(w)
            if w > max_incident[u]:
                max_incident[u] = w
            if w > max_incident[v]:
                max_incident[v] = w

        scaled = {}
        for (u, v), w in fused_weights.items():
            best_alt = max(max_incident[u], max_incident[v])
            if best_alt == w:
                # Current edge IS the best for at least one node.
                # Find second-best among all incident edges.
                alt_u = sorted(incident_weights[u], reverse=True)
                alt_v = sorted(incident_weights[v], reverse=True)
                second_u = alt_u[1] if len(alt_u) > 1 else 0
                second_v = alt_v[1] if len(alt_v) > 1 else 0
                best_alt = max(second_u, second_v)
            if best_alt == 0:
                best_alt = 1.0
            scaled[(u, v)] = w / best_alt

        return scaled

    def get_max_weight_matching(self, weights: dict) -> set:
        """Solve maximum weight matching.

        Small graphs use Edmonds' Blossom algorithm (exact).  Large
        graphs fall back to a greedy heaviest-edge-first matching —
        Blossom is O(V^2 * E) and infeasible beyond ~10^4 edges.
        """
        if len(weights) <= 20000:
            G = nx.Graph()
            for (u, v), w in weights.items():
                G.add_edge(u, v, weight=w)
            return nx.max_weight_matching(G, maxcardinality=False, weight='weight')
        matched = set()
        matching = set()
        for (u, v), _w in sorted(weights.items(), key=lambda kv: kv[1],
                                 reverse=True):
            if u == v or u in matched or v in matched:
                continue
            matched.add(u)
            matched.add(v)
            matching.add((u, v))
        return matching

    def cover_graph(self, matching: set, weights: dict) -> nx.Graph:
        """Build cover graph from matching, eliminating cycles.

        Adds infinite-weight intra-contig edges (_b ↔ _e) to prevent
        contigs from being split, then detects and breaks cycles.
        """
        cover = nx.Graph()
        # Add intra-contig edges (infinite weight)
        bases = set()
        for node in self._nodes:
            bases.add(node[:-2])
        for base in bases:
            b_node = base + "_b"
            e_node = base + "_e"
            if b_node in self._nodes and e_node in self._nodes:
                cover.add_edge(b_node, e_node, weight=float('inf'))

        # Add matching edges
        for u, v in matching:
            key = (u, v) if (u, v) in weights else (v, u)
            w = weights.get(key, 0)
            cover.add_edge(u, v, weight=w)

        # Detect and break cycles
        for cc in list(nx.connected_components(cover)):
            sub = cover.subgraph(cc).copy()
            if sub.number_of_nodes() == sub.number_of_edges():
                edges = [(u, v, d['weight']) for u, v, d in sub.edges(data=True)
                         if d['weight'] != float('inf')]
                if edges:
                    to_remove = min(edges, key=lambda e: e[2])
                    cover.remove_edge(to_remove[0], to_remove[1])

        return cover
