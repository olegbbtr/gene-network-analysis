"""Edge-level checks of a graph against curated TF->target pairs."""

import numpy as np
from scipy.sparse import csr_array

from .graph import edge_weights


def edge_metrics(
    adjacency: csr_array, sources: np.ndarray, targets: np.ndarray, signs: np.ndarray
) -> dict:
    """Compare a graph with reference pairs given as node indices.

    Only pairs whose regulator sends edges in this graph count. Chance rates assume
    each such regulator's edges point to uniformly random other genes.
    """
    out_degree = np.diff(adjacency.indptr)
    usable = out_degree[sources] > 0
    sources, targets, signs = sources[usable], targets[usable], signs[usable]
    regulators = np.unique(sources)
    possible = len(regulators) * (adjacency.shape[0] - 1)
    if not len(sources):
        return {"reference_pairs": 0}
    forward = edge_weights(adjacency, sources, targets)
    backward = edge_weights(adjacency, targets, sources)
    present = forward != 0
    rows = adjacency[regulators].tocsr()
    graph_edges = int(rows.nnz)
    chance_recall = graph_edges / possible
    chance_precision = len(sources) / possible
    any_edge = present | (backward != 0)
    result = {
        "reference_pairs": int(len(sources)),
        "reference_regulators": int(len(regulators)),
        "graph_edges_from_regulators": graph_edges,
        "recall": float(present.mean()),
        "chance_recall": float(chance_recall),
        "recall_enrichment": float(present.mean() / chance_recall),
        "direction": {
            "pairs_with_an_edge": int(any_edge.sum()),
            "forward_only": int((present & (backward == 0)).sum()),
            "reverse_only": int((~present & (backward != 0)).sum()),
            "both": int((present & (backward != 0)).sum()),
        },
        "recall_by_reference_sign": {
            name: float(present[signs == value].mean()) if (signs == value).any() else None
            for name, value in (("activation", 1), ("repression", -1), ("unknown", 0))
        },
        "precision_by_graph_sign": None,
        "sign_agreement": None,
    }
    negative = rows.data < 0
    if negative.any():
        result["precision_by_graph_sign"] = {
            name: {
                "edges": int(edges),
                "curated": int(curated),
                "precision": float(curated / edges) if edges else None,
                "enrichment": float(curated / edges / chance_precision) if edges else None,
            }
            for name, edges, curated in (
                ("positive", (~negative).sum(), (forward > 0).sum()),
                ("negative", negative.sum(), (forward < 0).sum()),
            )
        }
        signed = present & (signs != 0)
        if signed.any():
            graph_sign = np.sign(forward[signed])
            reference_sign = signs[signed]
            positive_share = float((graph_sign > 0).mean())
            activation_share = float((reference_sign > 0).mean())
            result["sign_agreement"] = {
                "pairs": int(signed.sum()),
                "agreement": float((graph_sign == reference_sign).mean()),
                "chance": positive_share * activation_share
                + (1 - positive_share) * (1 - activation_share),
            }
    return result
