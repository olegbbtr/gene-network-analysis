"""Graph transformations. Every step keeps the matrix shape, so node indices stay shared."""

import heapq

import numpy as np
from scipy.sparse import csr_array, vstack

from .graph import edge_weights


def _rebuild(adjacency: csr_array, rows, columns, weights) -> csr_array:
    keep = weights != 0
    result = csr_array(
        (weights[keep], (rows[keep], columns[keep])), shape=adjacency.shape, dtype=np.float64
    )
    result.sort_indices()
    return result


# Preprocessing: graph -> graph


def reverse_weights(adjacency: csr_array):
    """Edge arrays plus, for each edge A->B, the weight of B->A (0 when absent)."""
    edges = adjacency.tocoo()
    rows, columns, weights = edges.row, edges.col, edges.data
    return rows, columns, weights, edge_weights(adjacency, columns, rows)


def resolve_reciprocals(adjacency: csr_array) -> csr_array:
    """Same-sign pairs A->B, B->A keep the stronger edge, minus the weaker weight (2016 rule)."""
    rows, columns, weights, back = reverse_weights(adjacency)
    paired = weights * back > 0
    stronger = (np.abs(weights) > np.abs(back)) | (
        (np.abs(weights) == np.abs(back)) & (rows < columns)
    )
    result = np.where(paired, np.where(stronger, weights - back, 0.0), weights)
    return _rebuild(adjacency, rows, columns, result)


def orient_reciprocals(adjacency: csr_array, ratio: float = 2.0) -> csr_array:
    """Pairs A->B, B->A keep an edge only if it is at least `ratio` times the other one."""
    rows, columns, weights, back = reverse_weights(adjacency)
    keep = np.abs(weights) >= ratio * np.abs(back)
    return _rebuild(adjacency, rows, columns, np.where(keep, weights, 0.0))


def top_outgoing(adjacency: csr_array, k: int = 50) -> csr_array:
    """Keep each node's k strongest outgoing edges by absolute weight."""
    edges = adjacency.tocsr().tocoo()
    order = np.lexsort((-np.abs(edges.data), edges.row))
    starts = adjacency.indptr[edges.row[order]]
    rank = np.empty(edges.nnz, dtype=np.int64)
    rank[order] = np.arange(edges.nnz) - starts
    return _rebuild(adjacency, edges.row, edges.col, np.where(rank < k, edges.data, 0.0))


def keep_sign(adjacency: csr_array, sign: int = 1) -> csr_array:
    """Keep only activating (sign 1) or repressing (sign -1) edges."""
    edges = adjacency.tocoo()
    keep = np.sign(edges.data) == sign
    return _rebuild(adjacency, edges.row, edges.col, np.where(keep, edges.data, 0.0))


def drop_sinks(adjacency: csr_array) -> csr_array:
    """One pass: remove every edge into a node without outgoing edges."""
    active = np.diff(adjacency.indptr) > 0
    result = adjacency.multiply(active[None, :]).tocsr()
    result.eliminate_zeros()
    return result


def prune_triangles(adjacency: csr_array, block: int = 512) -> csr_array:
    """Drop A->C when some A->B->C has the same sign on all three edges (2016 rule).

    Decisions use the input snapshot, so the result does not depend on visit order.
    """
    matrix = adjacency.copy()
    matrix.setdiag(0)
    matrix.eliminate_zeros()
    positive = (matrix > 0).astype(np.float64).tocsr()
    negative = (matrix < 0).astype(np.float64).tocsr()
    explained = []
    for start in range(0, matrix.shape[0], block):
        rows = slice(start, start + block)
        hits = (positive[rows] @ positive).multiply(positive[rows]) + (
            negative[rows] @ negative
        ).multiply(negative[rows])
        explained.append(hits.tocsr())
    mask = vstack(explained).tocsr()
    mask.data = np.ones_like(mask.data)
    result = (matrix - matrix.multiply(mask)).tocsr()
    result.eliminate_zeros()
    return result


def absolute(adjacency: csr_array) -> csr_array:
    return abs(adjacency).tocsr()


# De-cycling: graph -> (DAG, level per node; lower level = further upstream)


def greedy_layering(adjacency: csr_array) -> tuple[csr_array, np.ndarray]:
    """2016 layering, visiting the strongest edges first (the notebook visited the weakest)."""
    edges = adjacency.tocoo()
    order = np.lexsort((edges.col, edges.row, -np.abs(edges.data)))
    layers = np.zeros(adjacency.shape[0], dtype=np.int64)
    assigned = np.zeros(len(layers), dtype=bool)
    kept = np.zeros(edges.nnz, dtype=bool)
    for position in order:
        source, target = edges.row[position], edges.col[position]
        if source == target:
            continue
        if not assigned[source] and not assigned[target]:
            layers[source], layers[target] = 0, 1
            assigned[source] = assigned[target] = True
        elif assigned[source] and not assigned[target]:
            layers[target] = layers[source] + 1
            assigned[target] = True
        elif assigned[target] and not assigned[source]:
            layers[source] = layers[target] - 1
            assigned[source] = True
        kept[position] = layers[source] < layers[target]
    dag = _rebuild(adjacency, edges.row, edges.col, np.where(kept, edges.data, 0.0))
    return dag, layers


def minimum_feedback_order(adjacency: csr_array) -> tuple[csr_array, np.ndarray]:
    """Eades-Lin-Smyth ordering on absolute weights; edges pointing backwards are dropped.

    Repeatedly take a sink (to the end), else a source (to the front), else the node with
    the largest outgoing minus incoming weight (to the front).
    """
    forward = abs(adjacency).tocsr()
    forward.setdiag(0)
    forward.eliminate_zeros()
    backward = forward.T.tocsr()
    size = forward.shape[0]
    out_count = np.diff(forward.indptr).tolist()
    in_count = np.diff(backward.indptr).tolist()
    out_weight = np.asarray(forward.sum(axis=1)).ravel().tolist()
    in_weight = np.asarray(backward.sum(axis=1)).ravel().tolist()
    removed = [False] * size

    def key(node: int) -> tuple[int, float]:
        if out_count[node] == 0:
            return (0, 0.0)
        if in_count[node] == 0:
            return (1, 0.0)
        return (2, in_weight[node] - out_weight[node])

    heap = [(*key(node), node) for node in range(size)]
    heapq.heapify(heap)
    head, tail = [], []
    while heap:
        category, priority, node = heapq.heappop(heap)
        if removed[node] or (category, priority) != key(node):
            continue
        removed[node] = True
        (tail if category == 0 else head).append(node)
        for matrix, counts, weights in (
            (forward, in_count, in_weight),
            (backward, out_count, out_weight),
        ):
            start, end = matrix.indptr[node], matrix.indptr[node + 1]
            for neighbor, weight in zip(
                matrix.indices[start:end].tolist(), matrix.data[start:end].tolist(), strict=True
            ):
                if not removed[neighbor]:
                    counts[neighbor] -= 1
                    weights[neighbor] -= weight
                    heapq.heappush(heap, (*key(neighbor), neighbor))
        if len(heap) > 8 * size:
            heap = [(*key(n), n) for n in range(size) if not removed[n]]
            heapq.heapify(heap)
    position = np.empty(size, dtype=np.int64)
    position[head + tail[::-1]] = np.arange(size)
    edges = adjacency.tocoo()
    kept = position[edges.row] < position[edges.col]
    return _rebuild(adjacency, edges.row, edges.col, np.where(kept, edges.data, 0.0)), position
