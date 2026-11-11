"""Graph-only descriptors and edge-weight summaries for the poster comparison."""

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.linalg import eigsh
from scipy.stats import rankdata


def rank_positions(scores: np.ndarray) -> np.ndarray:
    """Top-first rank positions per column, in [0, 1], averaging tied positions."""
    return (rankdata(-scores, axis=0, method="average") - 1) / max(len(scores) - 1, 1)


def core_depth(links: csr_array) -> np.ndarray:
    """Core numbers of an undirected graph, using degree bins in O(nodes + edges)."""
    degree = np.diff(links.indptr).copy()
    if not len(degree):
        return degree
    order = np.argsort(degree, kind="stable")
    position = np.empty_like(order)
    position[order] = np.arange(len(order))
    counts = np.bincount(degree)
    starts = np.cumsum(counts) - counts
    for i in range(len(order)):
        node = order[i]
        for neighbour in links.indices[links.indptr[node] : links.indptr[node + 1]]:
            if degree[neighbour] <= degree[node]:
                continue
            level = degree[neighbour]
            old = position[neighbour]
            first = starts[level]
            other = order[first]
            order[old], order[first] = other, neighbour
            position[other], position[neighbour] = old, first
            starts[level] += 1
            degree[neighbour] -= 1
    return degree


def eigenvector_score(links: csr_array) -> np.ndarray:
    size = links.shape[0]
    if not size:
        return np.empty(0)
    if size == 1 or not links.nnz:
        return np.full(size, 1 / size)
    _, vectors = eigsh(
        links.astype(float), k=1, which="LA", v0=np.ones(size), tol=1e-8, maxiter=1000
    )
    scores = vectors[:, 0]
    if scores.sum() < 0:
        scores = -scores
    scores = np.maximum(scores, 0)
    return scores / scores.sum()


def weight_summaries(adjacency: csr_array) -> np.ndarray:
    """Per gene, outgoing then incoming nonzero-weight mean, SD, absolute mean and
    positive share. Missing edges do not enter the distributions; counts are not features.
    """
    result = []
    for matrix in (adjacency.tocsr(), adjacency.T.tocsr()):
        count = np.diff(matrix.indptr)
        denominator = np.maximum(count, 1)
        mean = np.asarray(matrix.sum(axis=1)).ravel() / denominator
        square = np.asarray(matrix.multiply(matrix).sum(axis=1)).ravel() / denominator
        result.extend(
            [
                mean,
                np.sqrt(np.maximum(square - mean**2, 0)),
                np.asarray(abs(matrix).sum(axis=1)).ravel() / denominator,
                np.asarray((matrix > 0).sum(axis=1)).ravel() / denominator,
            ]
        )
    return np.column_stack(result)
