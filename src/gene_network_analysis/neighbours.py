"""Scores from the TF labels of other genes: a gene's own label never enters its score."""

import numpy as np
from scipy.sparse import csr_array

from .steps import top_outgoing


def _without_diagonal(adjacency: csr_array) -> csr_array:
    weight = abs(adjacency).tocoo()
    keep = weight.row != weight.col
    return csr_array(
        (weight.data[keep], (weight.row[keep], weight.col[keep])), shape=adjacency.shape
    )


def neighbour_vote(adjacency: csr_array, labels: np.ndarray, direction: str = "both") -> np.ndarray:
    """Share of each gene's absolute edge weight that connects it to labelled genes."""
    weight = _without_diagonal(adjacency)
    matrix = {"out": weight, "in": weight.T, "both": weight + weight.T}[direction].tocsr()
    total = np.asarray(matrix.sum(axis=1)).ravel()
    labelled = matrix @ labels.astype(float)
    return np.divide(labelled, total, out=np.zeros_like(total), where=total > 0)


def target_profiles(adjacency: csr_array, rows: np.ndarray, targets: int, signed: bool):
    """Unit-length rows over each gene's strongest outgoing edges (weights set to ±1 or 1)."""
    profiles = top_outgoing(adjacency, k=targets)[rows].tocsr() if targets else adjacency[rows]
    profiles = csr_array(profiles)
    profiles.data = np.sign(profiles.data) if signed else np.ones_like(profiles.data)
    norm = np.sqrt(np.diff(profiles.indptr)).astype(float)
    return csr_array(profiles.multiply(1 / np.where(norm > 0, norm, 1)[:, None]))


def similarity_vote(
    adjacency: csr_array,
    rows: np.ndarray,
    labels: np.ndarray,
    *,
    targets: int = 50,
    neighbours: int = 20,
    signed: bool = False,
    block: int = 1024,
) -> np.ndarray:
    """Similarity-weighted share of TFs among each gene's most similar target sets.

    labels are for rows; a gene is never its own neighbour.
    """
    profiles = target_profiles(adjacency, rows, targets, signed)
    transposed = profiles.T.tocsr()
    labels = labels.astype(float)
    count = min(neighbours, len(rows) - 1)
    scores = np.zeros(len(rows))
    for start in range(0, len(rows), block):
        stop = min(start + block, len(rows))
        similarity = (profiles[start:stop] @ transposed).toarray()
        similarity[np.arange(stop - start), np.arange(start, stop)] = -np.inf
        nearest = np.argpartition(-similarity, count - 1, axis=1)[:, :count]
        weight = np.clip(np.take_along_axis(similarity, nearest, axis=1), 0, None)
        total = weight.sum(axis=1)
        voted = (weight * labels[nearest]).sum(axis=1)
        scores[start:stop] = np.divide(voted, total, out=np.zeros_like(total), where=total > 0)
    return scores
