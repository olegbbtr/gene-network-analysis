from array import array
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components


@dataclass(frozen=True)
class Graph:
    nodes: tuple[str, ...]
    adjacency: csr_array


def iter_edges(path: Path):
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected source, target, weight")
            source, target, raw_weight = fields
            weight = float(raw_weight)
            if not np.isfinite(weight) or weight == 0:
                raise ValueError(f"{path}:{line_number}: weight must be finite and nonzero")
            yield source, target, weight


def intern_pair(index: dict[str, int], source: str, target: str) -> tuple[int, int]:
    for node in (source, target):
        if node not in index:
            index[node] = len(index)
    return index[source], index[target]


def read_graph(path: Path) -> Graph:
    index: dict[str, int] = {}
    rows, columns, weights = array("i"), array("i"), array("d")
    for source, target, weight in iter_edges(path):
        row, column = intern_pair(index, source, target)
        rows.append(row)
        columns.append(column)
        weights.append(weight)
    if not weights:
        raise ValueError("the edge list is empty")
    matrix = csr_array(
        (np.asarray(weights), (np.asarray(rows), np.asarray(columns))),
        shape=(len(index), len(index)),
    )
    if matrix.nnz != len(weights):
        raise ValueError("duplicate directed edges: refusing to silently sum their weights")
    return Graph(tuple(index), matrix)


def pagerank(
    adjacency: csr_array,
    *,
    alpha: float = 0.6,
    reverse: bool = False,
    weighted: bool = False,
    tolerance: float = 1e-12,
    max_iterations: int = 1000,
) -> tuple[np.ndarray, int]:
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    matrix = adjacency.T.tocsr() if reverse else adjacency.copy()
    if matrix.shape[0] == 0 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("PageRank requires a nonempty square matrix")
    matrix.data = np.abs(matrix.data) if weighted else np.ones(matrix.nnz)
    degree = np.asarray(matrix.sum(axis=1)).ravel()
    inverse = np.divide(1.0, degree, out=np.zeros_like(degree), where=degree != 0)
    transition = matrix.multiply(inverse[:, None]).T.tocsr()
    size = matrix.shape[0]
    ranks = np.full(size, 1.0 / size)
    for iteration in range(1, max_iterations + 1):
        updated = alpha * (transition @ ranks + ranks[degree == 0].sum() / size)
        updated += (1 - alpha) / size
        if np.abs(updated - ranks).sum() < tolerance:
            return updated / updated.sum(), iteration
        ranks = updated
    raise RuntimeError(f"PageRank did not converge within {max_iterations} iterations")


def component_summary(adjacency: csr_array, threshold: float = 0) -> dict:
    matrix = adjacency.copy()
    matrix.data = (np.abs(matrix.data) > threshold).astype(np.float64)
    matrix.eliminate_zeros()
    active = np.asarray(matrix.sum(axis=0) + matrix.sum(axis=1)).ravel() > 0
    count, labels = connected_components(matrix, directed=True, connection="strong")
    sizes = np.bincount(labels, minlength=count)
    weak_count = connected_components(matrix, directed=True, connection="weak", return_labels=False)
    return {
        "threshold": float(threshold),
        "edges": int(matrix.nnz),
        "active_nodes": int(active.sum()),
        "strong_components": int(count),
        "largest_strong_component": int(sizes.max()),
        "nontrivial_strong_components": int((sizes > 1).sum()),
        "weak_components": int(weak_count),
        "isolated_nodes": int((~active).sum()),
    }


def edge_weights(adjacency: csr_array, rows: np.ndarray, columns: np.ndarray) -> np.ndarray:
    """Weight of each edge rows[i] -> columns[i]; 0 where the edge is absent."""
    matrix = adjacency.tocsr()
    matrix.sort_indices()
    size = np.int64(matrix.shape[1])
    starts = np.repeat(np.arange(matrix.shape[0], dtype=np.int64), np.diff(matrix.indptr))
    keys = starts * size + matrix.indices
    wanted = np.asarray(rows, dtype=np.int64) * size + np.asarray(columns, dtype=np.int64)
    if not len(keys):
        return np.zeros(len(wanted))
    position = np.searchsorted(keys, wanted).clip(max=len(keys) - 1)
    return np.where(keys[position] == wanted, matrix.data[position], 0.0)
