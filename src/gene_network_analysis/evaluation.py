from pathlib import Path

import numpy as np
from scipy.stats import rankdata


def normalize_gene_id(value: str) -> str:
    value = value.strip()
    return value.split(".", 1)[0] if value.startswith("ENSG") else value


def read_reference(path: Path) -> set[str]:
    return {
        normalize_gene_id(line)
        for line in path.read_text().splitlines()
        if line.strip() and line.strip() != "Ensembl Gene ID"
    }


def resolve_labels(nodes: tuple[str, ...], reference: set[str]) -> tuple[np.ndarray, np.ndarray]:
    identities = [normalize_gene_id(node) for node in nodes]
    known = np.array([node.startswith("ENSG") for node in identities])
    resolved = [
        identity for identity, available in zip(identities, known, strict=True) if available
    ]
    if len(resolved) != len(set(resolved)):
        raise ValueError("multiple graph nodes resolve to the same gene")
    labels = np.array([identity in reference for identity in identities])
    return labels, known


def ranking_metrics(scores: np.ndarray, labels: np.ndarray, nodes: list[str], k: int = 100) -> dict:
    if len(scores) != len(labels) or len(nodes) != len(labels) or not len(labels):
        raise ValueError("scores, labels and node IDs must have the same nonzero length")
    if k <= 0 or not np.isfinite(scores).all():
        raise ValueError("invalid top-k or nonfinite scores")
    positives = int(labels.sum())
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        raise ValueError("evaluation requires both reference and nonreference genes")
    count = min(k, len(labels))
    order = np.lexsort((np.asarray(nodes), -scores))
    hits = int(labels[order[:count]].sum())
    expected = count * positives / len(labels)
    ranks = rankdata(scores, method="average")
    auc = (ranks[labels].sum() - positives * (positives + 1) / 2) / (positives * negatives)
    boundary = scores[order[count - 1]]
    higher = scores > boundary
    tied = scores == boundary
    slots = count - int(higher.sum())
    tie_expected = labels[higher].sum() + slots * labels[tied].mean()
    return {
        "genes": len(labels),
        "reference_genes": positives,
        "top_k": count,
        "top_k_reference_genes": hits,
        "expected_by_chance": float(expected),
        "fold_enrichment": hits / expected,
        "roc_auc": float(auc),
        "boundary_tie_size": int(tied.sum()),
        "tie_averaged_top_k_reference_genes": float(tie_expected),
        "tie_break": "ascending node ID",
    }


def random_rankings(labels: np.ndarray, k: int, runs: int = 2000, seed: int = 0) -> dict:
    """Top-k hits and AUC of uniformly random orderings, for p-values and the random row."""
    rng = np.random.default_rng(seed)
    size, positives = len(labels), int(labels.sum())
    negatives = size - positives
    count = min(k, size)
    hits = np.empty(runs, dtype=np.int64)
    aucs = np.empty(runs)
    for run in range(runs):
        ranks = rng.permutation(size)[:positives]
        hits[run] = int((ranks < count).sum())
        aucs[run] = ((size - ranks).sum() - positives * (positives + 1) / 2) / (
            positives * negatives
        )
    return {"runs": runs, "seed": seed, "top_k_reference_genes": hits, "roc_auc": aucs}


def p_value(null: np.ndarray, observed: float) -> float:
    return float((1 + (null >= observed - 1e-12).sum()) / (len(null) + 1))


def bootstrap_aucs(
    scores: dict[str, np.ndarray], labels: np.ndarray, runs: int = 1000, seed: int = 1
) -> dict[str, np.ndarray]:
    """ROC AUC of every score on the same resampled gene sets (ties count one half).

    Resamples without both reference and other genes are dropped.
    """
    from scipy.sparse import csr_array

    rng = np.random.default_rng(seed)
    size = len(labels)
    weights = np.stack(
        [np.bincount(rng.integers(0, size, size), minlength=size) for _ in range(runs)]
    )
    weights = weights.astype(np.float64)
    positive, negative = weights * labels, weights * ~labels
    result = {}
    for name, values in scores.items():
        _, group = np.unique(values, return_inverse=True)
        groups = csr_array((np.ones(size), (np.arange(size), group)))
        pos, neg = positive @ groups, negative @ groups
        below = np.cumsum(neg, axis=1) - neg
        pairs = pos.sum(axis=1) * neg.sum(axis=1)
        result[name] = (pos * (below + neg / 2)).sum(axis=1) / np.where(pairs, pairs, np.nan)
    usable = np.isfinite(next(iter(result.values()))) if result else None
    return {name: values[usable] for name, values in result.items()}
