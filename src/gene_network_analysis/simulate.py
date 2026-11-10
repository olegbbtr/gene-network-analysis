"""A regulatory network with known regulators and edges, for testing the pipelines.

The network is a layered DAG: regulators in layer 0 act on regulators of later layers
and on non-regulator genes. Expression follows a linear model with Gaussian noise.
The inferred graph z-scores each correlation (or ridge partial correlation) against the
target's values for all candidate regulators, so A->B and B->A get different weights, and
keeps |z| > threshold.
"""

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from .graph import write_graph


@dataclass(frozen=True)
class Settings:
    genes: int = 2000
    regulators: int = 200
    layers: int = 4
    mean_targets: float = 20.0
    regulator_target_share: float = 0.2
    repression: float = 0.3
    samples: int = 500
    noise: float = 1.0
    extra_senders: float = 0.3
    threshold: float = 2.35
    inference: str = "correlation"  # or "partial": ridge-regularised partial correlation
    ridge: float = 0.5
    seed: int = 0


def true_network(settings: Settings, rng: np.random.Generator) -> tuple[csr_array, np.ndarray]:
    """Signed DAG over genes; regulators are genes 0..regulators-1, ordered by layer."""
    layer = np.sort(rng.integers(0, settings.layers, settings.regulators))
    others = np.arange(settings.regulators, settings.genes)
    rows, columns, weights = [], [], []
    for regulator in range(settings.regulators):
        count = max(1, int(rng.lognormal(np.log(settings.mean_targets) - 0.5, 1.0)))
        lower = np.flatnonzero(layer > layer[regulator])
        from_regulators = (
            rng.random(count) < settings.regulator_target_share if len(lower) else np.zeros(count)
        )
        targets = np.where(
            from_regulators.astype(bool),
            rng.choice(lower, count) if len(lower) else 0,
            rng.choice(others, count),
        )
        targets = np.unique(targets)
        signs = np.where(rng.random(len(targets)) < settings.repression, -1.0, 1.0)
        rows += [regulator] * len(targets)
        columns += targets.tolist()
        weights += (signs * rng.uniform(0.5, 1.0, len(targets))).tolist()
    network = csr_array(
        (np.asarray(weights), (np.asarray(rows), np.asarray(columns))),
        shape=(settings.genes, settings.genes),
    )
    return network, layer


def expression(network: csr_array, settings: Settings, rng: np.random.Generator) -> np.ndarray:
    """Samples x genes; each gene is the weighted sum of its regulators plus noise."""
    dense = network.toarray()
    values = rng.normal(0, settings.noise, (settings.samples, settings.genes))
    for gene in range(settings.regulators):
        # Regulators are ordered by layer and only act on later genes, so each one is
        # final before its effect is added.
        values += np.outer(values[:, gene], dense[gene])
    return (values - values.mean(axis=0)) / values.std(axis=0)


def association(data: np.ndarray, method: str, ridge: float) -> np.ndarray:
    correlation = (data.T @ data) / len(data)
    if method == "partial":
        precision = np.linalg.inv(correlation + ridge * np.eye(len(correlation)))
        scale = np.sqrt(np.diag(precision))
        correlation = -precision / np.outer(scale, scale)
    elif method != "correlation":
        raise ValueError(f"unknown inference method: {method}")
    np.fill_diagonal(correlation, 0)
    return correlation


def infer(
    data: np.ndarray,
    senders: np.ndarray,
    threshold: float,
    method: str = "correlation",
    ridge: float = 0.5,
) -> csr_array:
    correlation = association(data, method, ridge)
    fisher = np.arctanh(np.clip(correlation[senders], -0.999999, 0.999999))
    z = (fisher - fisher.mean(axis=0)) / fisher.std(axis=0)
    z[np.arange(len(senders)), senders] = 0
    keep = np.abs(z) > threshold
    rows, columns = np.nonzero(keep)
    return csr_array(
        (z[rows, columns], (senders[rows], columns)), shape=(data.shape[1], data.shape[1])
    )


def simulate(settings: Settings, output: Path) -> dict:
    rng = np.random.default_rng(settings.seed)
    network, layer = true_network(settings, rng)
    data = expression(network, settings, rng)
    regulators = np.arange(settings.regulators)
    others = np.arange(settings.regulators, settings.genes)
    extra = rng.choice(others, int(settings.extra_senders * len(others)), replace=False)
    senders = np.sort(np.concatenate([regulators, extra]))
    graph = infer(data, senders, settings.threshold, settings.inference, settings.ridge)
    names = [f"ENSGSIM{index:06d}" for index in range(settings.genes)]
    output.mkdir(parents=True, exist_ok=True)
    write_graph(output / "graph.txt", graph, names)
    (output / "regulators.txt").write_text("".join(f"{names[i]}\n" for i in regulators))
    (output / "top-regulators.txt").write_text(
        "".join(f"{names[i]}\n" for i in regulators[layer == 0])
    )
    truth = network.tocoo()
    with (output / "true-edges.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["source_gene", "target_gene", "sign"])
        for u, v, w in zip(truth.row, truth.col, truth.data, strict=True):
            writer.writerow([names[u], names[v], int(np.sign(w))])
    return {
        "true_edges": int(network.nnz),
        "regulators": settings.regulators,
        "top_layer_regulators": int((layer == 0).sum()),
        "senders": int(len(senders)),
        "inferred_edges": int(graph.nnz),
        "inferred_negative": int((graph.data < 0).sum()),
    }
