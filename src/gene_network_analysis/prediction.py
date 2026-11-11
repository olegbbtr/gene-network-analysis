"""Do per-gene network descriptors predict gene properties beyond expression level?

For each gene list (for example core essential genes), a logistic regression is trained on
nested feature sets and scored on genes from held-out chromosomes, so co-located genes never
sit on both sides of a split:

- expression: the gene's mean and SD in the dataset(s);
- + degree: log in- and out-degree, the simplest network feature;
- + shape: the remaining label-free descriptors (signs, weight, net flow, clustering,
  share of edges on the gene's own chromosome arm);
- + neighbours: share of a gene's neighbours that are listed, using training labels only;
- + pagerank: plain and reverse PageRank on the directed graph (label-free);
- + personalised: PageRank restarting from training-fold listed genes, next to neighbours.

Gains between models are paired bootstraps over genes of the held-out predictions.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from .diagnostics import chromosome
from .gene_sets import neighbour_scores, personalised_pagerank
from .graph import Graph, pagerank, read_graph
from .public import read_expression
from .references import read_gene_table
from .structure import ARM, location_codes

DEGREE = ("log_out", "log_in")
SHAPE = ("log_positive_out", "log_negative_out", "log_weight_out", "net_flow", "clustering",
         "arm_share")  # fmt: skip
MODELS = {
    "expression": ("expression",),
    "degree": ("degree",),
    "shape": ("degree", "shape"),
    "neighbours": ("neighbours",),
    "expression+degree": ("expression", "degree"),
    "expression+shape": ("expression", "degree", "shape"),
    "expression+shape+neighbours": ("expression", "degree", "shape", "neighbours"),
    "expression+shape+pagerank": ("expression", "degree", "shape", "pagerank"),
    "expression+shape+neighbours+personalised": (
        "expression",
        "degree",
        "shape",
        "neighbours",
        "personalised",
    ),
}
GAINS = {
    "network_beyond_expression": ("expression+shape", "expression"),
    "shape_beyond_degree": ("expression+shape", "expression+degree"),
    "neighbours_beyond_shape": ("expression+shape+neighbours", "expression+shape"),
    "pagerank_beyond_shape": ("expression+shape+pagerank", "expression+shape"),
    "personalised_beyond_neighbours": (
        "expression+shape+neighbours+personalised",
        "expression+shape+neighbours",
    ),
}


def aligned(graph: Graph, genes: list[str]) -> csr_array:
    """The graph's adjacency reindexed to `genes`; genes absent from the graph get no edges."""
    index = {g: i for i, g in enumerate(graph.nodes)}
    edges = graph.adjacency.tocoo()
    position = np.array([index.get(g, -1) for g in genes])
    new = np.full(len(graph.nodes), -1)
    new[position[position >= 0]] = np.flatnonzero(position >= 0)
    keep = (new[edges.row] >= 0) & (new[edges.col] >= 0)
    return csr_array(
        (edges.data[keep], (new[edges.row[keep]], new[edges.col[keep]])),
        shape=(len(genes), len(genes)),
    )


def undirected(adjacency: csr_array) -> csr_array:
    present = (adjacency != 0).astype(np.float64)
    links = ((present + present.T) > 0).astype(np.float64)
    links.setdiag(0)
    links.eliminate_zeros()
    return links.tocsr()


def positive_links(adjacency: csr_array) -> csr_array:
    """Undirected 0/1 links between genes whose weights, summed over both directions, are
    positive."""
    links = ((adjacency + adjacency.T) > 0).astype(np.float64)
    links.setdiag(0)
    links.eliminate_zeros()
    return links.tocsr()


def descriptors(adjacency: csr_array, arms: np.ndarray) -> dict[str, np.ndarray]:
    size = adjacency.shape[0]
    present = (adjacency != 0).astype(np.float64)
    weight = abs(adjacency)
    out = np.asarray(present.sum(axis=1)).ravel()
    links = undirected(adjacency)
    degree = np.asarray(links.sum(axis=1)).ravel()
    triangles = np.asarray((links @ links).multiply(links).sum(axis=1)).ravel() / 2
    pairs = degree * (degree - 1) / 2
    edges = present.tocoo()
    same = (arms[edges.row] == arms[edges.col]) & (arms[edges.row] >= 0)
    local = np.bincount(edges.row, weights=same, minlength=size)
    return {
        "log_out": np.log1p(out),
        "log_in": np.log1p(np.asarray(present.sum(axis=0)).ravel()),
        "log_positive_out": np.log1p(np.asarray((adjacency > 0).sum(axis=1)).ravel()),
        "log_negative_out": np.log1p(np.asarray((adjacency < 0).sum(axis=1)).ravel()),
        "log_weight_out": np.log1p(np.asarray(weight.sum(axis=1)).ravel()),
        "net_flow": np.asarray(weight.sum(axis=1)).ravel() - np.asarray(weight.sum(axis=0)).ravel(),
        "clustering": np.divide(triangles, pairs, out=np.zeros(size), where=pairs > 0),
        "arm_share": np.divide(local, out, out=np.zeros(size), where=out > 0),
    }


def neighbour_share(links: list[csr_array], labels: np.ndarray) -> np.ndarray:
    """Per graph, the share of each gene's neighbours carrying a label (genes x graphs)."""
    return np.column_stack([neighbour_scores(m, labels[:, None])[:, 0] for m in links])


def personalised(links: list[csr_array], labels: np.ndarray) -> np.ndarray:
    """Per graph, log PageRank restarting from labelled genes, scaled so the mean is 1."""
    return np.column_stack(
        [np.log1p(len(labels) * personalised_pagerank(m, labels[:, None])[:, 0]) for m in links]
    )


LABELLED = {"neighbours": neighbour_share, "personalised": personalised}
# A walk restarting at a known gene scores that gene itself, so training genes get this
# feature from inner folds that hide them, as test genes are hidden.
CROSS_FITTED = {"personalised"}


def labelled_feature(
    part: str,
    function,
    links: list[csr_array],
    labels: np.ndarray,
    groups: np.ndarray,
    train: np.ndarray,
    folds: int,
) -> np.ndarray:
    known = np.zeros(len(labels), dtype=bool)
    known[train] = labels[train]
    column = function(links, known)
    if part in CROSS_FITTED:
        for inner_train, inner_test in GroupKFold(folds - 1).split(train, groups=groups[train]):
            inner = np.zeros(len(labels), dtype=bool)
            inner[train[inner_train]] = labels[train[inner_train]]
            column[train[inner_test]] = function(links, inner)[train[inner_test]]
    return column


def cross_validated(
    blocks: dict[str, np.ndarray],
    parts: tuple[str, ...],
    labels: np.ndarray,
    groups: np.ndarray,
    links: list[csr_array],
    folds: int = 5,
    labelled: dict | None = None,
) -> np.ndarray:
    """Held-out probability per gene; labelled parts use training-fold labels only.

    `labelled` maps a part to f(links, known labels) -> genes x columns; default LABELLED.
    """
    labelled = LABELLED if labelled is None else labelled
    predictions = np.zeros(len(labels))
    for train, test in GroupKFold(folds).split(labels, labels, groups):
        columns = [
            labelled_feature(p, labelled[p], links, labels, groups, train, folds)
            if p in labelled
            else blocks[p]
            for p in parts
        ]
        features = np.column_stack(columns)
        scaler = StandardScaler().fit(features[train])
        model = LogisticRegression(max_iter=5000).fit(
            scaler.transform(features[train]), labels[train]
        )
        predictions[test] = model.predict_proba(scaler.transform(features[test]))[:, 1]
    return predictions


def paired_gain(
    labels: np.ndarray, better: np.ndarray, worse: np.ndarray, runs: int = 1000, seed: int = 0
) -> dict:
    rng = np.random.default_rng(seed)
    gains = []
    for _ in range(runs):
        sample = rng.integers(0, len(labels), len(labels))
        if labels[sample].all() or not labels[sample].any():
            continue
        y = labels[sample]
        gains.append(roc_auc_score(y, better[sample]) - roc_auc_score(y, worse[sample]))
    observed = roc_auc_score(labels, better) - roc_auc_score(labels, worse)
    low, high = np.quantile(gains, [0.025, 0.975])
    return {"gain": float(observed), "ci95": [float(low), float(high)]}


@dataclass(frozen=True)
class NetworkFeatures:
    degree: np.ndarray
    shape: np.ndarray
    pagerank: np.ndarray
    links: csr_array


def pagerank_features(adjacency: csr_array) -> np.ndarray:
    """Log plain and reverse PageRank of the unweighted directed graph, mean 1 before log."""
    size = adjacency.shape[0]
    return np.column_stack(
        [
            np.log(size * pagerank(adjacency, alpha=0.85, reverse=r, tolerance=1e-10)[0])
            for r in (False, True)
        ]
    )


def network_features(graph: Graph, genes: list[str], arms: np.ndarray) -> NetworkFeatures:
    adjacency = aligned(graph, genes)
    features = descriptors(adjacency, arms)
    return NetworkFeatures(
        np.column_stack([features[n] for n in DEGREE]),
        np.column_stack([features[n] for n in SHAPE]),
        pagerank_features(adjacency),
        undirected(adjacency),
    )


def evaluate(
    expression: np.ndarray,
    networks: list[NetworkFeatures],
    labels: np.ndarray,
    groups: np.ndarray,
    runs: int = 1000,
) -> dict:
    """AUC of every model in MODELS and the gains in GAINS for one set of networks."""
    blocks = {
        "expression": expression,
        "degree": np.column_stack([n.degree for n in networks]),
        "shape": np.column_stack([n.shape for n in networks]),
        "pagerank": np.column_stack([n.pagerank for n in networks]),
    }
    links = [n.links for n in networks]
    predictions = {
        name: cross_validated(blocks, parts, labels, groups, links)
        for name, parts in MODELS.items()
    }
    return {
        "auc": {name: float(roc_auc_score(labels, p)) for name, p in predictions.items()},
        "gains": {
            name: paired_gain(labels, predictions[a], predictions[b], runs)
            for name, (a, b) in GAINS.items()
        },
    }


def read_gene_list(path: Path, symbols: dict[str, str]) -> set[str]:
    """Ensembl IDs of the symbols in a list's first column; the header line is skipped.

    DepMap's `SYMBOL (Entrez ID)` entries are read as SYMBOL.
    """
    lines = path.read_text().splitlines()[1:]
    names = [line.split("\t")[0].split(" (", 1)[0] for line in lines]
    return {symbols[name] for name in names if name in symbols}


def predict(root: Path, gene_table: Path, lists: dict[str, Path], runs: int = 1000) -> dict:
    """Every gene list against each network alone, all tumour, all healthy and all networks.

    `root` is the output of `public.build`: one folder per dataset and `expression.tsv`.
    """
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    genes, expression = read_expression(root / "expression.tsv")
    arms = location_codes(tuple(genes), table, ARM)
    groups = np.array([chromosome(table.get(g, {}).get("location", "")) for g in genes])
    stats = json.loads((root / "stats.json").read_text())
    names = [name for name in stats if name != "settings"]
    features = {
        name: network_features(read_graph(root / name / "graph.txt"), genes, arms) for name in names
    }

    def columns(dataset: str) -> list[np.ndarray]:
        base = dataset.removesuffix("-copy-number-removed")
        return [expression[f"{base}_mean"], expression[f"{base}_sd"]]

    plain = [n for n in names if not n.endswith("-copy-number-removed")]
    sets = {name: [name] for name in names} | {
        "all-tumour": [n for n in plain if stats[n]["state"] == "tumour"],
        "all-healthy": [n for n in plain if stats[n]["state"] == "healthy"],
        "all": plain,
    }
    report: dict = {"genes": len(genes), "networks": {n: stats[n] for n in names}, "lists": {}}
    for list_name, path in lists.items():
        listed = read_gene_list(path, symbols)
        labels = np.array([g in listed for g in genes])
        report["lists"][list_name] = {
            "listed": int(labels.sum()),
            "results": {
                set_name: evaluate(
                    np.column_stack([c for n in members for c in columns(n)]),
                    [features[n] for n in members],
                    labels,
                    groups,
                    runs,
                )
                for set_name, members in sets.items()
            },
        }
    return report
