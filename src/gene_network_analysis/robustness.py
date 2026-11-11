"""Three checks on the essential-gene result, all on held-out chromosomes.

- discovery: models trained on CEGv2 only are scored on the genes outside CEGv2 against
  DepMap's common essentials, so the test labels come from a list the model never saw;
- core machinery: the same comparison after removing the genes of the large housekeeping
  complexes (CORE_PATHWAYS, CORE_GROUP) from training and testing, to see whether the
  network only recognises the ribosome, spliceosome and proteasome modules;
- transfer: a model fitted on one dataset's network and expression scores another
  dataset's genes, each dataset's features standardised on its own training genes.

The sets and models were fixed before any of these were run. Two descriptions were added
afterwards: `near_misses` looks at the top new predictions that are not DepMap common
essentials (are they needed by some cell lines?), and `descriptor_aucs` gives the AUC of
every label-free descriptor on its own, to show which way each one points.
"""

import csv
import json
import re
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from .applications import Context, load
from .gene_sets import DEPENDENCY, hgnc_groups, read_gmt
from .graph import Graph
from .prediction import (
    DEGREE,
    SHAPE,
    NetworkFeatures,
    cross_validated,
    neighbour_share,
    network_features,
    paired_gain,
    read_gene_list,
)
from .structure import ARM, location_codes

NETWORK = ("expression", "degree", "shape", "neighbours")
MODELS = {"expression": ("expression",), "expression+network": NETWORK}
CORE_PATHWAYS = ("Translation", "Metabolism of RNA", "DNA Replication")
CORE_GROUP = re.compile("proteasome", re.IGNORECASE)
TOP = 100
LISTED_TOP = 30


def blocks(
    context: Context, features: dict[str, NetworkFeatures], members: list[str]
) -> dict[str, np.ndarray]:
    return {
        "expression": np.column_stack(
            [context.expression[f"{n}_{s}"] for n in members for s in ("mean", "sd")]
        ),
        "degree": np.column_stack([features[n].degree for n in members]),
        "shape": np.column_stack([features[n].shape for n in members]),
        "pagerank": np.column_stack([features[n].pagerank for n in members]),
    }


def held_out(
    parts: dict[str, np.ndarray],
    links: list[csr_array],
    labels: np.ndarray,
    groups: np.ndarray,
    keep: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Held-out probabilities per model; with `keep`, the other genes are removed first,
    so they are neither trained on, scored, nor counted as neighbours."""
    if keep is not None:
        parts = {name: block[keep] for name, block in parts.items()}
        links = [m[keep][:, keep] for m in links]
        labels, groups = labels[keep], groups[keep]
    return {
        name: cross_validated(parts, model, labels, groups, links) for name, model in MODELS.items()
    }


def compare(labels: np.ndarray, predictions: dict[str, np.ndarray], runs: int) -> dict:
    return {
        "listed": int(labels.sum()),
        "auc": {name: float(roc_auc_score(labels, p)) for name, p in predictions.items()},
        "gain": paired_gain(
            labels, predictions["expression+network"], predictions["expression"], runs
        ),
    }


def top_hits(scores: np.ndarray, labels: np.ndarray, top: int = TOP) -> int:
    return int(labels[np.argsort(-scores, kind="stable")[:top]].sum())


def discovery(
    known: np.ndarray, target: np.ndarray, predictions: dict[str, np.ndarray], runs: int
) -> dict:
    """Score predictions trained on `known` against `target` among the genes not `known`."""
    new = ~known
    result = compare(target[new], {name: p[new] for name, p in predictions.items()}, runs)
    result["base_rate"] = float(target[new].mean())
    result["top_hits"] = {name: top_hits(p[new], target[new]) for name, p in predictions.items()}
    return result


def read_effects(path: Path, symbols: dict[str, str], genes: list[str]) -> np.ndarray:
    """DepMap CRISPR gene effects as cell lines x genes; NaN where a gene was not screened."""
    index = {g: i for i, g in enumerate(genes)}
    with path.open(newline="") as stream:
        reader = csv.reader(stream)
        columns = [index.get(symbols.get(c.split(" (", 1)[0])) for c in next(reader)[1:]]
        rows = np.array([[float(v) if v else np.nan for v in row[1:]] for row in reader])
    values = np.full((len(rows), len(genes)), np.nan)
    kept = [j for j, column in enumerate(columns) if column is not None]
    values[:, [columns[j] for j in kept]] = rows[:, kept]
    return values


def near_misses(
    scores: np.ndarray, known: np.ndarray, target: np.ndarray, effects: np.ndarray
) -> dict:
    """Top new predictions that `target` does not list: do some cell lines need them?

    A cell line needs a gene when its effect is at most DEPENDENCY. The missed genes are
    compared with every other gene outside both lists.
    """
    screened = np.isfinite(effects).sum(axis=0)
    needed = (effects <= DEPENDENCY).sum(axis=0) / np.maximum(screened, 1)
    new = np.flatnonzero(~known)
    top = np.zeros(len(known), dtype=bool)
    top[new[np.argsort(-scores[new], kind="stable")[:TOP]]] = True
    missed = top & ~target

    def describe(chosen: np.ndarray) -> dict:
        shares = needed[chosen & (screened > 0)]
        return {
            "genes": int(chosen.sum()),
            "screened": int(len(shares)),
            "median_share_of_lines": float(np.median(shares)),
            "needed_by_a_quarter_of_lines": float((shares >= 0.25).mean()),
        }

    return {"missed": describe(missed), "others": describe(~known & ~target & ~missed)}


def descriptor_aucs(
    features: NetworkFeatures, labels: np.ndarray, keep: np.ndarray | None = None
) -> dict[str, float]:
    """AUC of each label-free descriptor on its own, among the genes in `keep` if given."""
    columns = np.column_stack([features.degree, features.shape])
    if keep is not None:
        columns, labels = columns[keep], labels[keep]
    return {
        name: float(roc_auc_score(labels, columns[:, j]))
        for j, name in enumerate((*DEGREE, *SHAPE))
    }


def core_machinery(context: Context, reactome: dict[str, set[str]]) -> np.ndarray:
    core = set().union(*(reactome[name] for name in CORE_PATHWAYS))
    for name, members in hgnc_groups(context.table).items():
        if CORE_GROUP.search(name):
            core |= members
    return np.array([g in core for g in context.genes])


def transfer(
    parts: dict[str, dict[str, np.ndarray]],
    links: dict[str, csr_array],
    labels: np.ndarray,
    groups: np.ndarray,
    folds: int = 5,
) -> dict[str, dict[str, float]]:
    """AUC of the NETWORK model fitted on dataset `source` and scored on `target`.

    Every dataset is standardised on its own training genes, and its neighbour share uses
    the training-fold labels on its own links, so only the fitted coefficients move.
    """
    names = list(parts)
    predictions = {s: {t: np.zeros(len(labels)) for t in names} for s in names}
    for train, test in GroupKFold(folds).split(labels, labels, groups):
        known = np.zeros(len(labels), dtype=bool)
        known[train] = labels[train]
        scaled = {}
        for name in names:
            columns = np.column_stack(
                [*(parts[name][p] for p in NETWORK[:-1]), neighbour_share([links[name]], known)]
            )
            scaled[name] = StandardScaler().fit(columns[train]).transform(columns)
        for source in names:
            model = LogisticRegression(max_iter=5000).fit(scaled[source][train], labels[train])
            for target in names:
                predictions[source][target][test] = model.predict_proba(scaled[target][test])[:, 1]
    return {s: {t: float(roc_auc_score(labels, predictions[s][t])) for t in names} for s in names}


def transfer_summary(
    matrix: dict[str, dict[str, float]], own_expression: dict[str, float], state: dict[str, str]
) -> dict:
    """`state` names each network's kind (tumour, healthy, ...): losses are also given
    per pair of kinds, source first."""
    names = list(matrix)
    pairs = [(s, t) for s in names for t in names if s != t]
    loss = [matrix[t][t] - matrix[s][t] for s, t in pairs]
    kinds = sorted(set(state.values()))
    return {
        "median_loss_by_state": {
            f"{a} to {b}": float(
                np.median(
                    [
                        x
                        for x, (s, t) in zip(loss, pairs, strict=True)
                        if (state[s], state[t]) == (a, b)
                    ]
                )
            )
            for a in kinds
            for b in kinds
            if any((state[s], state[t]) == (a, b) for s, t in pairs)
        },
        "pairs": len(pairs),
        "own_auc_range": [min(matrix[n][n] for n in names), max(matrix[n][n] for n in names)],
        "transferred_auc_range": [
            min(matrix[s][t] for s, t in pairs),
            max(matrix[s][t] for s, t in pairs),
        ],
        "loss_median": float(np.median(loss)),
        "loss_max": float(max(loss)),
        "pairs_above_own_expression": sum(matrix[s][t] > own_expression[t] for s, t in pairs),
    }


def run(root: Path, gene_table: Path, references: Path, output: Path, runs: int = 500) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    arms = location_codes(tuple(genes), context.table, ARM)
    names = context.plain
    features = {
        n: network_features(Graph(tuple(genes), context.adjacency[n]), genes, arms) for n in names
    }
    ceg = read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols)
    depmap = read_gene_list(
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", context.symbols
    )
    known = np.array([g in ceg for g in genes])
    target = np.array([g in depmap for g in genes])
    reactome = read_gmt(references / "reactome/ReactomePathways.gmt.zip", context.symbols)
    core = core_machinery(context, reactome)
    effects = read_effects(references / "depmap-24q2/CRISPRGeneEffect.csv", context.symbols, genes)
    report: dict = {
        "genes": len(genes),
        "cegv2": int(known.sum()),
        "depmap_outside_cegv2": int((target & ~known).sum()),
        "core_machinery": {
            "genes": int(core.sum()),
            "cegv2": int((core & known).sum()),
            "depmap": int((core & target).sum()),
        },
        "descriptor_auc": {n: descriptor_aucs(features[n], known) for n in names},
        "discovery": {},
        "without_core_machinery": {"cegv2": {}, "depmap": {}},
    }
    sets = {name: [name] for name in names} | {"all": names}
    for set_name, members in sets.items():
        parts = blocks(context, features, members)
        links = [features[n].links for n in members]
        predictions = held_out(parts, links, known, context.groups)
        report["discovery"][set_name] = discovery(known, target, predictions, runs)
        if set_name == "all":
            new = np.flatnonzero(~known)
            order = new[np.argsort(-predictions["expression+network"][new], kind="stable")]
            report["discovery_top_genes"] = [
                {"symbol": context.table[genes[i]]["symbol"], "depmap": bool(target[i])}
                for i in order[:LISTED_TOP]
            ]
            report["discovery_misses"] = near_misses(
                predictions["expression+network"], known, target, effects
            )
        for list_name, labels in (("cegv2", known), ("depmap", target)):
            kept = held_out(parts, links, labels, context.groups, ~core)
            report["without_core_machinery"][list_name][set_name] = compare(
                labels[~core], kept, runs
            )
        print(
            set_name,
            json.dumps(
                {
                    "discovery": report["discovery"][set_name],
                    "without_core": {
                        k: v[set_name] for k, v in report["without_core_machinery"].items()
                    },
                }
            ),
            flush=True,
        )
    matrix = transfer(
        {n: blocks(context, features, [n]) for n in names},
        {n: features[n].links for n in names},
        known,
        context.groups,
    )
    own_expression = {
        n: float(
            roc_auc_score(
                known,
                cross_validated(
                    blocks(context, features, [n]), ("expression",), known, context.groups, []
                ),
            )
        )
        for n in names
    }
    report["transfer"] = {
        "summary": transfer_summary(
            matrix, own_expression, {n: context.stats[n]["state"] for n in names}
        ),
        "own_expression_auc": own_expression,
        "auc": matrix,
    }
    print("transfer", json.dumps(report["transfer"]["summary"]), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    (output / "robustness.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
