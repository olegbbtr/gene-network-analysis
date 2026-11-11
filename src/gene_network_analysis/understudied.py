"""Does the network help most for the genes that STRING knows least about?

STRING draws on the literature, curated databases and experiments, so genes that have
been studied little have few links in it. Genes are split into three groups at the
terciles of their number of STRING links (combined score >= 700, among the 12,000 genes).
Two models of essential genes are fitted on all genes with held-out chromosomes, and their
predictions are compared inside each group:

- expression + STRING: the share of a gene's STRING partners that are listed in the
  training fold, and the log of its number of STRING links;
- the same plus the co-expression network (degree, shape, neighbour share).

The groups and models were fixed before the first run.
"""

import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from .applications import load
from .baselines import read_string
from .graph import Graph
from .prediction import (
    cross_validated,
    neighbour_share,
    network_features,
    paired_gain,
    read_gene_list,
)
from .references import GeneMap
from .robustness import blocks
from .structure import ARM, location_codes

STRING = ("expression", "string_links", "string")
MODELS = {"string": STRING, "string+network": (*STRING, "degree", "shape", "neighbours")}
GROUPS = ("fewest", "middle", "most")


def terciles(values: np.ndarray) -> np.ndarray:
    """Group 0, 1 or 2 per gene, cut at the terciles; equal values share a group."""
    low, high = np.quantile(values, [1 / 3, 2 / 3])
    return (values > low).astype(int) + (values > high)


def by_group(
    labels: np.ndarray, predictions: dict[str, np.ndarray], group: np.ndarray, runs: int
) -> dict:
    result = {}
    for number, name in enumerate(GROUPS):
        chosen = group == number
        result[name] = {
            "genes": int(chosen.sum()),
            "listed": int(labels[chosen].sum()),
            "auc": {
                model: float(roc_auc_score(labels[chosen], p[chosen]))
                for model, p in predictions.items()
            },
            "gain": paired_gain(
                labels[chosen],
                predictions["string+network"][chosen],
                predictions["string"][chosen],
                runs,
            ),
        }
    return result


def run(root: Path, gene_table: Path, references: Path, output: Path, runs: int = 500) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    arms = location_codes(tuple(genes), context.table, ARM)
    names = context.plain
    features = {
        n: network_features(Graph(tuple(genes), context.adjacency[n]), genes, arms) for n in names
    }
    string = read_string(
        references / "string-v12/9606.protein.info.v12.0.txt.gz",
        references / "string-v12/9606.protein.links.v12.0.txt.gz",
        genes,
        GeneMap(references / "hgnc/hgnc_complete_set.txt").by_symbol,
    )
    partners = np.asarray(string.sum(axis=1)).ravel()
    group = terciles(partners)
    labelled = {
        "neighbours": neighbour_share,
        "string": lambda _, known: neighbour_share([string], known),
    }
    lists = {
        "cegv2": references / "gene-lists/CEGv2.txt",
        "depmap": references / "depmap-24q2/CRISPRInferredCommonEssentials.csv",
    }
    report: dict = {
        "string_links_per_group": {
            name: [int(partners[group == i].min()), int(partners[group == i].max())]
            for i, name in enumerate(GROUPS)
        },
        "lists": {},
    }
    for list_name, path in lists.items():
        listed = read_gene_list(path, context.symbols)
        labels = np.array([g in listed for g in genes])
        report["lists"][list_name] = {}
        for set_name, members in ({name: [name] for name in names} | {"all": names}).items():
            parts = blocks(context, features, members)
            parts["string_links"] = np.log1p(partners)[:, None]
            links = [features[n].links for n in members]
            predictions = {
                name: cross_validated(
                    parts, model, labels, context.groups, links, labelled=labelled
                )
                for name, model in MODELS.items()
            }
            result = by_group(labels, predictions, group, runs)
            report["lists"][list_name][set_name] = result
            print(list_name, set_name, json.dumps(result), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    (output / "understudied.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
