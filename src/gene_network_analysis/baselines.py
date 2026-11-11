"""Does a co-expression network predict listed genes better than simpler or existing sources?

- correlation: a gene's mean correlation across samples with training-fold listed genes,
  with no cutoff and no network;
- string: share of a gene's STRING neighbours (combined score >= STRING_SCORE) that are
  listed in the training fold. STRING mixes experiments, curated databases, text mining
  and its own co-expression, so it is a strong, partly circular reference.

Per dataset, held-out-chromosome logistic models compare expression plus a baseline with
expression plus the network (degree, shape, neighbour share) and with both.
"""

import gzip
import json
from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from sklearn.metrics import roc_auc_score

from .diagnostics import chromosome
from .graph import read_graph
from .prediction import (
    cross_validated,
    neighbour_share,
    network_features,
    paired_gain,
    read_gene_list,
)
from .public import DATASETS, read_expression, read_samples, rows, standardise
from .references import GeneMap, read_gene_table
from .structure import ARM, location_codes

STRING_SCORE = 700
NETWORK = ("degree", "shape", "neighbours")
MODELS = {
    "expression": ("expression",),
    "expression+correlation": ("expression", "correlation"),
    "expression+string": ("expression", "string"),
    "expression+network": ("expression", *NETWORK),
    "expression+correlation+network": ("expression", "correlation", *NETWORK),
    "expression+string+network": ("expression", "string", *NETWORK),
}
GAINS = {
    "network_beyond_correlation": ("expression+correlation+network", "expression+correlation"),
    "correlation_beyond_network": ("expression+correlation+network", "expression+network"),
    "network_beyond_string": ("expression+string+network", "expression+string"),
    "string_beyond_network": ("expression+string+network", "expression+network"),
}


def mean_correlation(samples: np.ndarray, known: np.ndarray) -> np.ndarray:
    """Each gene's mean correlation with the known genes other than itself.

    `samples` is samples x genes, standardised per gene.
    """
    total = samples.T @ samples[:, known].sum(axis=1) / len(samples)
    others = known.sum() - known
    return np.divide(total - known, others, out=np.zeros(len(total)), where=others > 0)


def read_string(
    info: Path, links: Path, genes: list[str], lookup: Callable[[str], str | None]
) -> csr_array:
    """Undirected STRING links at or above STRING_SCORE between `genes`.

    `lookup` maps a STRING preferred name to a gene ID, or None.
    """
    index = {g: i for i, g in enumerate(genes)}
    with gzip.open(info, "rt") as stream:
        next(stream)
        protein = {}
        for line in stream:
            fields = line.split("\t")
            gene = lookup(fields[1])
            if gene in index:
                protein[fields[0]] = index[gene]
    rows_, columns = [], []
    with gzip.open(links, "rt") as stream:
        next(stream)
        for line in stream:
            first, second, score = line.split()
            if int(score) >= STRING_SCORE and first in protein and second in protein:
                rows_.append(protein[first])
                columns.append(protein[second])
    matrix = csr_array((np.ones(len(rows_)), (rows_, columns)), shape=(len(genes), len(genes)))
    matrix = ((matrix + matrix.T) > 0).astype(np.float64)
    matrix.setdiag(0)
    matrix.eliminate_zeros()
    return matrix.tocsr()


def compare(
    blocks: dict[str, np.ndarray],
    samples: list[np.ndarray],
    string: csr_array,
    links: list[csr_array],
    labels: np.ndarray,
    groups: np.ndarray,
    runs: int,
) -> dict:
    labelled = {
        "neighbours": neighbour_share,
        "correlation": lambda _, known: np.column_stack(
            [mean_correlation(s, known) for s in samples]
        ),
        "string": lambda _, known: neighbour_share([string], known),
    }
    predictions = {
        name: cross_validated(blocks, parts, labels, groups, links, labelled=labelled)
        for name, parts in MODELS.items()
    }
    return {
        "auc": {name: float(roc_auc_score(labels, p)) for name, p in predictions.items()},
        "gains": {
            name: paired_gain(labels, predictions[a], predictions[b], runs)
            for name, (a, b) in GAINS.items()
        },
    }


def run(
    root: Path, gene_table: Path, references: Path, lists: dict[str, Path], output: Path
) -> dict:
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    genes, expression = read_expression(root / "expression.tsv")
    arms = location_codes(tuple(genes), table, ARM)
    groups = np.array([chromosome(table.get(g, {}).get("location", "")) for g in genes])
    string = read_string(
        references / "string-v12/9606.protein.info.v12.0.txt.gz",
        references / "string-v12/9606.protein.links.v12.0.txt.gz",
        genes,
        GeneMap(references / "hgnc/hgnc_complete_set.txt").by_symbol,
    )
    samples, features = {}, {}
    for dataset in DATASETS:
        matrix, _ = read_samples(dataset, references, symbols)
        samples[dataset.name] = standardise(rows(matrix, genes))
        features[dataset.name] = network_features(
            read_graph(root / dataset.name / "graph.txt"), genes, arms
        )
    sets = {name: [name] for name in samples} | {"all": list(samples)}
    report: dict = {
        "string_links": int(string.nnz // 2),
        "string_genes": int((string.sum(axis=1) > 0).sum()),
        "lists": {},
    }
    for list_name, path in lists.items():
        listed = read_gene_list(path, symbols)
        labels = np.array([g in listed for g in genes])
        report["lists"][list_name] = {"listed": int(labels.sum()), "results": {}}
        for set_name, members in sets.items():
            blocks = {
                "expression": np.column_stack(
                    [expression[f"{n}_{s}"] for n in members for s in ("mean", "sd")]
                ),
                "degree": np.column_stack([features[n].degree for n in members]),
                "shape": np.column_stack([features[n].shape for n in members]),
            }
            result = compare(
                blocks,
                [samples[n] for n in members],
                string,
                [features[n].links for n in members],
                labels,
                groups,
                runs=500,
            )
            report["lists"][list_name]["results"][set_name] = result
            print(list_name, set_name, {k: round(v, 3) for k, v in result["auc"].items()})
    output.mkdir(parents=True, exist_ok=True)
    (output / "baselines.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
