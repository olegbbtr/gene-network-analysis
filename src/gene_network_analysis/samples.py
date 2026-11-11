"""How many samples does a network need before its gene-level uses appear?

For a few large datasets, networks are rebuilt from random subsets of samples with the same
inference as `public.build`. Each is scored on essential genes (gain over expression, as in
`prediction`) and on gene function (median held-out AUC over Reactome pathways).
"""

import json
from pathlib import Path

import numpy as np

from .diagnostics import chromosome
from .gene_sets import held_out_scores, membership, read_gmt, restrict, set_aucs
from .graph import Graph
from .prediction import evaluate, network_features, read_gene_list
from .public import DATASETS, read_expression, read_samples, rows, zscores
from .references import read_gene_table
from .simulate import infer
from .structure import ARM, location_codes

USED = ("tcga-brca", "tcga-luad", "gtex-thyroid", "gtex-lung")
SIZES = (25, 50, 100, 200)
SEEDS = (0, 1)
THRESHOLD = 2.35


def score(
    values: np.ndarray,
    genes: list[str],
    arms: np.ndarray,
    groups: np.ndarray,
    essential: np.ndarray,
    pathways: np.ndarray,
) -> dict:
    """`values` is genes x samples; genes constant in a subset get no links."""
    adjacency = infer(zscores(values), np.arange(len(genes)), THRESHOLD)
    features = network_features(Graph(tuple(genes), adjacency), genes, arms)
    expression = np.column_stack([values.mean(axis=1), values.std(axis=1)])
    result = evaluate(expression, [features], essential, groups, runs=200)
    function = set_aucs(held_out_scores(features.links, pathways, groups), pathways)
    return {
        "edges": int(adjacency.nnz),
        "essential_auc": {
            k: result["auc"][k] for k in ("expression", "neighbours", "expression+shape+neighbours")
        },
        "essential_gains": {
            k: result["gains"][k] for k in ("network_beyond_expression", "neighbours_beyond_shape")
        },
        "function_median_auc": float(np.median(function)),
    }


def run(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    genes, _ = read_expression(root / "expression.tsv")
    arms = location_codes(tuple(genes), table, ARM)
    groups = np.array([chromosome(table.get(g, {}).get("location", "")) for g in genes])
    listed = read_gene_list(references / "gene-lists/CEGv2.txt", symbols)
    essential = np.array([g in listed for g in genes])
    reactome = restrict(read_gmt(references / "reactome/ReactomePathways.gmt.zip", symbols), genes)
    pathways = membership(genes, reactome)
    report: dict = {"pathways": len(reactome), "datasets": {}}
    for dataset in (d for d in DATASETS if d.name in USED):
        matrix, _ = read_samples(dataset, references, symbols)
        values = rows(matrix, genes)
        total = values.shape[1]
        entries = []
        for size in [s for s in SIZES if s < total]:
            for seed in SEEDS:
                chosen = np.random.default_rng(seed).choice(total, size, replace=False)
                entry = {"samples": size, "seed": seed} | score(
                    values[:, chosen], genes, arms, groups, essential, pathways
                )
                entries.append(entry)
                print(dataset.name, json.dumps(entry))
        entries.append(
            {"samples": total, "seed": None}
            | score(values, genes, arms, groups, essential, pathways)
        )
        print(dataset.name, json.dumps(entries[-1]))
        report["datasets"][dataset.name] = entries
    output.mkdir(parents=True, exist_ok=True)
    (output / "samples.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
