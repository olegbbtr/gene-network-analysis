"""Degree baselines, label controls and nested regularisation for the poster models."""

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from .applications import load
from .feature_comparison import GLOBAL, GROUP_MODELS, feature_blocks
from .pagerank_discovery import predict_variants
from .prediction import cross_validated, nested_cross_validated, paired_gain, read_gene_list
from .report import revision, sha256
from .robustness import top_hits

PENALTIES = (0.001, 0.01, 0.1, 1.0, 10.0)
AUDIT_MODELS = {
    "degree": ("degree",),
    "degree+positions": ("degree", "degree_rank"),
    "pagerank": ("pagerank",),
    "pagerank_positions": ("pagerank_rank",),
    "core": ("core",),
    "eigenvector": ("eigenvector",),
    "global": GLOBAL,
    "degree+global": ("degree", *GLOBAL),
    "expression+network": GROUP_MODELS["expression+network"],
}


def shuffled_within_groups(labels, groups, seed):
    shuffled = labels.copy()
    rng = np.random.default_rng(seed)
    for group in np.unique(groups):
        selected = np.flatnonzero(groups == group)
        shuffled[selected] = rng.permutation(labels[selected])
    return shuffled


def metrics(known, target, predictions):
    new = ~known
    return {
        name: {
            "known_auc": float(roc_auc_score(known, scores)),
            "discovery_auc": float(roc_auc_score(target[new], scores[new])),
            "top_hits": top_hits(scores[new], target[new]),
        }
        for name, scores in predictions.items()
    }


def correlations(parts, names):
    rows = []
    for i, name in enumerate(names):
        columns = {
            "pagerank_forward": parts["pagerank"][:, 2 * i],
            "pagerank_reverse": parts["pagerank"][:, 2 * i + 1],
            "position_forward": parts["pagerank_rank"][:, 2 * i],
            "position_reverse": parts["pagerank_rank"][:, 2 * i + 1],
            "core": parts["core"][:, i],
            "eigenvector": parts["eigenvector"][:, i],
        }
        for feature, values in columns.items():
            row = {"network": name, "feature": feature}
            for block, suffixes in (("degree", ("out", "in")), ("expression", ("mean", "sd"))):
                for j, suffix in enumerate(suffixes):
                    rho = float(spearmanr(values, parts[block][:, 2 * i + j]).statistic)
                    row[f"{block}_{suffix}"] = rho if np.isfinite(rho) else None
            rows.append(row)
    return rows


def run(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
    from .graph_features import rank_positions

    context = load(root, gene_table)
    parts, links = feature_blocks(context)
    parts["core"] = parts["graph_theory"][:, ::2]
    parts["eigenvector"] = parts["graph_theory"][:, 1::2]
    parts["degree_rank"] = rank_positions(parts["degree"])
    known = np.isin(
        context.genes, list(read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols))
    )
    target = np.isin(
        context.genes,
        list(
            read_gene_list(
                references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", context.symbols
            )
        ),
    )
    fixed = predict_variants(parts, links, known, context.groups, AUDIT_MODELS)
    controls = []
    for seed in range(10):
        shuffled = shuffled_within_groups(known, context.groups, seed)
        scores = cross_validated(parts, GLOBAL, shuffled, context.groups, links)
        controls.append({"seed": seed, "auc": float(roc_auc_score(shuffled, scores))})
        print("shuffled", seed, flush=True)
    tuned, choices = {}, {}
    for name, model in GROUP_MODELS.items():
        tuned[name], choices[name] = nested_cross_validated(
            parts, model, known, context.groups, links, penalties=PENALTIES
        )
        print("tuned", name, flush=True)
    inputs = [
        gene_table,
        root / "expression.tsv",
        root / "stats.json",
        references / "gene-lists/CEGv2.txt",
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv",
        *(root / name / "graph.txt" for name in context.plain),
    ]
    report = {
        "implementation": revision(),
        "inputs": {str(path): sha256(path) for path in inputs},
        "settings": {
            "models": AUDIT_MODELS,
            "tuned_models": GROUP_MODELS,
            "global_columns": sum(parts[key].shape[1] for key in GLOBAL),
            "networks": context.plain,
            "outer_folds": 5,
            "inner_folds": 4,
            "C_grid": PENALTIES,
            "selection_metric": "pooled inner-held-out CEGv2 ROC AUC; first C wins ties",
            "scaling": "fit on each inner/outer training fold only",
            "neighbours": "recomputed from each inner/outer training fold's labels only",
            "null": "ten fixed seeds; permute CEGv2 labels within chromosomes; C=1",
            "scope": "unseen gene labels on observed graphs; held-out nodes remain in topology",
        },
        "fixed": metrics(known, target, fixed),
        "shuffled_label_controls": controls,
        "correlations": correlations(parts, context.plain),
        "nested": metrics(known, target, tuned),
        "selected_C": choices,
        "global_beyond_degree": paired_gain(known, fixed["degree+global"], fixed["degree"]),
        "nested_all_beyond_gene_local": {
            "cegv2": paired_gain(
                known, tuned["expression+network+global"], tuned["expression+network"]
            ),
            "depmap": paired_gain(
                target[~known],
                tuned["expression+network+global"][~known],
                tuned["expression+network"][~known],
            ),
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "feature-audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ("fixed", "nested", "shuffled_label_controls")}))
    return report
