"""Gene, neighbourhood and structural descriptors, with fixed PageRank replacements."""

import importlib.metadata
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from .applications import load
from .graph import Graph
from .graph_features import core_depth, eigenvector_score, rank_positions, weight_summaries
from .pagerank_discovery import predict_variants
from .prediction import network_features, paired_gain, read_gene_list
from .report import revision, sha256
from .robustness import NETWORK, blocks, discovery, top_indices
from .structure import ARM, location_codes

MAIN_MODELS = {
    "expression": "Expression only",
    "network": "Network features only",
    "expression+network": "Expression + network features",
}
GROUP_LABELS = {
    "expression": "Gene only",
    "network": "Local graph only",
    "global": "Global graph only",
    "expression+network": "Gene + local graph",
    "expression+network+global": "Gene + local graph + global graph",
}
GLOBAL = ("pagerank", "pagerank_rank", "graph_theory")
GROUP_MODELS = {
    "expression": ("expression",),
    "network": ("degree", "shape", "neighbours"),
    "global": GLOBAL,
    "expression+network": NETWORK,
    "expression+network+global": (*NETWORK, *GLOBAL),
}
MODELS = {
    "expression": ("expression",),
    "network": ("degree", "shape", "neighbours"),
    "expression+network": NETWORK,
    "expression+network+pagerank": (*NETWORK, "pagerank"),
    "expression+network+pagerank-rank": (*NETWORK, "pagerank_rank"),
    "network+graph-theory": ("degree", "shape", "neighbours", "graph_theory"),
    "expression+network+graph-theory": (*NETWORK, "graph_theory"),
    "expression+network+graph-theory+pagerank-rank": (*NETWORK, "graph_theory", "pagerank_rank"),
    "network+graph-theory+pagerank-rank": (
        "degree",
        "shape",
        "neighbours",
        "graph_theory",
        "pagerank_rank",
    ),
    "expression+network+weight-summaries": (*NETWORK, "weights"),
}
for removed in ("degree", "shape", "neighbours"):
    kept = tuple(part for part in NETWORK if part != removed)
    MODELS[f"without-{removed}"] = kept
    MODELS[f"without-{removed}+pagerank"] = (*kept, "pagerank")
    MODELS[f"without-{removed}+pagerank-rank"] = (*kept, "pagerank_rank")


def run(
    root: Path,
    gene_table: Path,
    references: Path,
    output: Path,
    runs: int = 1000,
    *,
    grouped: bool = False,
) -> dict:
    models = GROUP_MODELS if grouped else MODELS
    main_models = GROUP_LABELS if grouped else MAIN_MODELS
    context = load(root, gene_table)
    names, genes = context.plain, context.genes
    arms = location_codes(tuple(genes), context.table, ARM)
    features, extra, weights = {}, {}, {}
    for name in names:
        matrix = context.adjacency[name]
        features[name] = network_features(Graph(tuple(genes), matrix), genes, arms)
        links = features[name].links
        extra[name] = np.column_stack(
            [np.log1p(core_depth(links)), np.log1p(len(genes) * eigenvector_score(links))]
        )
        weights[name] = weight_summaries(matrix)
        print("features", name, flush=True)
    known = np.isin(
        genes, list(read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols))
    )
    target = np.isin(
        genes,
        list(
            read_gene_list(
                references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", context.symbols
            )
        ),
    )
    parts = blocks(context, features, names)
    parts["pagerank_rank"] = rank_positions(parts["pagerank"])
    parts["graph_theory"] = np.column_stack([extra[name] for name in names])
    parts["weights"] = np.column_stack([weights[name] for name in names])
    predictions = predict_variants(
        parts, [features[name].links for name in names], known, context.groups, models
    )
    result = discovery(known, target, predictions, runs)
    new = np.flatnonzero(~known)
    result["gain_over_network"] = {
        name: paired_gain(target[new], scores[new], predictions["expression+network"][new], runs)
        for name, scores in predictions.items()
        if name not in {"expression", "expression+network"}
    }
    result["top_candidates"] = {
        name: [
            {
                "gene": genes[i],
                "symbol": context.table[genes[i]]["symbol"],
                "depmap": bool(target[i]),
                "score": float(scores[i]),
            }
            for i in new[top_indices(scores[new])]
        ]
        for name, scores in predictions.items()
    }
    known_auc = {name: float(roc_auc_score(known, p)) for name, p in predictions.items()}
    report = {
        "implementation": revision(),
        "versions": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "scikit-learn")},
        "genes": len(genes),
        "training_essentials": int(known.sum()),
        "settings": {
            "networks": names,
            "models": models,
            "main_models": main_models,
            "folds": 5,
            "split": "GroupKFold by chromosome",
            "bootstrap_runs": runs,
            "bootstrap_seed": 0,
            "pagerank": "forward/reverse unweighted scores, alpha=0.85; log(genes * score)",
            "pagerank_rank": "top-first rank per column, ties averaged, scaled to [0,1]",
            "graph_theory": "log(1+core number), log(1+genes*eigenvector score); undirected",
            "neighbours": "training-fold CEGv2 neighbour share; no held-out labels",
            "weights": "incoming/outgoing mean, SD, absolute mean and positive share",
            "scope": (
                "combined 15 public networks; variants fixed before this run; "
                "no personalised PageRank"
            ),
        },
        "inputs": {
            str(path): sha256(path)
            for path in [
                gene_table,
                root / "expression.tsv",
                root / "stats.json",
                references / "gene-lists/CEGv2.txt",
                references / "depmap-24q2/CRISPRInferredCommonEssentials.csv",
                *(root / name / "graph.txt" for name in names),
            ]
        },
        "known_essential_auc": known_auc,
        "main_comparison": [
            {
                "key": key,
                "label": label,
                "known_auc": known_auc[key],
                "discovery_auc": result["auc"][key],
                "top_hits": result["top_hits"][key],
            }
            for key, label in main_models.items()
        ],
        "discovery": result,
    }
    output.mkdir(parents=True, exist_ok=True)
    filename = "feature-groups.json" if grouped else "feature-comparison.json"
    (output / filename).write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "main_comparison": report["main_comparison"],
                "top_hits": result["top_hits"],
                "gains": result["gain_over_network"],
            }
        )
    )
    return report
