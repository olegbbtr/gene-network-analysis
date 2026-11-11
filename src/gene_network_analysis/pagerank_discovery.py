"""Fixed PageRank additions to the CEGv2-to-DepMap discovery model."""

import importlib.metadata
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from sklearn.metrics import roc_auc_score

from .applications import load
from .graph import Graph
from .graph_features import rank_positions
from .prediction import cross_validated, network_features, paired_gain, read_gene_list
from .report import revision, sha256
from .robustness import NETWORK, blocks, discovery, top_indices
from .structure import ARM, location_codes

MODELS = {
    "expression": ("expression",),
    "expression+network": NETWORK,
    "expression+network+pagerank": (*NETWORK, "pagerank"),
    "expression+network+personalised": (*NETWORK, "personalised"),
}
JOINT_MODELS = {
    "expression": ("expression",),
    "network": ("degree", "shape", "neighbours"),
    "expression+network": NETWORK,
    "expression+network+pagerank": (*NETWORK, "pagerank"),
    "expression+network+pagerank-rank": (*NETWORK, "pagerank_rank"),
    "expression+network+pagerank-both": (*NETWORK, "pagerank", "pagerank_rank"),
    "network+pagerank-both": ("degree", "shape", "neighbours", "pagerank", "pagerank_rank"),
}


def predict_variants(
    parts: dict[str, np.ndarray],
    links: list[csr_array],
    known: np.ndarray,
    groups: np.ndarray,
    models: dict[str, tuple[str, ...]] | None = None,
) -> dict[str, np.ndarray]:
    predictions = {}
    for name, model in (MODELS if models is None else models).items():
        predictions[name] = cross_validated(parts, model, known, groups, links)
        print("fitted", name, flush=True)
    return predictions


def run(
    root: Path,
    gene_table: Path,
    references: Path,
    output: Path,
    runs: int = 1000,
    *,
    models: dict[str, tuple[str, ...]] | None = None,
    output_name: str = "pagerank-discovery",
) -> dict:
    models = MODELS if models is None else models
    context = load(root, gene_table)
    names = context.plain
    genes = context.genes
    arms = location_codes(tuple(genes), context.table, ARM)
    features = {}
    for name in names:
        features[name] = network_features(Graph(tuple(genes), context.adjacency[name]), genes, arms)
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
    links = [features[name].links for name in names]
    predictions = predict_variants(parts, links, known, context.groups, models)
    result = discovery(known, target, predictions, runs)
    new = np.flatnonzero(~known)
    baseline = predictions["expression+network"][new]
    result["gain_over_network"] = {
        name: paired_gain(target[new], predictions[name][new], baseline, runs)
        for name in models
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
    inputs = [
        gene_table,
        root / "expression.tsv",
        root / "stats.json",
        references / "gene-lists/CEGv2.txt",
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv",
        *(root / name / "graph.txt" for name in names),
    ]
    report = {
        "implementation": revision(),
        "versions": {
            name: importlib.metadata.version(name) for name in ("numpy", "scipy", "scikit-learn")
        },
        "inputs": {str(path): sha256(path) for path in inputs},
        "settings": {
            "networks": names,
            "models": models,
            "folds": 5,
            "split": "GroupKFold by chromosome",
            "bootstrap_runs": runs,
            "bootstrap_seed": 0,
            "pagerank": "log(genes * score), unweighted forward and reverse, alpha=0.85",
            "pagerank_rank": "top-first positions per column, tied ranks averaged, in [0,1]",
            "personalised": "training-fold CEGv2 seeds, undirected links, alpha=0.85, 50 steps",
            "personalised_training": "inner chromosome folds hide each training gene's label",
            "scope": "combined 15 original public networks; no tuning on DepMap outcomes",
        },
        "genes": len(genes),
        "training_essentials": int(known.sum()),
        "known_essential_auc": {
            name: float(roc_auc_score(known, p)) for name, p in predictions.items()
        },
        "discovery": result,
    }
    if not any("personalised" in parts for parts in models.values()):
        for key in ("personalised", "personalised_training"):
            report["settings"].pop(key)
    output.mkdir(parents=True, exist_ok=True)
    (output / f"{output_name}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "top_candidates"}))
    return report
