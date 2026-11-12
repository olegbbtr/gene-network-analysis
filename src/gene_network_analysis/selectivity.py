"""Which genes look more essential in tumour networks than in matched normal networks?

The label-free all-four model (without the percentile link counts, which exist only for the
main networks) is fitted on CEGv2 twice with the nested chromosome folds of `label_free`:
once on features of the 8 TCGA tumour networks, once on features of their 8 TCGA
adjacent-normal networks (same processing, unlike GTEx). Each gene's selectivity is the
difference of its held-out score percentiles, tumour minus normal.

Check, fixed before the first run: within the `targets.tumour_selective` universe, the ROC
AUC of the selectivity score for tumour-selective essential genes (needed by at least 25%
of DepMap cancer lines, by at most one of 6 healthy groups), with a gene-bootstrap 95%
interval, next to the AUC of the tumour and normal scores alone.
"""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

from .applications import load
from .consensus import together
from .feature_comparison import feature_blocks
from .gene_sets import membership, read_gmt, restrict
from .label_free import ALL, MODELS, PAIR_SIZES, PENALTIES, REACTOME, coherence, essential_labels
from .prediction import SHAPE, nested_cross_validated
from .references import read_gene_table
from .report import revision
from .targets import common_essentials, tumour_selective

TUMOURS = ["luad", "lusc", "brca", "coad", "thca", "prad", "stad", "lihc"]
BOOTSTRAP = 1000
TOP = 25


def parts_for(context, names: list[str], reactome) -> dict[str, np.ndarray]:
    sub = replace(
        context,
        adjacency={n: context.adjacency[n] for n in names},
        links={n: context.links[n] for n in names},
    )
    parts, links = feature_blocks(sub)
    arm = [i for i in range(parts["shape"].shape[1]) if i % len(SHAPE) == SHAPE.index("arm_share")]
    parts["arm_share"] = parts["shape"][:, arm]
    parts["shape"] = np.delete(parts["shape"], arm, axis=1)
    genes = context.genes
    parts["pathway_coherence"] = coherence(
        links,
        [context.adjacency[n] for n in names],
        together(genes, reactome),
        membership(genes, restrict(reactome, genes, *PAIR_SIZES)).any(axis=1).astype(float),
    )
    return parts


def run(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    known, _ = essential_labels(references, genes, symbols)
    reactome = read_gmt(references / REACTOME, context.symbols)
    model = tuple(p for p in MODELS[ALL] if p != "cutoff_counts")
    scores = {}
    for side, names in {
        "tumour": [f"tcga-{t}" for t in TUMOURS],
        "normal": [f"tcga-{t}-normal" for t in TUMOURS],
    }.items():
        parts = parts_for(context, names, reactome)
        scores[side] = nested_cross_validated(
            parts, model, known, context.groups, [], penalties=PENALTIES
        )[0]
        print(side, "CEGv2 AUC", roc_auc_score(known, scores[side]), flush=True)
    percentile = {k: rankdata(v) / len(v) for k, v in scores.items()}
    selectivity = percentile["tumour"] - percentile["normal"]
    common = common_essentials(references, genes, symbols)
    target = tumour_selective(references, genes, symbols, common)
    keep, labels = target["universe"], target["labels"]
    rng = np.random.default_rng(0)
    index = np.flatnonzero(keep)
    samples = []
    while len(samples) < BOOTSTRAP:
        chosen = rng.choice(index, len(index))
        if labels[chosen].any():
            samples.append(roc_auc_score(labels[chosen], selectivity[chosen]))
    order = np.argsort(-selectivity)
    report = {
        "implementation": revision(),
        "cegv2_auc": {k: float(roc_auc_score(known, v)) for k, v in scores.items()},
        "tumour_selective": {
            "genes": int(keep.sum()),
            "positives": int(labels[keep].sum()),
            "auc_selectivity": float(roc_auc_score(labels[keep], selectivity[keep])),
            "auc_selectivity_ci95": np.quantile(samples, [0.025, 0.975]).tolist(),
            "auc_tumour_score": float(roc_auc_score(labels[keep], scores["tumour"][keep])),
            "auc_normal_score": float(roc_auc_score(labels[keep], scores["normal"][keep])),
        },
        "top": [
            {
                "symbol": table.get(genes[i], {}).get("symbol", genes[i]),
                "selectivity": float(selectivity[i]),
                "tumour_percentile": float(percentile["tumour"][i]),
                "normal_percentile": float(percentile["normal"][i]),
                "in_cegv2": bool(known[i]),
                "tumour_selective": bool(labels[i]) if keep[i] else None,
            }
            for i in order[:TOP]
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "selectivity.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report
