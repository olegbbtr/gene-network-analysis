"""What are the hubs of a co-expression network, under the standard definition too?

Papers that call hubs key regulators mostly use WGCNA, where a gene's connectivity is the sum
over all other genes of |correlation|^6 (unsigned) or ((1 + correlation) / 2)^12 (signed),
its default powers. Per dataset, that connectivity and this repository's link counts are
scored against what a gene is known to be or do: a TF (Lambert 2018), cell-essential
(CEGv2, DepMap), constrained in people (gnomAD LOEUF) and a knockdown with a broad effect on
the transcriptome (Perturb-seq, genome-wide K562 screen). An AUC of 0.5 means hubs are no
more often of that kind; the top TOP hubs are also counted against their expected number.

The powers, lists and TOP were fixed before the first run.
"""

import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from .applications import load
from .constraint import labels as constraint_labels
from .constraint import read_constraint
from .perturbation import BROAD, SCREENS, read_screen
from .prediction import read_gene_list
from .public import DATASETS, read_samples, rows, standardise

UNSIGNED_POWER, SIGNED_POWER = 6, 12
TOP = 100


def connectivity(samples: np.ndarray) -> dict[str, np.ndarray]:
    """WGCNA connectivity per gene; `samples` is samples x genes, standardised per gene."""
    data = samples.astype(np.float32)
    correlation = np.clip(data.T @ data / len(data), -1, 1)
    np.fill_diagonal(correlation, 0)
    signed = ((1 + correlation) / 2) ** SIGNED_POWER
    np.fill_diagonal(signed, 0)
    return {
        "wgcna_unsigned": (np.abs(correlation) ** UNSIGNED_POWER).sum(axis=1),
        "wgcna_signed": signed.sum(axis=1),
    }


def scored(scores: dict[str, np.ndarray], kinds: dict[str, tuple[np.ndarray, np.ndarray]]) -> dict:
    """`kinds` maps a name to (label, genes it is defined for)."""
    result: dict = {}
    for name, (label, known) in kinds.items():
        top = np.argsort(-scores["wgcna_unsigned"][known], kind="stable")[:TOP]
        result[name] = {
            "auc": {s: float(roc_auc_score(label[known], v[known])) for s, v in scores.items()},
            "among_top_hubs": int(label[known][top].sum()),
            "expected_among_top_hubs": float(TOP * label[known].mean()),
        }
    return result


def run(root: Path, gene_table: Path, references: Path, data: Path, output: Path) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    index = {gene: i for i, gene in enumerate(genes)}
    everyone = np.ones(len(genes), dtype=bool)

    def listed(found: set[str]) -> np.ndarray:
        return np.array([g in found for g in genes])

    constrained, usable, _ = constraint_labels(
        read_constraint(references / "gnomad-v2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"),
        genes,
    )
    screen = read_screen(references / "perturb-seq-replogle-2022" / SCREENS["k562-genome-wide"])
    knocked = [
        (index[g], c) for g, c in zip(screen.targets, screen.changed, strict=True) if g in index
    ]
    perturbed = np.zeros(len(genes), dtype=bool)
    changed = np.zeros(len(genes))
    for gene, count in knocked:
        perturbed[gene], changed[gene] = True, count
    kinds = {
        "tf": (listed(set((data / "lambert-2018-tfs.txt").read_text().split())), everyone),
        "cegv2": (
            listed(read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols)),
            everyone,
        ),
        "depmap": (
            listed(
                read_gene_list(
                    references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", context.symbols
                )
            ),
            everyone,
        ),
        "constrained": (constrained, usable),
        "broad_knockdown": (
            perturbed & (changed >= np.quantile(changed[perturbed], 1 - BROAD)),
            perturbed,
        ),
    }
    report: dict = {
        "powers": {"unsigned": UNSIGNED_POWER, "signed": SIGNED_POWER},
        "kinds": {name: int(label[known].sum()) for name, (label, known) in kinds.items()},
        "datasets": {},
    }
    for dataset in DATASETS:
        matrix, _ = read_samples(dataset, references, context.symbols)
        scores = connectivity(standardise(rows(matrix, genes)))
        links = context.links[dataset.name]
        scores["links"] = np.asarray(links.sum(axis=1)).ravel()
        scores["out_degree"] = np.diff(context.adjacency[dataset.name].indptr).astype(float)
        report["datasets"][dataset.name] = scored(scores, kinds)
        print(dataset.name, json.dumps(report["datasets"][dataset.name]), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    (output / "hubs.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
