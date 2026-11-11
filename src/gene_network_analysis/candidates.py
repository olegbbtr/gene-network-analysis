"""Audit top-ranked genes outside both essential lists against DepMap gene effects."""

import json
from pathlib import Path

import numpy as np

from .gene_sets import DEPENDENCY
from .prediction import read_gene_list
from .references import read_gene_table
from .report import revision, sha256
from .robustness import read_effects


def summarise(
    ranking: list[dict], known: set[str], common: set[str], effects: np.ndarray
) -> list[dict]:
    if effects.shape[1] != len(ranking):
        raise ValueError("gene-effect columns must match the ranking order")
    result = []
    for column, gene in enumerate(ranking):
        if gene["depmap"] != (gene["gene"] in common):
            raise ValueError("DepMap membership differs between prediction and reference")
        if gene["gene"] in known or gene["gene"] in common:
            continue
        values = effects[:, column]
        measured = values[np.isfinite(values)]
        needed = int((measured <= DEPENDENCY).sum())
        result.append(
            {
                **gene,
                "rank_outside_cegv2": column + 1,
                "screened_lines": len(measured),
                "dependent_lines": needed,
                "dependent_share": needed / len(measured) if len(measured) else None,
                "median_gene_effect": float(np.median(measured)) if len(measured) else None,
                "in_cegv2": False,
                "in_depmap_common": False,
            }
        )
    return result


def run(predictions: Path, gene_table: Path, references: Path, output: Path) -> dict:
    source = json.loads(predictions.read_text())
    model = "expression+network"
    ranking = source["discovery"]["top_candidates"][model]
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    ceg_path = references / "gene-lists/CEGv2.txt"
    common_path = references / "depmap-24q2/CRISPRInferredCommonEssentials.csv"
    effects_path = references / "depmap-24q2/CRISPRGeneEffect.csv"
    known, common = (read_gene_list(p, symbols) for p in (ceg_path, common_path))
    effects = read_effects(effects_path, symbols, [row["gene"] for row in ranking])
    candidates = summarise(ranking, known, common, effects)
    report = {
        "implementation": revision(),
        "inputs": {
            str(p): sha256(p)
            for p in (predictions, gene_table, ceg_path, common_path, effects_path)
        },
        "settings": {
            "model": model,
            "selection": "all top-100 predictions outside CEGv2 and DepMap common essentials",
            "dependency_threshold": DEPENDENCY,
            "scope": "absence from these lists is not absence from the literature",
        },
        "screened_cell_lines": effects.shape[0],
        "candidates": candidates,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "candidate-evidence.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(candidates, indent=2))
    return report
