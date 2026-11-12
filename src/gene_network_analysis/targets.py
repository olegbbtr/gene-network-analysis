"""Harder targets for the label-free feature groups of `label_free`.

Each target is a binary label on a universe of genes. The same feature groups, nested
chromosome-fold penalty selection, models (each group alone, all four, and all four
without gene, graph or biological features) and leave-one-out impact as `label_free` are
used. All definitions below were fixed before the first run.

- selective_dependency: genes outside CEGv2 and DepMap 24Q2 common essentials, screened in
  at least MIN_LINES DepMap lines. Positive when at least SELECTIVE_SHARE of the screened
  lines need the gene (gene effect at most DEPENDENCY).
- knockdown_impact: genes knocked down in the K562 genome-wide Perturb-seq screen (Replogle
  2022; one row per gene as in `perturbation`), outside CEGv2 and DepMap common essentials,
  whose knockdown removed at least half of the target's expression. Positive when the
  knockdown changed at least as many genes as the top BROAD share of all K562 genome-wide
  knockdowns of network genes (the `perturbation` threshold).
- first_function: genes with no focused experimental Gene Ontology annotation dated before
  CUTOFF (PubMed ID, about August 2018). Positive when they gained one from CUTOFF on.
  Experimental codes: EXP, IDA, IPI, IMP, IGI (high-throughput codes left out); a reference
  is focused when the GAF cites it for at most FOCUSED genes, so large screens and
  interactome maps do not count. Any aspect (process, function, component); NOT rows are
  left out. Because well-studied genes attract more papers, the gene group of this target
  also holds log(1 + distinct PubMed IDs citing the gene before CUTOFF, any evidence).

Added on 2026-11-12 after the three targets above, because the graph added little to them;
chosen where earlier analyses (Q29, Q30) suggested structure matters, and fixed before
their first run here:

- constraint: genes outside CEGv2 and DepMap common essentials with at least MIN_EXPECTED
  expected loss-of-function variants in gnomAD v2.1.1. Positive when LOEUF < LOEUF_CUTOFF.
  The gene group also holds log expected variant count, so gene length is not credited to
  the network.
- understudied_essential: genes in the lowest tercile of STRING v12 links (score >= 700,
  among the 12,000 genes; as in `understudied`). Positive when a DepMap common essential.

Added later on 2026-11-12, fixed before its first run:

- tumour_selective: DepMap 24Q2 lines are split into cancer lines and the 10 lines DepMap
  calls Non-Cancerous; the five RPE1 clones are averaged into one, giving HEALTHY_GROUPS
  healthy groups. Genes screened in every healthy group and in at least MIN_LINES cancer
  lines. Positive when at least TUMOUR_SHARE of the cancer lines need the gene (effect at
  most DEPENDENCY) and at most HEALTHY_ALLOWED healthy groups do.
"""

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed
from sklearn.metrics import average_precision_score, roc_auc_score

from . import cutoff_counts
from .baselines import read_string
from .constraint import labels as constraint_labels
from .constraint import read_constraint
from .gene_sets import DEPENDENCY
from .label_free import ALL, LABELS, LAYOUT, MODELS, PENALTIES, REACTOME, WITHOUT, build, removal
from .perturbation import BROAD, SCREENS, read_screen
from .prediction import nested_cross_validated, read_gene_list
from .prospective import CUTOFF
from .references import GeneMap, read_gene_table
from .report import revision, sha256
from .robustness import read_effects
from .understudied import terciles

MIN_LINES = 1000
SELECTIVE_SHARE = 0.05
KNOCKDOWN = -0.5
FOCUSED = 20
CODES = frozenset({"EXP", "IDA", "IPI", "IMP", "IGI"})
TOP = 100
TUMOUR_SHARE = 0.25
HEALTHY_ALLOWED = 1
HEALTHY_GROUPS = 6


def common_essentials(references: Path, genes: list[str], symbols: dict[str, str]) -> np.ndarray:
    listed = read_gene_list(references / "gene-lists/CEGv2.txt", symbols)
    listed |= read_gene_list(references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", symbols)
    return np.array([g in listed for g in genes])


def selective_dependency(references, genes, symbols, common) -> dict:
    effects = read_effects(references / "depmap-24q2/CRISPRGeneEffect.csv", symbols, genes)
    screened = np.isfinite(effects).sum(axis=0)
    share = (effects <= DEPENDENCY).sum(axis=0) / np.maximum(screened, 1)
    universe = ~common & (screened >= MIN_LINES)
    return {"universe": universe, "labels": share >= SELECTIVE_SHARE, "extra": None}


def knockdown_impact(references, genes, symbols, common) -> dict:
    index = {g: i for i, g in enumerate(genes)}
    screen = read_screen(references / "perturb-seq-replogle-2022" / SCREENS["k562-genome-wide"])
    rows = [r for r, g in enumerate(screen.targets) if g in index]
    perturbed = np.zeros(len(genes), dtype=bool)
    changed = np.zeros(len(genes))
    efficiency = np.full(len(genes), np.nan)
    for r in rows:
        i = index[screen.targets[r]]
        perturbed[i] = True
        changed[i] = screen.changed[r]
        efficiency[i] = screen.covariates[r, 0]
    threshold = float(np.quantile(changed[perturbed], 1 - BROAD))
    effective = perturbed & np.isfinite(efficiency) & (efficiency <= KNOCKDOWN)
    return {
        "universe": effective & ~common,
        "labels": changed >= threshold,
        "extra": None,
        "threshold_changed_genes": threshold,
    }


def first_function(references, genes, symbols, common) -> dict:
    lookup = GeneMap(references / "hgnc/hgnc_complete_set.txt").by_symbol
    rows = []
    cited = defaultdict(set)
    with gzip.open(references / "gene-ontology/goa_human.gaf.gz", "rt") as stream:
        for line in stream:
            if line.startswith("!"):
                continue
            row = line.rstrip("\n").split("\t")
            gene = lookup(row[2])
            pmids = [int(r[5:]) for r in row[5].split("|") if r.startswith("PMID:")]
            if not gene or not pmids or "NOT" in row[3]:
                continue
            for pmid in pmids:
                cited[pmid].add(gene)
            rows.append((gene, min(pmids), pmids, row[6]))
    focused = {pmid for pmid, members in cited.items() if len(members) <= FOCUSED}
    before, after, prior = set(), set(), defaultdict(set)
    for gene, first, pmids, evidence in rows:
        if first < CUTOFF:
            prior[gene].update(p for p in pmids if p < CUTOFF)
        if evidence in CODES and any(p in focused for p in pmids):
            (before if first < CUTOFF else after).add(gene)
    universe = np.array([g not in before for g in genes])
    labels = np.array([g in after for g in genes])
    papers = np.log1p([len(prior.get(g, ())) for g in genes])[:, None]
    return {"universe": universe, "labels": labels, "extra": papers}


def constraint(references, genes, symbols, common) -> dict:
    table = read_constraint(references / "gnomad-v2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
    constrained, usable, expected = constraint_labels(table, genes)
    return {"universe": usable & ~common, "labels": constrained, "extra": expected[:, None]}


def understudied_essential(references, genes, symbols, common) -> dict:
    string = read_string(
        references / "string-v12/9606.protein.info.v12.0.txt.gz",
        references / "string-v12/9606.protein.links.v12.0.txt.gz",
        genes,
        GeneMap(references / "hgnc/hgnc_complete_set.txt").by_symbol,
    )
    group = terciles(np.asarray(string.sum(axis=1)).ravel())
    depmap = read_gene_list(references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", symbols)
    return {
        "universe": group == 0,
        "labels": np.array([g in depmap for g in genes]),
        "extra": None,
    }


def tumour_selective(references, genes, symbols, common) -> dict:
    path = references / "depmap-24q2/CRISPRGeneEffect.csv"
    effects = read_effects(path, symbols, genes)
    with path.open() as stream:
        next(stream)
        lines = [row.split(",", 1)[0] for row in stream]
    with (references / "depmap-24q2/Model.csv").open(newline="") as stream:
        models = {r["ModelID"]: r for r in csv.DictReader(stream)}
    healthy = [
        i for i, m in enumerate(lines)
        if models.get(m, {}).get("OncotreePrimaryDisease") == "Non-Cancerous"
    ]  # fmt: skip
    names = [models[lines[i]].get("CellLineName", "") for i in healthy]
    rpe1 = [i for i, n in zip(healthy, names, strict=True) if n.startswith("RPE1")]
    groups = [[i] for i in healthy if i not in rpe1] + [rpe1]
    if len(groups) != HEALTHY_GROUPS:
        raise ValueError(f"expected {HEALTHY_GROUPS} healthy groups, found {len(groups)}")
    healthy_effect = np.column_stack([np.nanmean(effects[g], axis=0) for g in groups])
    cancer = np.setdiff1d(np.arange(len(lines)), healthy)
    cancer_effects = effects[cancer]
    screened = np.isfinite(cancer_effects).sum(axis=0)
    share = (cancer_effects <= DEPENDENCY).sum(axis=0) / np.maximum(screened, 1)
    healthy_needed = (healthy_effect <= DEPENDENCY).sum(axis=1)
    universe = (screened >= MIN_LINES) & np.isfinite(healthy_effect).all(axis=1)
    labels = (share >= TUMOUR_SHARE) & (healthy_needed <= HEALTHY_ALLOWED)
    return {"universe": universe, "labels": labels, "extra": None}


TARGETS = {
    "selective_dependency": selective_dependency,
    "knockdown_impact": knockdown_impact,
    "first_function": first_function,
    "constraint": constraint,
    "understudied_essential": understudied_essential,
    "tumour_selective": tumour_selective,
}


def evaluate(labels: np.ndarray, predictions: dict[str, np.ndarray]) -> dict:
    order = np.argsort(-predictions[ALL])
    return {
        "models": {
            k: {
                "auc": float(roc_auc_score(labels, p)),
                "average_precision": float(average_precision_score(labels, p)),
                "top_hits": int(labels[np.argsort(-p)[:TOP]].sum()),
            }
            for k, p in predictions.items()
        },
        "base_rate": float(labels.mean()),
        "chance_top_hits": TOP * float(labels.mean()),
        "leave_one_out": removal(labels, predictions),
        "top_genes": order[:TOP].tolist(),
    }


def run(root: Path, gene_table: Path, references: Path, output: Path, names=None) -> dict:
    stats = json.loads((root / "stats.json").read_text())
    networks = [n for n in stats if n != "settings"]
    paths = [
        gene_table,
        root / "expression.tsv",
        root / "stats.json",
        references / "gene-lists/CEGv2.txt",
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv",
        references / REACTOME,
        root / cutoff_counts.CACHE,
        *(root / n / "graph.txt" for n in networks if not n.endswith("-copy-number-removed")),
    ]
    inputs = {str(path): sha256(path) for path in paths} | {"layout": LAYOUT}
    genes, groups, parts = build(root, gene_table, references, inputs)
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    common = common_essentials(references, genes, symbols)
    path = output / "targets.json"
    report = json.loads(path.read_text()) if path.exists() else {}
    report["implementation"] = revision()
    for name in names or TARGETS:
        target = TARGETS[name](references, genes, symbols, common)
        keep = target["universe"]
        labels = target["labels"][keep]
        blocks = {p: block[keep] for p, block in parts.items()}
        if target["extra"] is not None:
            blocks["expression"] = np.column_stack([blocks["expression"], target["extra"][keep]])
        fitted = Parallel(n_jobs=len(MODELS))(
            delayed(nested_cross_validated)(
                blocks, model, labels, groups[keep], [], penalties=PENALTIES
            )
            for model in MODELS.values()
        )
        predictions = {k: result[0] for k, result in zip(MODELS, fitted, strict=True)}
        selected = {
            k: [f["C"] for f in result[1]] for k, result in zip(MODELS, fitted, strict=True)
        }
        result = evaluate(labels, predictions)
        kept = np.flatnonzero(keep)
        result["top_genes"] = [
            table.get(genes[kept[i]], {}).get("symbol", genes[kept[i]]) for i in result["top_genes"]
        ]
        result["genes"] = int(keep.sum())
        result["positives"] = int(labels.sum())
        result["labels"] = LABELS
        result["without"] = WITHOUT
        result["selected_C"] = selected
        result["C_grid"] = PENALTIES
        for key in ("threshold_changed_genes",):
            if key in target:
                result[key] = target[key]
        report[name] = result
        print(name, json.dumps({k: v for k, v in result.items() if k != "top_genes"}), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n")
    return report
