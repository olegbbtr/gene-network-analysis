"""Continuous version of the understudied-essential target, fixed before its first run.

The binary target in `targets` (DepMap common essential among genes with 0-4 STRING links)
has only 96 positives, so its leave-one-out intervals are wide. Here the same universe
(lowest STRING tercile, screened in at least MIN_LINES DepMap lines) is kept and the label is
the gene's mean DepMap gene effect across screened lines, negated so higher means more
needed. Models are ridge regressions on the label-free groups of `label_free`, with five
outer chromosome folds and four inner folds choosing the penalty by Spearman correlation;
scaling uses each training fold only. Score: Spearman correlation on held-out genes.
Leave-one-out: the drop in that correlation when gene, graph (local + global) or biological
features are removed, with gene-bootstrap 95% intervals over the fixed held-out predictions.
"""

import json
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from . import cutoff_counts
from .baselines import read_string
from .label_free import ALL, BLOCKS, LABELS, LAYOUT, MODELS, REACTOME, WITHOUT, build
from .references import GeneMap, read_gene_table
from .report import revision, sha256
from .robustness import read_effects
from .targets import MIN_LINES
from .understudied import terciles

ALPHAS = (0.1, 1.0, 10.0, 100.0, 1000.0)
BOOTSTRAP = 1000


def rho(y: np.ndarray, p: np.ndarray) -> float:
    return float(spearmanr(y, p).statistic)


def nested_ridge(x: np.ndarray, y: np.ndarray, groups: np.ndarray) -> np.ndarray:
    predictions = np.zeros(len(y))
    for train, test in GroupKFold(5).split(x, groups=groups):
        inner = np.zeros((len(ALPHAS), len(train)))
        for a, b in GroupKFold(4).split(train, groups=groups[train]):
            scaler = StandardScaler().fit(x[train[a]])
            for i, alpha in enumerate(ALPHAS):
                model = Ridge(alpha=alpha).fit(scaler.transform(x[train[a]]), y[train[a]])
                inner[i, b] = model.predict(scaler.transform(x[train[b]]))
        best = ALPHAS[int(np.argmax([rho(y[train], s) for s in inner]))]
        scaler = StandardScaler().fit(x[train])
        model = Ridge(alpha=best).fit(scaler.transform(x[train]), y[train])
        predictions[test] = model.predict(scaler.transform(x[test]))
    return predictions


def loss(y: np.ndarray, predictions: dict[str, np.ndarray], seed: int = 0) -> dict:
    def losses(chosen):
        full = rho(y[chosen], predictions[ALL][chosen])
        return {b: full - rho(y[chosen], predictions[m][chosen]) for b, m in WITHOUT.items()}

    values = losses(np.arange(len(y)))
    rng = np.random.default_rng(seed)
    samples = [losses(rng.integers(0, len(y), len(y))) for _ in range(BOOTSTRAP)]
    return {
        b: {
            "removed": BLOCKS[b],
            "rho_loss": values[b],
            "rho_loss_ci95": np.quantile([s[b] for s in samples], [0.025, 0.975]).tolist(),
        }
        for b in BLOCKS
    }


def run(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
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
    string = read_string(
        references / "string-v12/9606.protein.info.v12.0.txt.gz",
        references / "string-v12/9606.protein.links.v12.0.txt.gz",
        genes,
        GeneMap(references / "hgnc/hgnc_complete_set.txt").by_symbol,
    )
    low = terciles(np.asarray(string.sum(axis=1)).ravel()) == 0
    effects = read_effects(references / "depmap-24q2/CRISPRGeneEffect.csv", symbols, genes)
    screened = np.isfinite(effects).sum(axis=0)
    keep = low & (screened >= MIN_LINES)
    y = -np.nanmean(effects[:, keep], axis=0)
    blocks = {p: b[keep] for p, b in parts.items()}
    fitted = Parallel(n_jobs=len(MODELS))(
        delayed(nested_ridge)(np.column_stack([blocks[p] for p in model]), y, groups[keep])
        for model in MODELS.values()
    )
    predictions = dict(zip(MODELS, fitted, strict=True))
    report = {
        "implementation": revision(),
        "genes": int(keep.sum()),
        "models": {k: rho(y, p) for k, p in predictions.items()},
        "labels": LABELS,
        "leave_one_out": loss(y, predictions),
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "graded-understudied.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report
