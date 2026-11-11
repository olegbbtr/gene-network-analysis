"""Can copy-number changes in lung adenocarcinoma orient co-expression links?

A tumour that gains copies of gene A makes more of A. If A drives B, B's expression follows
A's copy number, while A's expression does not follow B's. For genes A and B on different
chromosomes the direction score is |corr(CN_A, E_B)| - |corr(CN_B, E_A)|; positive means
A -> B. Both genes must have copy number that moves their own expression (`CIS`).

- curated: how often the score points from TF to target in curated pairs, against the same
  TFs paired with random genes that are not their targets;
- regional: the same with each gene's copy number replaced by its chromosome arm's mean,
  so a direction that only reflects a whole amplified region scores the same; genes stay
  eligible by their own copy number;
- tfs: does the share of a gene's network links that point away from it find TFs?
"""

import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from .diagnostics import chromosome
from .evaluation import read_reference
from .public import DATASETS, read_expression, read_samples, rows, zscores
from .references import read_edges, read_gene_table
from .structure import ARM, location_codes

CIS = 0.3
MIN_LINKS = 10
CONTROLS = 10


def arm_means(copy_number: np.ndarray, arms: np.ndarray) -> np.ndarray:
    """Each gene's copy number replaced by the mean over genes on its arm, re-standardised."""
    regional = copy_number.copy()
    for arm in np.unique(arms[arms >= 0]):
        members = arms == arm
        regional[:, members] = copy_number[:, members].mean(axis=1, keepdims=True)
    return zscores(regional.T)


def cross(copy_number: np.ndarray, expression: np.ndarray) -> np.ndarray:
    """corr(CN_a, E_b) for every pair, as genes x genes."""
    return (copy_number.T @ expression / len(copy_number)).astype(np.float32)


def direction(effects: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.abs(effects[a, b]) - np.abs(effects[b, a])


def share_ci(values: np.ndarray, runs: int = 1000, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    means = [values[rng.integers(0, len(values), len(values))].mean() for _ in range(runs)]
    low, high = np.quantile(means, [0.025, 0.975])
    return {"share": float(values.mean()), "ci95": [float(low), float(high)]}


def curated(
    effects: np.ndarray,
    pairs: np.ndarray,
    eligible: np.ndarray,
    chromosomes: np.ndarray,
    linked: set[tuple[int, int]],
    seed: int = 0,
) -> dict:
    """Share of curated pairs oriented TF -> target, and of the same TFs with non-targets."""
    rng = np.random.default_rng(seed)
    a, b = pairs.T
    keep = eligible[a] & eligible[b] & (chromosomes[a] != chromosomes[b])
    a, b = a[keep], b[keep]
    targets = set(zip(a.tolist(), b.tolist(), strict=True))
    pool = np.flatnonzero(eligible)
    control_a, control_b = [], []
    for source in a:
        drawn = rng.choice(pool, 4 * CONTROLS)
        drawn = [
            g for g in drawn if chromosomes[g] != chromosomes[source] and (source, g) not in targets
        ][:CONTROLS]
        control_a += [source] * len(drawn)
        control_b += drawn
    oriented = direction(effects, a, b) > 0
    control = direction(effects, np.array(control_a, int), np.array(control_b, int)) > 0
    on_links = np.array([(i, j) in linked or (j, i) in linked for i, j in zip(a, b, strict=True)])
    return {
        "pairs": int(len(a)),
        "tf_to_target": share_ci(oriented) if len(a) else None,
        "tf_to_non_target": share_ci(control) if len(a) else None,
        "linked_pairs": int(on_links.sum()),
        "linked_tf_to_target": share_ci(oriented[on_links]) if on_links.any() else None,
    }


def random_auc_range(labels: np.ndarray, scores: np.ndarray, runs: int = 1000) -> list[float]:
    rng = np.random.default_rng(0)
    null = [roc_auc_score(rng.permutation(labels), scores) for _ in range(runs)]
    return [float(q) for q in np.quantile(null, [0.025, 0.975])]


def outward(
    effects: np.ndarray, links: np.ndarray, eligible: np.ndarray, chromosomes: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Per gene, its scored links pointing away from it and all its scored links."""
    a, b = links.T
    keep = eligible[a] & eligible[b] & (chromosomes[a] != chromosomes[b])
    a, b = a[keep], b[keep]
    outward = direction(effects, a, b) > 0
    size = len(eligible)
    total = np.bincount(a, minlength=size) + np.bincount(b, minlength=size)
    away = np.bincount(a, weights=outward, minlength=size) + np.bincount(
        b, weights=~outward, minlength=size
    )
    return away, total


def tf_scores(
    effects: np.ndarray,
    links: np.ndarray,
    eligible: np.ndarray,
    chromosomes: np.ndarray,
    tfs: np.ndarray,
) -> dict:
    """Share of each gene's network links pointing away from it, against the TF list."""
    away, total = outward(effects, links, eligible, chromosomes)
    scored = total >= MIN_LINKS
    share = away[scored] / total[scored]
    labels = tfs[scored]
    return {
        "genes": int(scored.sum()),
        "tfs": int(labels.sum()),
        "auc": float(roc_auc_score(labels, share)),
        "degree_auc": float(roc_auc_score(labels, total[scored])),
        "random_95": random_auc_range(labels, share),
    }


def read_luad(references: Path, genes: list[str], symbols: dict[str, str]):
    dataset = next(d for d in DATASETS if d.copy_number)
    expression, copy_number = read_samples(dataset, references, symbols)
    return zscores(rows(expression, genes)), zscores(rows(copy_number, genes)), dataset.name


def link_pairs(path: Path, index: dict[str, int]) -> np.ndarray:
    pairs = set()
    with path.open() as stream:
        for line in stream:
            source, target, _ = line.split("\t")
            i, j = index[source], index[target]
            pairs.add((min(i, j), max(i, j)))
    return np.array(sorted(pairs))


def run(root: Path, gene_table: Path, references: Path, data: Path, output: Path) -> dict:
    """`root` is the output of `public.build`; `data` holds `edges/` and the TF lists."""
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    genes, _ = read_expression(root / "expression.tsv")
    index = {g: i for i, g in enumerate(genes)}
    chromosomes = np.array([chromosome(table.get(g, {}).get("location", "")) for g in genes])
    arms = location_codes(tuple(genes), table, ARM)
    expression, copy_number, name = read_luad(references, genes, symbols)
    links = link_pairs(root / name / "graph.txt", index)
    linked = set(map(tuple, links.tolist()))
    curated_pairs = {
        source: np.array(
            [
                (index[s], index[t])
                for s, t, _ in read_edges(data / f"edges/{source}.tsv")
                if s in index and t in index and s != t
            ]
        )
        for source in ("collectri", "dorothea-abc", "trrust")
    }
    tf_lists = {
        list_name: np.array([g in read_reference(data / path) for g in genes])
        for list_name, path in (("lambert-2018", "lambert-2018-tfs.txt"), ("tfs", "tfs.txt"))
    }
    eligible = np.diag(cross(copy_number, expression)) >= CIS
    report: dict = {
        "dataset": name,
        "samples": int(len(expression)),
        "cis_threshold": CIS,
        "eligible_genes": int(eligible.sum()),
    }
    for level, cn in (("gene", copy_number), ("regional", arm_means(copy_number, arms))):
        effects = cross(cn, expression)
        report[level] = {
            "curated": {
                source: curated(effects, pairs, eligible, chromosomes, linked)
                for source, pairs in curated_pairs.items()
            },
            "tfs": {
                list_name: tf_scores(effects, links, eligible, chromosomes, labels)
                for list_name, labels in tf_lists.items()
            },
        }
        print(level, json.dumps(report[level]["curated"]["collectri"]["tf_to_target"]))
    output.mkdir(parents=True, exist_ok=True)
    (output / "copy_number.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
