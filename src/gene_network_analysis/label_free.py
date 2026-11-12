"""Essential-gene models whose features carry no essentiality information.

CEGv2 labels are only the training target on the training chromosomes and the answer key
on the held-out ones; no feature is computed from them. Four groups:

- gene: expression mean and SD per dataset;
- local: a gene's own links (degree, signs, weights, flow, clustering),
  plus out- and in-link counts above the 25th, 50th and 75th percentiles of each network's
  own |z| (`cutoff_counts`); the training-essential neighbour share is left out;
- global: forward/reverse PageRank scores and positions, k-core depth, eigenvector centrality;
- biological: what is known about the gene's partners. Whether its neighbours share a
  Reactome pathway of PAIR_SIZES genes with it (as in `consensus`): the share of its
  neighbours that do, per network; how many of its partners linked positively in every
  network do; whether it has a pathway at all. And the share of its links that stay on its
  own chromosome arm, per network.

Five models are fitted with the nested penalty selection of `feature_audit`: each group
alone and all four together. Impact is leave-one-out: the held-out AUC lost when a block is
removed from all four, for gene, graph (local and global removed together, because they
carry largely the same information and each covers for the other) and biological. That needs
three more fits. Intervals bootstrap genes over the fixed held-out predictions, on CEGv2 and
on DepMap essentials outside CEGv2.
"""

import json
from itertools import combinations
from math import factorial
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed
from scipy.sparse import csr_array
from sklearn.metrics import roc_auc_score

from . import cutoff_counts
from .applications import load
from .consensus import together
from .feature_audit import PENALTIES, metrics
from .feature_comparison import GLOBAL, feature_blocks
from .gene_sets import membership, read_gmt, restrict
from .prediction import SHAPE, nested_cross_validated, read_gene_list
from .references import read_gene_table
from .report import revision, sha256
from .robustness import top_indices

REACTOME = "reactome/ReactomePathways.gmt.zip"
PAIR_SIZES = (2, 200)
GROUPS = {
    "gene": ("expression",),
    "local": ("degree", "shape", "cutoff_counts"),
    "global": GLOBAL,
    "biological": ("pathway_coherence", "arm_share"),
}
LAYOUT = "biological group: pathway coherence + chromosome-arm share"
SUBSETS = [c for size in range(1, 5) for c in combinations(GROUPS, size)]
COMBINATIONS = {"+".join(c): tuple(p for group in c for p in GROUPS[group]) for c in SUBSETS}
ALL = "+".join(GROUPS)
LABELS = {
    "gene": "Gene only",
    "local": "Local graph only",
    "global": "Global graph only",
    "biological": "Biological only",
    ALL: "All four",
}
BLOCKS = {"gene": ("gene",), "graph": ("local", "global"), "biological": ("biological",)}
WITHOUT = {
    block: "+".join(g for g in GROUPS if g not in removed) for block, removed in BLOCKS.items()
}
MODELS = {k: COMBINATIONS[k] for k in [*LABELS, *WITHOUT.values()]}
CHANCE = 0.5
BOOTSTRAP = 1000


def shapley(auc: dict[str, float]) -> dict[str, float]:
    """Each group's AUC gain averaged over every order of adding the groups to chance.

    `auc` maps "+"-joined subsets in GROUPS order to AUC; the empty model scores CHANCE.
    The values sum to auc[ALL] - CHANCE.
    """
    groups = list(GROUPS)
    n = len(groups)

    def score(members) -> float:
        return auc["+".join(g for g in groups if g in members)] if members else CHANCE

    values = {}
    for group in groups:
        others = [g for g in groups if g != group]
        total = 0.0
        for size in range(n):
            weight = factorial(size) * factorial(n - size - 1) / factorial(n)
            for chosen in combinations(others, size):
                total += weight * (score({*chosen, group}) - score(set(chosen)))
        values[group] = total
    return values


def leave_one_out(auc: dict[str, float]) -> dict[str, float]:
    """Each group's AUC loss when it alone is removed from all four."""
    return {group: auc[ALL] - auc["+".join(g for g in GROUPS if g != group)] for group in GROUPS}


IMPACTS = {
    "shapley": (shapley, "AUC added, averaged over every order of adding groups to chance"),
    "leave_one_out": (leave_one_out, "AUC lost when the group alone is removed from all four"),
}


def impact(labels: np.ndarray, predictions: dict[str, np.ndarray], method, seed: int = 0) -> dict:
    """Per group, `method`'s AUC values and their shares of the four values' sum (sums to 1),
    with gene-bootstrap 95% intervals over the fixed held-out predictions."""

    def measure(chosen: np.ndarray) -> tuple[dict, dict]:
        values = method(
            {k: roc_auc_score(labels[chosen], p[chosen]) for k, p in predictions.items()}
        )
        total = sum(values.values())
        return values, {g: v / total for g, v in values.items()}

    values, shares = measure(np.arange(len(labels)))
    rng = np.random.default_rng(seed)
    samples = []
    while len(samples) < BOOTSTRAP:
        chosen = rng.integers(0, len(labels), len(labels))
        if labels[chosen].any() and not labels[chosen].all():
            samples.append(measure(chosen))
    return {
        group: {
            "auc": values[group],
            "auc_ci95": np.quantile([v[group] for v, _ in samples], [0.025, 0.975]).tolist(),
            "share": shares[group],
            "share_ci95": np.quantile([s[group] for _, s in samples], [0.025, 0.975]).tolist(),
        }
        for group in GROUPS
    }


def removal(labels: np.ndarray, predictions: dict[str, np.ndarray], seed: int = 0) -> dict:
    """Per block of BLOCKS, the AUC lost when it is removed from all four, with a
    gene-bootstrap 95% interval over the fixed held-out predictions."""

    def losses(chosen: np.ndarray) -> dict:
        full = roc_auc_score(labels[chosen], predictions[ALL][chosen])
        return {
            b: full - roc_auc_score(labels[chosen], predictions[m][chosen])
            for b, m in WITHOUT.items()
        }

    values = losses(np.arange(len(labels)))
    rng = np.random.default_rng(seed)
    samples = []
    while len(samples) < BOOTSTRAP:
        chosen = rng.integers(0, len(labels), len(labels))
        if labels[chosen].any() and not labels[chosen].all():
            samples.append(losses(chosen))
    return {
        b: {
            "removed": BLOCKS[b],
            "auc_loss": values[b],
            "auc_loss_ci95": np.quantile([v[b] for v in samples], [0.025, 0.975]).tolist(),
        }
        for b in BLOCKS
    }


def removals(known: np.ndarray, target: np.ndarray, tuned: dict[str, np.ndarray]) -> dict:
    new = ~known
    return {
        "method": "AUC lost when the block is removed from all four",
        "cegv2": removal(known, tuned),
        "depmap": removal(target[new], {k: p[new] for k, p in tuned.items()}),
    }


def impacts(known: np.ndarray, target: np.ndarray, tuned: dict[str, np.ndarray]) -> dict:
    new = ~known
    return {
        name: {
            "method": text,
            "cegv2": impact(known, tuned, method),
            "depmap": impact(target[new], {k: p[new] for k, p in tuned.items()}, method),
        }
        for name, (method, text) in IMPACTS.items()
    }


def essential_labels(references: Path, genes: list[str], symbols: dict[str, str]):
    known = np.isin(genes, list(read_gene_list(references / "gene-lists/CEGv2.txt", symbols)))
    target = np.isin(
        genes,
        list(
            read_gene_list(references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", symbols)
        ),
    )
    return known, target


def top_candidates(genes, table, scores, known, target) -> dict:
    """The combined model's top predictions outside CEGv2, in the `candidates` input shape."""
    new = np.flatnonzero(~known)
    return {
        ALL: [
            {
                "gene": genes[i],
                "symbol": table[genes[i]]["symbol"],
                "depmap": bool(target[i]),
                "score": float(scores[i]),
            }
            for i in new[top_indices(scores[new])]
        ]
    }


def rescore(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
    """Recompute the impact section from the saved held-out predictions, without refitting."""
    saved = np.load(root / "label-free-predictions.npz")
    genes = saved["genes"].tolist()
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    known, target = essential_labels(references, genes, symbols)
    path = output / "label-free-groups.json"
    report = json.loads(path.read_text())
    if all(k in saved for k in COMBINATIONS):
        report["impact"] = impacts(known, target, {k: saved[k] for k in COMBINATIONS})
    report["discovery"] = {
        "top_candidates": top_candidates(genes, table, saved[ALL], known, target)
    }
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["impact"]))
    return report


def consistent_partners(adjacency: list[csr_array]) -> csr_array:
    """Pairs linked positively (summed weight in both directions) in every network."""
    count = csr_array(adjacency[0].shape, dtype=np.int16)
    for a in adjacency:
        count = count + (a + a.T > 0).astype(np.int16)
    consistent = (count >= len(adjacency)).astype(np.float64).tocsr()
    consistent.setdiag(0)
    consistent.eliminate_zeros()
    return consistent


def coherence(
    links: list[csr_array], adjacency: list[csr_array], shared: np.ndarray, annotated: np.ndarray
) -> np.ndarray:
    columns = [
        np.asarray(m.multiply(shared).sum(axis=1)).ravel()
        / np.maximum(np.asarray(m.sum(axis=1)).ravel(), 1)
        for m in links
    ]
    partners = consistent_partners(adjacency)
    related = np.asarray(partners.multiply(shared).sum(axis=1)).ravel()
    return np.column_stack([*columns, np.log1p(related), annotated])


def build(root: Path, gene_table: Path, references: Path, inputs: dict[str, str]):
    """Feature blocks and chromosome groups, cached in `root` under the input hashes."""
    cache = root / "label-free-features.npz"
    if cache.exists():
        saved = np.load(cache)
        if json.loads(str(saved["inputs"])) == inputs:
            parts = {p: saved[p] for group in GROUPS.values() for p in group}
            return saved["genes"].tolist(), saved["groups"], parts
    context = load(root, gene_table)
    genes = context.genes
    parts, links = feature_blocks(context)
    arm = [i for i in range(parts["shape"].shape[1]) if i % len(SHAPE) == SHAPE.index("arm_share")]
    parts["arm_share"] = parts["shape"][:, arm]
    parts["shape"] = np.delete(parts["shape"], arm, axis=1)
    parts["cutoff_counts"] = cutoff_counts.load(root, genes, context.plain)
    reactome = read_gmt(references / REACTOME, context.symbols)
    parts["pathway_coherence"] = coherence(
        links,
        [context.adjacency[n] for n in context.plain],
        together(genes, reactome),
        membership(genes, restrict(reactome, genes, *PAIR_SIZES)).any(axis=1).astype(float),
    )
    parts = {p: parts[p] for group in GROUPS.values() for p in group}
    np.savez(
        cache, inputs=json.dumps(inputs), genes=np.array(genes), groups=context.groups, **parts
    )
    return genes, context.groups, parts


def run(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
    names = [n for n in json.loads((root / "stats.json").read_text()) if n != "settings"]
    paths = [
        gene_table,
        root / "expression.tsv",
        root / "stats.json",
        references / "gene-lists/CEGv2.txt",
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv",
        references / REACTOME,
        root / cutoff_counts.CACHE,
        *(root / n / "graph.txt" for n in names if not n.endswith("-copy-number-removed")),
    ]
    inputs = {str(path): sha256(path) for path in paths} | {"layout": LAYOUT}
    genes, groups, parts = build(root, gene_table, references, inputs)
    print("features", flush=True)
    table = read_gene_table(gene_table)
    known, target = essential_labels(
        references, genes, {row["symbol"]: gene for gene, row in table.items()}
    )
    fitted = Parallel(n_jobs=min(len(MODELS), 12), verbose=10)(
        delayed(nested_cross_validated)(
            {p: parts[p] for p in model}, model, known, groups, [], penalties=PENALTIES
        )
        for model in MODELS.values()
    )
    tuned = {name: result[0] for name, result in zip(MODELS, fitted, strict=True)}
    choices = {name: result[1] for name, result in zip(MODELS, fitted, strict=True)}
    np.savez_compressed(root / "label-free-predictions.npz", genes=np.array(genes), **tuned)
    new = ~known
    results = metrics(known, target, tuned)
    report = {
        "implementation": revision(),
        "inputs": inputs,
        "settings": {
            "models": MODELS,
            "columns": {
                name: sum(parts[p].shape[1] for p in model) for name, model in MODELS.items()
            },
            "pair_sizes": PAIR_SIZES,
            "cutoff_percentiles": cutoff_counts.PERCENTILES,
            "cutoffs": cutoff_counts.cutoffs(root),
            "outer_folds": 5,
            "inner_folds": 4,
            "C_grid": PENALTIES,
            "labels": "CEGv2 only as training target and held-out answer key; no label features",
        },
        "discovery_base_rate": float(target[new].mean()),
        "main_comparison": [{"key": k, "label": LABELS[k]} | results[k] for k in LABELS],
        "models": results,
        "selected_C": choices,
        "leave_one_out": removals(known, target, tuned),
        "discovery": {"top_candidates": top_candidates(genes, table, tuned[ALL], known, target)},
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "label-free-groups.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["main_comparison"]))
    return report
