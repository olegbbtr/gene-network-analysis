"""Do network hubs and links match what happens when a gene is knocked down?

A TF list says what a gene is; an intervention shows what it does. Replogle et al. (2022)
knocked down every expressed gene of K562 cells by CRISPRi, and the essential genes of K562
and RPE1 cells, and sequenced the cells carrying each knockdown (Perturb-seq). Their
pseudobulk files give, per knockdown, the z-score of every measured gene against control
cells and the number of differentially expressed genes (Anderson-Darling test).

Gene level. A knockdown has a broad effect when its number of differentially expressed
genes is in the top BROAD share of the screen. The models of `robustness` (held-out
chromosomes) ask whether a gene's place in the tissue networks predicts that, beyond its
expression in the tissue and three properties of the experiment: how far the knockdown
lowered its target, how many cells carried it and the target's expression in the cell line.
The comparison is repeated without cell-essential genes, whose loss disturbs any cell, and
the odds of a broad effect are given for known TFs against other genes.

Pair level. A knockdown of A affects B when |z| is in the top AFFECTED share of all pairs on
different chromosomes (CRISPRi can also silence neighbours of its target). Pairs are grouped
by the number of networks linking them, as in `consensus`: are linked pairs affected more
often, and does B fall when a positively linked A is knocked down?

Direction. The networks keep an edge A->B without B->A for some pairs. Where both knockdowns
were measured and at least one affects the other gene, does the kept direction agree with the
stronger effect?

One knockdown per gene is used: the one that lowers its target most, then the one in most
cells; knockdowns in fewer than MIN_CELLS cells are left out. The thresholds, covariates and
tests were fixed before the first run.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_array
from scipy.stats import fisher_exact

from .applications import load
from .consensus import link_counts
from .graph import Graph
from .prediction import network_features, read_gene_list
from .robustness import blocks, compare, descriptor_aucs, held_out
from .structure import ARM, location_codes

BROAD = 0.1
AFFECTED = 0.01
MIN_CELLS = 25
SCREENS = {
    "k562-genome-wide": "K562_gwps_normalized_bulk_01.h5ad",
    "k562-essential": "K562_essential_normalized_bulk_01.h5ad",
    "rpe1-essential": "rpe1_normalized_bulk_01.h5ad",
}


@dataclass(frozen=True)
class Screen:
    targets: list[str]  # Ensembl ID of the gene knocked down, one row each
    measured: list[str]  # Ensembl IDs of the measured genes
    effects: np.ndarray  # targets x measured: z-scores against control cells
    changed: np.ndarray  # differentially expressed genes per knockdown
    covariates: np.ndarray  # targets x 3: change of the target, log cells, control expression


def one_per_gene(names: list[str], knockdown: np.ndarray, cells: np.ndarray) -> dict[str, int]:
    """The row kept for each target gene, named by the Ensembl ID that ends a row's name."""
    best: dict[str, tuple[tuple[float, float], int]] = {}
    for row, name in enumerate(names):
        gene = name.rsplit("_", 1)[-1]
        if not gene.startswith("ENSG") or cells[row] < MIN_CELLS:
            continue
        key = (knockdown[row] if np.isfinite(knockdown[row]) else np.inf, -cells[row])
        if gene not in best or key < best[gene][0]:
            best[gene] = (key, row)
    return {gene: row for gene, (_, row) in best.items()}


def filled(values: np.ndarray) -> np.ndarray:
    """Missing values replaced by the median of the others."""
    return np.where(np.isfinite(values), values, np.nanmedian(values))


def read_screen(path: Path) -> Screen:
    with h5py.File(path) as store:
        names = [name.decode() for name in store["obs/gene_transcript"][:]]
        knockdown = store["obs/pct_expr"][:]
        cells = store["obs/num_cells_filtered"][:]
        control = store["obs/control_expr"][:]
        changed = store["obs/anderson_darling_counts"][:]
        measured = [gene.decode() for gene in store["var/gene_id"][:]]
        effects = store["X"][:]
    kept = one_per_gene(names, knockdown, cells)
    rows = np.array(list(kept.values()))
    covariates = np.column_stack([filled(knockdown), np.log(cells), filled(control)])
    return Screen(
        list(kept), measured, effects[rows], changed[rows].astype(float), covariates[rows]
    )


def tf_odds(broad: np.ndarray, tf: np.ndarray, keep: np.ndarray) -> dict:
    """Odds of a broad effect for TFs against the other kept genes (Fisher's exact test)."""
    table = [
        [int((broad & tf & keep).sum()), int((~broad & tf & keep).sum())],
        [int((broad & ~tf & keep).sum()), int((~broad & ~tf & keep).sum())],
    ]
    test = fisher_exact(table)
    return {
        "tfs": table[0][0] + table[0][1],
        "broad_tfs": table[0][0],
        "broad_others": table[1][0],
        "odds_ratio": float(test.statistic),
        "p": float(test.pvalue),
    }


def by_links(effects: np.ndarray, count: np.ndarray, valid: np.ndarray, size: int) -> dict:
    """Per number of linking networks: pairs, the share affected and, of those, the share
    where the measured gene falls. `effects`, `count` and `valid` are targets x measured."""
    magnitude = np.abs(effects)
    threshold = float(np.quantile(magnitude[valid], 1 - AFFECTED))
    affected = valid & (magnitude >= threshold)
    pairs = np.bincount(count[valid], minlength=size)
    hits = np.bincount(count[affected], minlength=size)
    falls = np.bincount(count[affected & (effects < 0)], minlength=size)
    return {
        "threshold": threshold,
        "falls_among_all_affected": float(falls.sum() / hits.sum()),
        "by_networks": [
            {
                "networks": k,
                "pairs": int(pairs[k]),
                "affected": float(hits[k] / pairs[k]) if pairs[k] else None,
                "falls_among_affected": float(falls[k] / hits[k]) if hits[k] else None,
            }
            for k in range(size)
        ],
    }


def direction(edges: np.ndarray, effects: np.ndarray, valid: np.ndarray, threshold: float) -> dict:
    """Among one-way edges i->j with an effect of at least `threshold` in either direction:
    how often is the knockdown of i the stronger one? All arrays are genes x genes, with
    `effects[i, j]` the effect of knocking down i on j."""
    magnitude = np.abs(effects)
    chosen = edges & ~edges.T & valid & (np.maximum(magnitude, magnitude.T) >= threshold)
    agree = (magnitude > magnitude.T)[chosen]
    if not len(agree):
        return {"pairs": 0}
    half = 1.96 * np.sqrt(0.25 / len(agree))
    return {
        "pairs": int(len(agree)),
        "agreement": float(agree.mean()),
        "chance_95": [0.5 - half, 0.5 + half],
    }


def run(root: Path, gene_table: Path, references: Path, data: Path, output: Path) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    index = {gene: i for i, gene in enumerate(genes)}
    arms = location_codes(tuple(genes), context.table, ARM)
    names = context.plain
    features = {
        n: network_features(Graph(tuple(genes), context.adjacency[n]), genes, arms) for n in names
    }
    counts = link_counts([context.adjacency[n] for n in names])
    essential = read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols)
    essential |= read_gene_list(
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", context.symbols
    )
    essential = np.array([g in essential for g in genes])
    listed = set((data / "lambert-2018-tfs.txt").read_text().split())
    tf = np.array([g in listed for g in genes])
    sets = {name: [name] for name in names} | {"all": names}
    report: dict = {"broad_share": BROAD, "affected_share": AFFECTED, "screens": {}}
    for screen_name, file in SCREENS.items():
        screen = read_screen(references / "perturb-seq-replogle-2022" / file)
        rows = np.array([r for r, g in enumerate(screen.targets) if g in index])
        target = np.array([index[screen.targets[r]] for r in rows])
        columns = np.array([c for c, g in enumerate(screen.measured) if g in index])
        measured = np.array([index[screen.measured[c]] for c in columns])
        perturbed = np.zeros(len(genes), dtype=bool)
        perturbed[target] = True
        changed = np.zeros(len(genes))
        changed[target] = screen.changed[rows]
        covariates = np.zeros((len(genes), screen.covariates.shape[1]))
        covariates[target] = screen.covariates[rows]
        broad = perturbed & (changed >= np.quantile(changed[perturbed], 1 - BROAD))
        kept = {"all_knockdowns": perturbed, "without_cell_essential": perturbed & ~essential}
        entry: dict = {
            "knockdowns": int(perturbed.sum()),
            "measured_genes": int(len(measured)),
            "broad": int(broad.sum()),
            "broad_from_changed_genes": float(changed[broad].min()),
            "broad_among_cell_essential": float(broad[perturbed & essential].mean()),
            "broad_among_others": float(broad[perturbed & ~essential].mean()),
            "tf_odds": {name: tf_odds(broad, tf, keep) for name, keep in kept.items()},
            "descriptor_auc": {
                name: {n: descriptor_aucs(features[n], broad, keep) for n in names}
                for name, keep in kept.items()
            },
            "models": {name: {} for name in kept},
        }
        for set_name, members in sets.items():
            parts = blocks(context, features, members)
            parts["expression"] = np.column_stack([parts["expression"], covariates])
            links = [features[n].links for n in members]
            for name, keep in kept.items():
                predictions = held_out(parts, links, broad, context.groups, keep)
                entry["models"][name][set_name] = compare(broad[keep], predictions, runs=500)
            print(
                screen_name,
                set_name,
                json.dumps({name: entry["models"][name][set_name] for name in kept}),
                flush=True,
            )
        effects = screen.effects[np.ix_(rows, columns)]
        apart = context.groups[target][:, None] != context.groups[measured][None, :]
        entry["pairs"] = {
            kind: by_links(effects, counts[kind][np.ix_(target, measured)], apart, len(names) + 1)
            for kind in ("positive", "negative")
        }
        both = np.intersect1d(target, measured)
        row_of = {gene: r for r, gene in zip(rows, target, strict=True)}
        column_of = {gene: c for c, gene in zip(columns, measured, strict=True)}
        square = screen.effects[np.ix_([row_of[g] for g in both], [column_of[g] for g in both])]
        valid = context.groups[both][:, None] != context.groups[both][None, :]
        threshold = entry["pairs"]["positive"]["threshold"]
        entry["direction"] = {
            n: direction(
                csr_array(context.adjacency[n])[both][:, both].toarray() != 0,
                square,
                valid,
                threshold,
            )
            for n in names
        }
        report["screens"][screen_name] = entry
        print(
            screen_name,
            json.dumps({k: v for k, v in entry.items() if k not in ("models", "descriptor_auc")}),
            flush=True,
        )
    output.mkdir(parents=True, exist_ok=True)
    (output / "perturbation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
