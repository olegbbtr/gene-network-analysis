"""Do a gene's co-expression partners respond to its knockdown, and to it in particular?

`perturbation` counts affected pairs over a whole screen, so a few knockdowns that change
thousands of genes, or a growth programme that falls under any stress, could produce its
result. Here the knockdown is the unit.

Within each knockdown the measured genes are ranked, once by how far they fall and once by
how strongly they change in either direction, so every knockdown has the same distribution
whatever its breadth. A gene's partners are the measured genes on other chromosomes linked
to it positively in at least MAJORITY networks (also: in all networks; and, for negative
links, in at least NEGATIVE networks). For every knockdown with at least MIN_PARTNERS
partners:

- response: the mean rank of its partners under its own knockdown (0.5 = like any gene);
- specificity: the share of the other knockdowns in the screen under which the same partner
  genes rank lower than under its own. 0.5 means the partners respond to this knockdown no
  more than to any other; the share of knockdowns at or above 1 - EXTREME is compared with
  EXTREME. Knockdowns of a gene's own partners stay in the comparison, which works against
  the result.

Strata: all knockdowns; without cell-essential targets; without the housekeeping complexes
(`robustness.core_machinery`) among targets and partners; knockdowns without a broad effect;
knockdowns of known TFs. All thresholds and strata were fixed before the first run.
"""

import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, wilcoxon

from .applications import load
from .consensus import link_counts
from .gene_sets import read_gmt
from .perturbation import BROAD, SCREENS, read_screen
from .prediction import read_gene_list
from .robustness import core_machinery

MAJORITY = 8
NEGATIVE = 4
MIN_PARTNERS = 10
EXTREME = 0.05
LISTED = 25


def within_ranks(values: np.ndarray) -> np.ndarray:
    """Per row, each value's rank as a share of the row (ties averaged), between 0 and 1."""
    return ((rankdata(values, axis=1) - 0.5) / values.shape[1]).astype(np.float32)


def partner_means(ranks: np.ndarray, partners: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """`means[k, a]`: mean rank under knockdown k of the partners of knockdown a (NaN when a
    has fewer than MIN_PARTNERS partners); and the number of partners of each knockdown.
    Both inputs are knockdowns x measured genes."""
    sizes = partners.sum(axis=1)
    sums = ranks @ partners.T.astype(np.float32)
    enough = (sizes >= MIN_PARTNERS)[None, :]
    means = np.divide(sums, sizes[None, :], out=np.full(sums.shape, np.nan), where=enough)
    return means, sizes


def specificity(means: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Each knockdown's own partner mean, and the share of other knockdowns giving its
    partners a lower mean (ties count half)."""
    own = np.diag(means).copy()
    below = (means < own[None, :]).sum(axis=0)
    ties = (means == own[None, :]).sum(axis=0) - 1
    return own, (below + 0.5 * ties) / (len(means) - 1)


def summarise(own: np.ndarray, share: np.ndarray, chosen: np.ndarray) -> dict:
    chosen = chosen & np.isfinite(own)
    if chosen.sum() < MIN_PARTNERS:
        return {"knockdowns": int(chosen.sum())}
    return {
        "knockdowns": int(chosen.sum()),
        "median_partner_rank": float(np.median(own[chosen])),
        "partners_respond": float((own[chosen] > 0.5).mean()),
        "p": float(wilcoxon(own[chosen] - 0.5).pvalue),
        "median_specificity": float(np.median(share[chosen])),
        "most_specific": float((share[chosen] >= 1 - EXTREME).mean()),
        "most_specific_by_chance": EXTREME,
    }


def run(root: Path, gene_table: Path, references: Path, data: Path, output: Path) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    index = {gene: i for i, gene in enumerate(genes)}
    names = context.plain
    counts = link_counts([context.adjacency[n] for n in names])
    essential = read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols)
    essential |= read_gene_list(
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", context.symbols
    )
    essential = np.array([g in essential for g in genes])
    core = core_machinery(
        context, read_gmt(references / "reactome/ReactomePathways.gmt.zip", context.symbols)
    )
    listed = set((data / "lambert-2018-tfs.txt").read_text().split())
    tf = np.array([g in listed for g in genes])
    definitions = {
        "positive_majority": ("positive", MAJORITY),
        "positive_all": ("positive", len(names)),
        "negative": ("negative", NEGATIVE),
    }
    report: dict = {
        "majority": MAJORITY,
        "negative_networks": NEGATIVE,
        "min_partners": MIN_PARTNERS,
        "screens": {},
    }
    for screen_name, file in SCREENS.items():
        screen = read_screen(references / "perturb-seq-replogle-2022" / file)
        rows = np.array([r for r, g in enumerate(screen.targets) if g in index])
        target = np.array([index[screen.targets[r]] for r in rows])
        columns = np.array([c for c, g in enumerate(screen.measured) if g in index])
        measured = np.array([index[screen.measured[c]] for c in columns])
        effects = screen.effects[np.ix_(rows, columns)]
        ranks = {"fall": within_ranks(-effects), "change": within_ranks(np.abs(effects))}
        apart = context.groups[target][:, None] != context.groups[measured][None, :]
        changed = screen.changed[rows]
        everyone = np.ones(len(target), dtype=bool)
        strata = {
            "all": (everyone, None),
            "without_cell_essential": (~essential[target], None),
            "without_housekeeping": (~core[target], ~core[measured]),
            "not_broad": (changed < np.quantile(changed, 1 - BROAD), None),
            "tfs": (tf[target], None),
        }
        entry: dict = {"knockdowns": int(len(target)), "measured_genes": int(len(measured))}
        for name, (kind, least) in definitions.items():
            linked = (counts[kind][np.ix_(target, measured)] >= least) & apart
            entry[name] = {}
            for stratum, (chosen, kept_partners) in strata.items():
                partners = linked if kept_partners is None else linked & kept_partners[None, :]
                result: dict = {}
                for measure, values in ranks.items():
                    means, sizes = partner_means(values, partners)
                    own, share = specificity(means)
                    result[measure] = summarise(own, share, chosen)
                    if (name, stratum, measure) == ("positive_majority", "all", "fall"):
                        order = np.argsort(-np.where(np.isfinite(own), share, -1), kind="stable")
                        entry["most_specific_knockdowns"] = [
                            {
                                "symbol": context.table[genes[target[i]]]["symbol"],
                                "partners": int(sizes[i]),
                                "partner_rank": float(own[i]),
                                "specificity": float(share[i]),
                                "tf": bool(tf[target[i]]),
                                "cell_essential": bool(essential[target[i]]),
                            }
                            for i in order[:LISTED]
                        ]
                entry[name][stratum] = result
            print(screen_name, name, json.dumps(entry[name]), flush=True)
        report["screens"][screen_name] = entry
    output.mkdir(parents=True, exist_ok=True)
    (output / "partners.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
