"""Are links seen in more tissues more often known relationships between genes?

Every pair of genes is counted by the number of networks (0-15) that link it in either
direction, and separately by the number that link it positively and negatively. For each
count, the share of pairs with a known relationship is compared with the share among all
pairs:

- reactome: both genes in one Reactome pathway of at most LARGEST genes;
- family: both in one HGNC gene group of at most LARGEST genes;
- string: a STRING v12 link (combined score >= 700; STRING is partly built from
  co-expression, so this one is not independent);
- regulation: a curated TF-target pair (CollecTRI), in either order.

Tumour links often join neighbours on a chromosome (copy number), so the same table is
given for pairs on different chromosomes, and for those of them outside the large
housekeeping complexes (`robustness.core_machinery`), which co-vary in every tissue.
"""

import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from .applications import load
from .baselines import read_string
from .gene_sets import hgnc_groups, membership, read_gmt, restrict
from .references import GeneMap, read_edges
from .robustness import core_machinery

LARGEST = 200


def together(genes: list[str], sets: dict[str, set[str]]) -> np.ndarray:
    """Genes x genes: the two genes share a set of 2 to LARGEST of these genes."""
    chosen = restrict(sets, genes, smallest=2, largest=LARGEST)
    members = csr_array(membership(genes, chosen).astype(np.float32))
    return (members @ members.T).toarray() > 0


def listed_pairs(genes: list[str], pairs: list[tuple[str, str]]) -> np.ndarray:
    index = {g: i for i, g in enumerate(genes)}
    matrix = np.zeros((len(genes), len(genes)), dtype=bool)
    for first, second in pairs:
        if first in index and second in index:
            matrix[index[first], index[second]] = matrix[index[second], index[first]] = True
    return matrix


def link_counts(adjacencies: list[csr_array]) -> dict[str, np.ndarray]:
    """Genes x genes: in how many networks a pair is linked at all, positively, negatively.

    A pair linked both ways takes the sign of the summed weights.
    """
    size = adjacencies[0].shape[0]
    counts = {
        name: np.zeros((size, size), dtype=np.int8) for name in ("any", "positive", "negative")
    }
    for adjacency in adjacencies:
        both = (adjacency + adjacency.T).tocoo()
        for name, chosen in (
            ("any", both.data != 0),
            ("positive", both.data > 0),
            ("negative", both.data < 0),
        ):
            counts[name][both.row[chosen], both.col[chosen]] += 1
    return counts


def tally(count: np.ndarray, known: dict[str, np.ndarray], pairs: np.ndarray, size: int) -> dict:
    """Per number of networks, the pairs among `pairs` and the share of each known kind."""
    totals = np.bincount(count[pairs], minlength=size)
    hits = {name: np.bincount(count[pairs & m], minlength=size) for name, m in known.items()}
    return {
        "base": {"pairs": int(totals.sum())}
        | {name: float(h.sum() / totals.sum()) for name, h in hits.items()},
        "by_networks": [
            {"networks": k, "pairs": int(totals[k])}
            | {name: float(h[k] / totals[k]) if totals[k] else None for name, h in hits.items()}
            for k in range(size)
        ],
    }


def run(root: Path, gene_table: Path, references: Path, data: Path, output: Path) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    names = context.plain
    reactome = read_gmt(references / "reactome/ReactomePathways.gmt.zip", context.symbols)
    known = {
        "reactome": together(genes, reactome),
        "family": together(genes, hgnc_groups(context.table)),
        "string": read_string(
            references / "string-v12/9606.protein.info.v12.0.txt.gz",
            references / "string-v12/9606.protein.links.v12.0.txt.gz",
            genes,
            GeneMap(references / "hgnc/hgnc_complete_set.txt").by_symbol,
        ).toarray()
        > 0,
        "regulation": listed_pairs(
            genes, [(a, b) for a, b, _ in read_edges(data / "edges/collectri.tsv")]
        ),
    }
    size = len(genes)
    upper = np.triu(np.ones((size, size), dtype=bool), 1)
    apart = upper & (context.groups[:, None] != context.groups[None, :])
    core = core_machinery(context, reactome)
    pairs = {
        "all_pairs": upper,
        "different_chromosomes": apart,
        "different_chromosomes_without_core": apart & ~core[:, None] & ~core[None, :],
    }
    counts = link_counts([context.adjacency[n] for n in names])
    report: dict = {
        "networks": names,
        "largest_set": LARGEST,
        "links": {
            kind: {
                name: tally(count, known, chosen, len(names) + 1) for name, chosen in pairs.items()
            }
            for kind, count in counts.items()
        },
    }
    for state in ("healthy", "tumour"):
        members = [context.adjacency[n] for n in names if context.stats[n]["state"] == state]
        counts[state] = link_counts(members)["positive"]
    chosen = pairs["different_chromosomes_without_core"]
    healthy, tumour = counts["healthy"][chosen], counts["tumour"][chosen]
    shape = (int(healthy.max()) + 1, int(tumour.max()) + 1)
    cell = healthy.astype(np.int64) * shape[1] + tumour
    totals = np.bincount(cell, minlength=shape[0] * shape[1])
    hits = np.bincount(cell[known["reactome"][chosen]], minlength=shape[0] * shape[1])
    report["positive_links_healthy_by_tumour"] = {
        "pairs": totals.reshape(shape).tolist(),
        "reactome_pairs": hits.reshape(shape).tolist(),
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "consensus.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
