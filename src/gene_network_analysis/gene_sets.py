"""Gene sets and the neighbour-share engine that scores them on a network.

A gene's score for a set is the share of its network neighbours that are set members.
`held_out_scores` hides the members on held-out chromosomes, so a gene's own label never
enters its score; `set_aucs` turns the scores into one ROC AUC per set.

Readers turn public sources into gene sets of Ensembl IDs: HGNC gene groups, Reactome
pathways (GMT), DepMap lineage-selective dependencies and GWAS Catalog traits.
"""

import csv
import gzip
import io
import re
import zipfile
from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

GWAS_SIGNIFICANCE = 5e-8
DEPENDENCY = -0.5
SELECTIVITY = -0.25
MIN_MODELS = 5


def membership(genes: list[str], sets: dict[str, set[str]]) -> np.ndarray:
    """Genes x sets boolean matrix, in the order of `sets`."""
    index = {g: i for i, g in enumerate(genes)}
    matrix = np.zeros((len(genes), len(sets)), dtype=bool)
    for column, members in enumerate(sets.values()):
        matrix[[index[g] for g in members if g in index], column] = True
    return matrix


def restrict(
    sets: dict[str, set[str]], genes: list[str], smallest: int = 20, largest: int = 500
) -> dict[str, set[str]]:
    """Sets cut to `genes`, keeping those with `smallest` to `largest` members."""
    universe = set(genes)
    cut = {name: members & universe for name, members in sets.items()}
    return {name: m for name, m in cut.items() if smallest <= len(m) <= largest}


def neighbour_scores(links: csr_array, labels: np.ndarray) -> np.ndarray:
    degree = np.asarray(links.sum(axis=1)).ravel()
    return (links @ labels.astype(np.float64)) / np.maximum(degree, 1)[:, None]


def personalised_pagerank(
    links: csr_array, labels: np.ndarray, alpha: float = 0.85, iterations: int = 50
) -> np.ndarray:
    """Random walk with restart to each set's members, following edge direction.

    One column per set. Walkers at genes without outgoing links are lost.
    """
    degree = np.asarray(links.sum(axis=1)).ravel()
    walk = (links.multiply(1 / np.maximum(degree, 1)[:, None])).T.tocsr()
    restart = labels.astype(np.float64)
    restart /= np.maximum(restart.sum(axis=0), 1)
    scores = restart.copy()
    for _ in range(iterations):
        scores = alpha * (walk @ scores) + (1 - alpha) * restart
    return scores


def condensed_pagerank(
    links: csr_array, labels: np.ndarray, alpha: float = 0.85, iterations: int = 50
) -> np.ndarray:
    """Personalised PageRank on the DAG of strongly connected groups, one walk per set.

    Each group becomes one node, restarting in proportion to its members; a gene's score
    is its group's score divided by the group's size.
    """
    _, group = connected_components(links, directed=True, connection="strong")
    size = group.max() + 1
    collapse = csr_array((np.ones(len(group)), (np.arange(len(group)), group)))
    between = (collapse.T @ links @ collapse).tocsr()
    between.setdiag(0)
    between.eliminate_zeros()
    restart = collapse.T @ labels.astype(np.float64)
    scores = personalised_pagerank((between != 0).astype(np.float64), restart, alpha, iterations)
    members = np.bincount(group, minlength=size)
    return scores[group] / members[group][:, None]


def held_out_scores(
    links: csr_array,
    labels: np.ndarray,
    groups: np.ndarray,
    folds: int = 5,
    scorer: Callable[[csr_array, np.ndarray], np.ndarray] = neighbour_scores,
) -> np.ndarray:
    """Scores for every gene and set, from training-fold members only."""
    scores = np.zeros(labels.shape)
    for train, test in GroupKFold(folds).split(labels, groups=groups):
        known = np.zeros_like(labels)
        known[train] = labels[train]
        scores[test] = scorer(links, known)[test]
    return scores


def set_aucs(scores: np.ndarray, labels: np.ndarray) -> np.ndarray:
    return np.array([roc_auc_score(labels[:, j], scores[:, j]) for j in range(labels.shape[1])])


def hgnc_groups(gene_table: dict[str, dict]) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    for gene, row in gene_table.items():
        for group in row["gene_group"]:
            sets.setdefault(group, set()).add(gene)
    return sets


def read_gmt(path: Path, symbols: dict[str, str]) -> dict[str, set[str]]:
    """A GMT file (or a ZIP holding one): name, ID, then member symbols per line."""
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            text = archive.read(archive.namelist()[0]).decode()
    else:
        text = path.read_text()
    sets = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) > 2:
            sets[fields[0]] = {symbols[s] for s in fields[2:] if s in symbols}
    return sets


def _depmap_symbol(column: str) -> str:
    return column.split(" (", 1)[0]


def depmap_lineages(
    effects: Path, models: Path, symbols: dict[str, str], lineages: dict[str, str]
) -> tuple[list[str | None], dict[str, tuple[np.ndarray, np.ndarray]], dict[str, int]]:
    """Gene IDs (None if unmapped), each lineage's mean gene effect inside and outside it,
    and model counts. Lineages with fewer than MIN_MODELS models get no means.

    `lineages` maps DepMap OncotreeLineage to a name.
    """
    with models.open(newline="") as stream:
        lineage_of = {row["ModelID"]: row["OncotreeLineage"] for row in csv.DictReader(stream)}
    with effects.open(newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        ids, rows = [], []
        for row in reader:
            ids.append(row[0])
            rows.append([float(v) if v else np.nan for v in row[1:]])
    values = np.array(rows)
    gene_ids = [symbols.get(_depmap_symbol(c)) for c in header[1:]]
    labels = np.array([lineage_of.get(i, "") for i in ids])
    means, counts = {}, {}
    for lineage, name in lineages.items():
        inside = labels == lineage
        counts[name] = int(inside.sum())
        if inside.sum() >= MIN_MODELS:
            means[name] = (
                np.nanmean(values[inside], axis=0),
                np.nanmean(values[~inside], axis=0),
            )
    return gene_ids, means, counts


def depmap_selective(
    effects: Path,
    models: Path,
    common_essentials: Path,
    symbols: dict[str, str],
    lineages: dict[str, str],
) -> tuple[dict[str, set[str]], dict[str, int]]:
    """Per lineage, genes the lineage's models depend on more than other models do.

    A gene is selected when its mean effect in the lineage is at most DEPENDENCY, at least
    |SELECTIVITY| below its mean in all other models, and it is not a common essential.
    `lineages` maps DepMap OncotreeLineage to a set name. Returns sets and model counts.
    """
    with common_essentials.open() as stream:
        common = {symbols.get(_depmap_symbol(line.strip())) for line in stream} - {None}
    gene_ids, means, counts = depmap_lineages(effects, models, symbols, lineages)
    sets = {}
    for name, (mean_in, mean_out) in means.items():
        chosen = (mean_in <= DEPENDENCY) & (mean_in - mean_out <= SELECTIVITY)
        sets[name] = {
            g for g, keep in zip(gene_ids, chosen, strict=True) if keep and g and g not in common
        }
    return sets, counts


def held_out_means(
    links: csr_array, values: np.ndarray, known: np.ndarray, groups: np.ndarray, folds: int = 5
) -> np.ndarray:
    """Each gene's mean value over its neighbours with training-fold values; NaN if none."""
    estimates = np.full(len(values), np.nan)
    for train, test in GroupKFold(folds).split(values, groups=groups):
        mask = np.zeros(len(values))
        mask[train] = known[train]
        filled = np.where(mask > 0, values, 0.0)
        count = links @ mask
        total = links @ filled
        estimates[test] = np.divide(
            total, count, out=np.full(len(values), np.nan), where=count > 0
        )[test]
    return estimates


def _open_text(path: Path):
    if path.suffix == ".zip":
        archive = zipfile.ZipFile(path)
        return io.TextIOWrapper(archive.open(archive.namelist()[0]), encoding="utf-8")
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def gwas_traits(path: Path, symbols: dict[str, str], traits: set[str]) -> dict[str, set[str]]:
    """Mapped genes of genome-wide significant GWAS Catalog associations per EFO trait."""
    csv.field_size_limit(1 << 30)
    sets: dict[str, set[str]] = {trait: set() for trait in traits}
    with _open_text(path) as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        trait_column = header.index("MAPPED_TRAIT")
        gene_column = header.index("MAPPED_GENE")
        p_column = header.index("P-VALUE")
        for row in reader:
            try:
                if float(row[p_column]) > GWAS_SIGNIFICANCE:
                    continue
            except (ValueError, IndexError):
                continue
            hits = [t for t in row[trait_column].split(", ") if t in traits]
            if not hits:
                continue
            for symbol in re.split(r"[,;]\s*| - ", row[gene_column]):
                gene = symbols.get(symbol.strip())
                if gene:
                    for trait in hits:
                        sets[trait].add(gene)
    return sets
