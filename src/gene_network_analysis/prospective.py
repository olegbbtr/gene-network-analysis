"""Do a gene's network neighbours predict functions that were published later?

Gene Ontology biological-process annotations are split by the PubMed ID of their
reference, a stand-in for publication date: PubMed IDs are given in order of entry, and
CUTOFF was reached in about August 2018, after every expression sample used here had been
collected.

- before: any annotation citing a PubMed ID below CUTOFF. The scores see only these.
- after: experimental annotations (EXPERIMENTAL; expression-pattern evidence is left out)
  citing a PubMed ID from CUTOFF on.
- A gene is a positive for a term if it gained the term after and nothing else supports
  it: no annotation before and no undated one (electronic, phylogenetic and database
  inferences carry no PubMed ID). Genes with such support, or with only non-experimental
  support after, are left out for that term; every other gene is a negative.

Terms are propagated through is_a and part_of; the scored terms have SMALLEST to LARGEST
genes annotated before. A gene's score for a term is the share of its network neighbours
annotated with it before, averaged over the networks. Two rankings that do not depend on
the term are scored the same way: how many terms the gene had before (how well studied it
was) and how many links it has.

Each positive gets the percentile of its score among the term's negatives (ties count
half), so a term's mean percentile is its ROC AUC and 0.5 is chance. Results are given for
all genes, for genes that received any experimental annotation after (this controls for
which genes were studied), and for those of them with no annotation before. Intervals are
bootstraps over genes, since a gene's terms are nested.

The cutoff, evidence codes, term sizes and universes were fixed before the first run.
"""

import gzip
import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from .applications import load
from .gene_sets import neighbour_scores
from .ontology import ancestors, read_ontology
from .references import GeneMap

CUTOFF = 30_000_000
EXPERIMENTAL = frozenset({"EXP", "IDA", "IPI", "IMP", "IGI", "HTP", "HDA", "HMP", "HGI"})
KINDS = ("before", "after", "later", "undated")
SMALLEST, LARGEST = 20, 500
MIN_POSITIVES = 10
TOP_SHARE = 0.01
EXAMPLES = 40


def kind(reference: str, evidence: str) -> str:
    """When and how a GAF row states its annotation: one of KINDS."""
    pmids = [int(r[5:]) for r in reference.split("|") if r.startswith("PMID:")]
    if not pmids:
        return "undated"
    if min(pmids) < CUTOFF:
        return "before"
    return "after" if evidence in EXPERIMENTAL else "later"


def read_dated(path: Path, lookup, parents: dict[str, set[str]]) -> dict[str, dict[str, set[str]]]:
    """Per kind, the propagated biological-process terms of each gene ID."""
    cache: dict[str, set[str]] = {}
    terms: dict[str, dict[str, set[str]]] = {name: {} for name in KINDS}
    with gzip.open(path, "rt") as stream:
        for line in stream:
            row = line.rstrip("\n").split("\t")
            if line.startswith("!") or row[8] != "P" or "NOT" in row[3] or row[4] not in parents:
                continue
            gene = lookup(row[2])
            if gene:
                found = terms[kind(row[5], row[6])].setdefault(gene, set())
                found.update(ancestors(row[4], parents, cache))
    return terms


def read_names(path: Path) -> dict[str, str]:
    names, term = {}, None
    with path.open() as stream:
        for line in stream:
            if line.startswith("id: "):
                term = line[4:].strip()
            elif line.startswith("name: ") and term:
                names[term] = line[6:].strip()
                term = None
    return names


@dataclass(frozen=True)
class Split:
    terms: list[str]
    before: np.ndarray  # genes x terms
    positive: np.ndarray
    negative: np.ndarray
    studied: np.ndarray  # per gene: any experimental annotation after
    unannotated: np.ndarray  # per gene: no annotation before
    annotations: np.ndarray  # per gene: number of terms before


def membership(terms: dict[str, set[str]], genes: list[str], names: list[str]) -> np.ndarray:
    index = {term: j for j, term in enumerate(names)}
    matrix = np.zeros((len(genes), len(names)), dtype=bool)
    for i, gene in enumerate(genes):
        matrix[i, [index[t] for t in terms.get(gene, ()) if t in index]] = True
    return matrix


def split(dated: dict[str, dict[str, set[str]]], genes: list[str]) -> Split:
    counts = Counter(t for g in genes for t in dated["before"].get(g, ()))
    names = sorted(t for t, n in counts.items() if SMALLEST <= n <= LARGEST)
    before, after, later, undated = (membership(dated[k], genes, names) for k in KINDS)
    known = before | undated
    return Split(
        names,
        before,
        after & ~known,
        ~(known | after | later),
        np.array([g in dated["after"] for g in genes]),
        np.array([g not in dated["before"] for g in genes]),
        np.array([len(dated["before"].get(g, ())) for g in genes]),
    )


def percentiles(scores: np.ndarray, positive: np.ndarray, negative: np.ndarray) -> np.ndarray:
    """Share of the negatives scoring below each positive; ties count half."""
    ordered = np.sort(scores[negative])
    below = np.searchsorted(ordered, scores[positive], side="left")
    upto = np.searchsorted(ordered, scores[positive], side="right")
    return (below + 0.5 * (upto - below)) / len(ordered)


def top_hits(scores: np.ndarray, positive: np.ndarray, k: int) -> float:
    """Expected positives among the k highest scores when ties are broken at random."""
    threshold = np.partition(scores, len(scores) - k)[len(scores) - k]
    above, tied = scores > threshold, scores == threshold
    return float(positive[above].sum() + positive[tied].sum() * (k - above.sum()) / tied.sum())


@dataclass(frozen=True)
class Scored:
    gene: np.ndarray  # one entry per positive gene-term pair
    term: np.ndarray
    percentile: np.ndarray
    hits: float  # positives among each term's top TOP_SHARE of candidates
    expected: float


def score_pairs(scores: np.ndarray, data: Split, universe: np.ndarray) -> Scored:
    """`scores` is genes x terms, or one value per gene used for every term."""
    genes, terms, values = [], [], []
    hits = expected = 0.0
    for j in range(len(data.terms)):
        positive = data.positive[:, j] & universe
        negative = data.negative[:, j] & universe
        if not positive.any() or not negative.any():
            continue
        column = scores[:, j] if scores.ndim == 2 else scores
        found = np.flatnonzero(positive)
        genes.append(found)
        terms.append(np.full(len(found), j))
        values.append(percentiles(column, positive, negative))
        candidates = positive | negative
        k = math.ceil(TOP_SHARE * candidates.sum())
        hits += top_hits(column[candidates], positive[candidates], k)
        expected += k * positive.sum() / candidates.sum()
    if not genes:
        raise ValueError("no term has both a positive and a negative gene")
    return Scored(
        np.concatenate(genes), np.concatenate(terms), np.concatenate(values), hits, expected
    )


def gene_bootstrap(gene: np.ndarray, values: np.ndarray, runs: int, seed: int = 0) -> list[float]:
    """95% interval of the mean of `values` when whole genes are resampled."""
    _, index = np.unique(gene, return_inverse=True)
    sums, counts = np.bincount(index, weights=values), np.bincount(index)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(sums), (runs, len(sums)))
    means = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    return [float(q) for q in np.quantile(means, [0.025, 0.975])]


def summarise(scored: Scored, runs: int) -> dict:
    counts = np.bincount(scored.term)
    sums = np.bincount(scored.term, weights=scored.percentile)
    aucs = sums[counts >= MIN_POSITIVES] / counts[counts >= MIN_POSITIVES]
    return {
        "mean_percentile": float(scored.percentile.mean()),
        "ci95": gene_bootstrap(scored.gene, scored.percentile, runs),
        "terms_with_enough_positives": int(len(aucs)),
        "median_term_auc": float(np.median(aucs)) if len(aucs) else None,
        "share_of_terms_above_0.6": float((aucs > 0.6).mean()) if len(aucs) else None,
        "top_share_hits": scored.hits,
        "top_share_expected": scored.expected,
        "top_share_fold": scored.hits / scored.expected,
    }


def difference(first: Scored, second: Scored, runs: int) -> dict:
    gain = first.percentile - second.percentile
    return {"gain": float(gain.mean()), "ci95": gene_bootstrap(first.gene, gain, runs)}


def neighbour_means(links: list[csr_array], labels: np.ndarray) -> np.ndarray:
    total = np.zeros(labels.shape)
    for matrix in links:
        total += neighbour_scores(matrix, labels)
    return total / len(links)


def examples(
    scored: Scored, data: Split, chosen: np.ndarray, symbols: list[str], names: dict[str, str]
) -> list:
    """The best-ranked new term of each gene in `chosen`, for the genes ranked highest."""
    rows, seen = [], set()
    for i in np.argsort(-scored.percentile, kind="stable"):
        gene = int(scored.gene[i])
        if gene in seen or not chosen[gene]:
            continue
        seen.add(gene)
        column = scored.term[i]
        negatives = int(data.negative[:, column].sum())
        rows.append(
            {
                "symbol": symbols[gene],
                "term": data.terms[column],
                "name": names.get(data.terms[column], ""),
                "annotated_before": int(data.before[:, column].sum()),
                "rank": 1 + round((1 - float(scored.percentile[i])) * negatives),
                "candidates": negatives + 1,
            }
        )
        if len(rows) == EXAMPLES:
            break
    return rows


def run(root: Path, gene_table: Path, references: Path, output: Path, runs: int = 1000) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    parents = read_ontology(references / "gene-ontology/go-basic.obo")
    dated = read_dated(
        references / "gene-ontology/goa_human.gaf.gz",
        GeneMap(references / "hgnc/hgnc_complete_set.txt").by_symbol,
        parents,
    )
    data = split(dated, genes)
    links = {n: context.links[n] for n in context.plain}
    degree = np.mean([np.asarray(m.sum(axis=1)).ravel() for m in links.values()], axis=0)
    scores = {
        "neighbours": neighbour_means(list(links.values()), data.before),
        "terms_before": data.annotations.astype(np.float64),
        "links": degree,
    }
    universes = {
        "all": np.ones(len(genes), dtype=bool),
        "studied_after": data.studied,
        "first_annotated_after": data.studied & data.unannotated,
    }
    report: dict = {
        "cutoff_pmid": CUTOFF,
        "terms": len(data.terms),
        "genes": {
            "all": len(genes),
            "annotated_before": int((~data.unannotated).sum()),
            "studied_after": int(data.studied.sum()),
            "first_annotated_after": int((data.studied & data.unannotated).sum()),
        },
        "universes": {},
        "single_networks": {},
    }
    symbols = [context.table.get(g, {}).get("symbol", g) for g in genes]
    names = read_names(references / "gene-ontology/go-basic.obo")
    for name, universe in universes.items():
        scored = {s: score_pairs(values, data, universe) for s, values in scores.items()}
        entry: dict = {
            "pairs": int(len(scored["neighbours"].gene)),
            "genes_with_a_new_term": int(len(np.unique(scored["neighbours"].gene))),
            "terms_with_a_new_gene": int(len(np.unique(scored["neighbours"].term))),
        }
        entry |= {s: summarise(found, runs) for s, found in scored.items()}
        for baseline in ("terms_before", "links"):
            entry[f"neighbours_minus_{baseline}"] = difference(
                scored["neighbours"], scored[baseline], runs
            )
        if name == "all":
            report["examples"] = {
                u: examples(scored["neighbours"], data, chosen, symbols, names)
                for u, chosen in universes.items()
                if u != "studied_after"
            }
        report["universes"][name] = entry
        print(name, json.dumps(entry), flush=True)
    for name, matrix in links.items():
        single = neighbour_scores(matrix, data.before)
        report["single_networks"][name] = {
            u: float(score_pairs(single, data, universe).percentile.mean())
            for u, universe in universes.items()
        }
        print(name, json.dumps(report["single_networks"][name]), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    (output / "prospective.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
