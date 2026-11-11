"""Do Gene Ontology annotations agree with the pathways guessed for uncharacterised genes?

Each guess in `uncharacterised.tsv` names a Reactome pathway for a C*orf* gene that has no
Reactome pathway. A gene's GO biological-process terms (propagated through is_a and
part_of; specific terms only) are compared with each pathway's members: the similarity is
the mean share of members carrying each of the gene's terms. The guessed pathway's rank
among all pathways is then compared with chance (median percentile 0.5).

Annotations that could repeat the network or Reactome are left out: electronic (IEA),
expression-pattern (IEP) and no-data (ND) evidence, NOT qualifiers, and those assigned by
Reactome.
"""

import csv
import gzip
import json
from pathlib import Path

import numpy as np

from .gene_sets import read_gmt
from .references import GeneMap

EXCLUDED_EVIDENCE = {"IEA", "IEP", "ND"}
SPECIFIC = 0.05


def read_ontology(path: Path) -> dict[str, set[str]]:
    """Parents of each biological-process term through is_a and part_of."""
    parents: dict[str, set[str]] = {}
    term, namespace, links, obsolete = None, None, set(), False

    def flush():
        if term and namespace == "biological_process" and not obsolete:
            parents[term] = links

    with path.open() as stream:
        for line in stream:
            line = line.rstrip("\n")
            if line.startswith("["):
                flush()
                term, namespace, links, obsolete = None, None, set(), False
                continue
            if line.startswith("id: GO:"):
                term = line[4:]
            elif line.startswith("namespace: "):
                namespace = line[11:]
            elif line.startswith("is_a: "):
                links.add(line[6:16])
            elif line.startswith("relationship: part_of "):
                links.add(line[22:32])
            elif line == "is_obsolete: true":
                obsolete = True
        flush()
    return parents


def ancestors(term: str, parents: dict[str, set[str]], cache: dict[str, set[str]]) -> set[str]:
    if term not in cache:
        found = {term}
        for parent in parents.get(term, ()):
            found |= ancestors(parent, parents, cache)
        cache[term] = found
    return cache[term]


def read_annotations(path: Path, lookup, parents: dict[str, set[str]]) -> dict[str, set[str]]:
    """Propagated biological-process terms per gene ID."""
    cache: dict[str, set[str]] = {}
    terms: dict[str, set[str]] = {}
    with gzip.open(path, "rt") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            if row[0].startswith("!") or row[8] != "P" or "NOT" in row[3]:
                continue
            if row[6] in EXCLUDED_EVIDENCE or row[14] == "Reactome" or row[4] not in parents:
                continue
            gene = lookup(row[2])
            if gene:
                terms.setdefault(gene, set()).update(ancestors(row[4], parents, cache))
    return terms


def percentiles(
    guesses: dict[str, str], terms: dict[str, set[str]], pathways: dict[str, set[str]]
) -> dict[str, float]:
    """Per gene, the share of pathways scoring below the guessed one (ties count half)."""
    frequency: dict[str, int] = {}
    for found in terms.values():
        for term in found:
            frequency[term] = frequency.get(term, 0) + 1
    specific = {t for t, n in frequency.items() if n <= SPECIFIC * len(terms)}
    names = [p for p in pathways if any(m in terms for m in pathways[p])]
    profiles = {}
    for name in names:
        members = [terms[m] for m in pathways[name] if m in terms]
        counts: dict[str, int] = {}
        for found in members:
            for term in found & specific:
                counts[term] = counts.get(term, 0) + 1
        profiles[name] = (counts, len(members))
    result = {}
    for gene, guess in guesses.items():
        own = terms.get(gene, set()) & specific
        if not own or guess not in profiles:
            continue
        scores = np.array(
            [sum(profiles[n][0].get(t, 0) for t in own) / profiles[n][1] / len(own) for n in names]
        )
        mine = scores[names.index(guess)]
        result[gene] = float(
            ((scores < mine).sum() + 0.5 * ((scores == mine).sum() - 1)) / (len(names) - 1)
        )
    return result


def run(references: Path, guesses_path: Path, output: Path) -> dict:
    genes = GeneMap(references / "hgnc/hgnc_complete_set.txt")
    with guesses_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    guesses = {row["gene"]: row["pathway"] for row in rows}
    with (references / "hgnc/hgnc_complete_set.txt").open(newline="") as stream:
        symbols = {
            row["symbol"]: row["ensembl_gene_id"]
            for row in csv.DictReader(stream, delimiter="\t")
            if row["ensembl_gene_id"]
        }
    parents = read_ontology(references / "gene-ontology/go-basic.obo")
    terms = read_annotations(
        references / "gene-ontology/goa_human.gaf.gz", genes.by_symbol, parents
    )
    pathways = read_gmt(references / "reactome/ReactomePathways.gmt.zip", symbols)
    ranks = percentiles(guesses, terms, pathways)
    values = np.array(list(ranks.values()))
    rng = np.random.default_rng(0)
    means = [values[rng.integers(0, len(values), len(values))].mean() for _ in range(2000)]
    precision = {row["gene"]: float(row["expected_precision"]) for row in rows}
    confident = np.array([ranks[g] for g in ranks if precision[g] >= 0.5])
    report = {
        "guesses": len(guesses),
        "with_specific_go_terms": len(ranks),
        "mean_percentile": float(values.mean()),
        "mean_percentile_ci95": [float(q) for q in np.quantile(means, [0.025, 0.975])],
        "share_in_top_tenth": float((values >= 0.9).mean()),
        "confident_guesses": int(len(confident)),
        "confident_mean_percentile": float(confident.mean()) if len(confident) else None,
        "genes": {
            g: {"pathway": guesses[g], "percentile": ranks[g]}
            for g in sorted(ranks, key=lambda g: -ranks[g])
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "ontology.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
