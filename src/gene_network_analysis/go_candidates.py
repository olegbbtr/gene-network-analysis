"""Functions guessed today for genes with no Gene Ontology process annotation yet.

The method of `prospective`, applied to the current annotation instead of the pre-2018 one:
every biological-process annotation in the GAF counts (any date, any evidence), propagated
through is_a and part_of; terms with SMALLEST to LARGEST annotated genes are scored. A gene's
score for a term is the share of its network neighbours annotated with it, averaged over the
networks. Genes with no process annotation at all are the candidates. A guess (gene, term) is kept
when the gene ranks in the top TOP_SHARE of the genes not annotated with that term (ties
count half), the level at which `prospective` found guesses later confirmed 4.2 times more
often than chance. Each candidate's guess is its kept term with the largest share of
neighbours carrying it; candidates are listed by that share, since rare terms put almost
any non-zero share at rank 1.
"""

import json
from pathlib import Path

import numpy as np

from .applications import load
from .ontology import read_ontology
from .prospective import (
    KINDS,
    LARGEST,
    SMALLEST,
    membership,
    neighbour_means,
    read_dated,
    read_names,
)
from .references import GeneMap

TOP = 25
TOP_SHARE = 0.01


def run(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    parents = read_ontology(references / "gene-ontology/go-basic.obo")
    names = read_names(references / "gene-ontology/go-basic.obo")
    dated = read_dated(
        references / "gene-ontology/goa_human.gaf.gz",
        GeneMap(references / "hgnc/hgnc_complete_set.txt").by_symbol,
        parents,
    )
    current: dict[str, set[str]] = {}
    for kind in KINDS:
        for gene, terms in dated[kind].items():
            current.setdefault(gene, set()).update(terms)
    counts: dict[str, int] = {}
    for gene in genes:
        for term in current.get(gene, ()):
            counts[term] = counts.get(term, 0) + 1
    terms = sorted(t for t, n in counts.items() if SMALLEST <= n <= LARGEST)
    annotated = membership(current, genes, terms)
    unannotated = np.array([not current.get(g) for g in genes])
    links = [context.links[n] for n in context.plain]
    degree = np.mean([np.asarray(m.sum(axis=1)).ravel() for m in links], axis=0)
    scores = neighbour_means(links, annotated.astype(float))
    candidates = np.flatnonzero(unannotated & (degree > 0))
    rows = []
    for i in candidates:
        best_rank, best = None, None
        for j in np.flatnonzero(scores[i] > 0):
            negative = ~annotated[:, j]
            column = scores[negative, j]
            rank = 1 + (column > scores[i, j]).sum() + 0.5 * ((column == scores[i, j]).sum() - 1)
            if rank > TOP_SHARE * negative.sum():
                continue
            if best is None or scores[i, j] > scores[i, best[0]]:
                best_rank, best = rank, (j, int(negative.sum()))
        if best is not None:
            j, n = best
            rows.append(
                {
                    "symbol": context.table.get(genes[i], {}).get("symbol", genes[i]),
                    "gene": genes[i],
                    "term": terms[j],
                    "name": names.get(terms[j], ""),
                    "rank": float(best_rank),
                    "candidates": n,
                    "score": float(scores[i, j]),
                    "mean_links": float(degree[i]),
                }
            )
    rows.sort(key=lambda r: -r["score"])
    report = {
        "genes_without_process_annotation": int(unannotated.sum()),
        "scored": len(rows),
        "terms": len(terms),
        "top": rows[:TOP],
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "go-candidates.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report
