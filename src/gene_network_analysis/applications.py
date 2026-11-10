"""Uses of per-gene network neighbourhoods, on the networks written by `public.build`.

- function: do neighbours predict HGNC gene groups and Reactome pathways?
- uncharacterised: top pathway guesses for C*orf* genes, with precision expected from
  the same held-out scores on known genes;
- matched tissue: do DepMap lineage-selective dependencies and tissue GWAS traits score
  best on the network of the matching tissue?
- rewiring: do cancer genes change their closest partners between healthy and tumour
  tissue more than expression and degree explain?

Every score hides set members on held-out chromosomes (`gene_sets.held_out_scores`).
"""

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from sklearn.metrics import roc_auc_score

from .diagnostics import chromosome
from .gene_sets import (
    depmap_selective,
    gwas_traits,
    held_out_scores,
    hgnc_groups,
    membership,
    neighbour_scores,
    personalised_pagerank,
    read_gmt,
    restrict,
    set_aucs,
)
from .graph import read_graph
from .pipelines import CLEAN_2016
from .prediction import aligned, cross_validated, paired_gain, read_gene_list, undirected
from .public import read_expression
from .references import read_gene_table

LINEAGES = {
    "Lung": "lung",
    "Breast": "breast",
    "Bowel": "colon",
    "Thyroid": "thyroid",
    "Prostate": "prostate",
    "Esophagus/Stomach": "stomach",
    "Liver": "liver",
}
GWAS_TRAITS = {
    "forced expiratory volume": "lung",
    "chronic obstructive pulmonary disease": "lung",
    "lung carcinoma": "lung",
    "breast carcinoma": "breast",
    "mammographic density measurement": "breast",
    "colorectal cancer": "colon",
    "ulcerative colitis": "colon",
    "hypothyroidism": "thyroid",
    "prostate carcinoma": "prostate",
    "benign prostatic hyperplasia": "prostate",
    "gastroesophageal reflux disease": "stomach",
    "peptic ulcer disease": "stomach",
    "gastric carcinoma": "stomach",
    "aspartate aminotransferase measurement": "liver",
    "hepatocellular carcinoma": "liver",
    "body height": None,
}
UNCHARACTERISED = re.compile(r"^C(\d+|X|Y)orf\d+$")
PARTNERS = 50


@dataclass
class Context:
    genes: list[str]
    table: dict[str, dict]
    symbols: dict[str, str]
    groups: np.ndarray
    stats: dict
    expression: dict[str, np.ndarray]
    adjacency: dict[str, csr_array]
    links: dict[str, csr_array]

    @property
    def plain(self) -> list[str]:
        return [n for n in self.links if not n.endswith("-copy-number-removed")]


def load(root: Path, gene_table: Path) -> Context:
    table = read_gene_table(gene_table)
    genes, expression = read_expression(root / "expression.tsv")
    stats = json.loads((root / "stats.json").read_text())
    names = [n for n in stats if n != "settings"]
    adjacency = {n: aligned(read_graph(root / n / "graph.txt"), genes) for n in names}
    return Context(
        genes,
        table,
        {row["symbol"]: gene for gene, row in table.items()},
        np.array([chromosome(table.get(g, {}).get("location", "")) for g in genes]),
        stats,
        expression,
        adjacency,
        {n: undirected(a) for n, a in adjacency.items()},
    )


def network_scores(context: Context, labels: np.ndarray) -> dict[str, np.ndarray]:
    """Held-out scores per network, plus `combined`: their mean over plain networks."""
    scores = {n: held_out_scores(m, labels, context.groups) for n, m in context.links.items()}
    scores["combined"] = np.mean([scores[n] for n in context.plain], axis=0)
    return scores


def summarise(aucs: np.ndarray) -> dict:
    return {
        "sets": int(len(aucs)),
        "median_auc": float(np.median(aucs)),
        "share_above_0.6": float((aucs > 0.6).mean()),
        "share_above_0.7": float((aucs > 0.7).mean()),
    }


def function(context: Context, sets: dict[str, dict[str, set[str]]]) -> dict:
    report = {}
    for source, raw in sets.items():
        chosen = restrict(raw, context.genes)
        labels = membership(context.genes, chosen)
        scores = network_scores(context, labels)
        aucs = {n: set_aucs(s, labels) for n, s in scores.items()}
        order = np.argsort(-aucs["combined"])
        names = list(chosen)
        report[source] = {
            "networks": {n: summarise(a) for n, a in aucs.items()},
            "best": [(names[i], float(aucs["combined"][i])) for i in order[:15]],
            "worst": [(names[i], float(aucs["combined"][i])) for i in order[-10:]],
        }
    return report


def uncharacterised(context: Context, sets: dict[str, set[str]], bins: int = 10) -> dict:
    """Top pathway per gene from held-out combined scores, with binned expected precision."""
    chosen = restrict(sets, context.genes)
    names = list(chosen)
    labels = membership(context.genes, chosen)
    scores = network_scores(context, labels)["combined"]
    top = scores.argmax(axis=1)
    best = scores[np.arange(len(top)), top]
    annotated = labels.any(axis=1)
    hit = labels[np.arange(len(top)), top]
    edges = np.quantile(best[annotated], np.linspace(0, 1, bins + 1))
    edges[-1] = np.inf
    level = np.clip(np.searchsorted(edges, best, side="right") - 1, 0, bins - 1)
    precision = np.array([hit[annotated & (level == b)].mean() for b in range(bins)])
    rows = []
    for i, gene in enumerate(context.genes):
        symbol = context.table.get(gene, {}).get("symbol", "")
        if UNCHARACTERISED.match(symbol) and not annotated[i]:
            rows.append(
                {
                    "gene": gene,
                    "symbol": symbol,
                    "pathway": names[top[i]],
                    "score": float(best[i]),
                    "expected_precision": float(precision[level[i]]),
                }
            )
    rows.sort(key=lambda r: -r["score"])
    return {
        "pathways": len(names),
        "annotated_genes": int(annotated.sum()),
        "top1_precision_known_genes": float(hit[annotated].mean()),
        "precision_by_score_decile": [float(p) for p in precision],
        "genes": rows,
    }


def matched(context: Context, sets: dict[str, set[str]], tissue: dict[str, str | None]) -> dict:
    """AUC of each set on each network; matched-tissue networks against the others."""
    names = [n for n in sets if len(sets[n] & set(context.genes)) >= 10]
    labels = membership(context.genes, {n: sets[n] for n in names})
    scores = network_scores(context, labels)
    report = {}
    for column, name in enumerate(names):
        y = labels[:, column]
        aucs = {n: float(roc_auc_score(y, s[:, column])) for n, s in scores.items()}
        entry: dict = {"tissue": tissue.get(name), "members": int(y.sum()), "auc": aucs}
        for state in ("healthy", "tumour"):
            pool = [n for n in context.plain if context.stats[n]["state"] == state]
            own = [n for n in pool if context.stats[n]["tissue"] == tissue.get(name)]
            other = [n for n in pool if n not in own]
            if not own:
                continue
            ranking = sorted(pool, key=lambda n: -aucs[n])
            entry[state] = {
                "matched": own,
                "matched_rank": min(ranking.index(n) for n in own) + 1,
                "networks": len(pool),
                "matched_vs_others": paired_gain(
                    y,
                    np.mean([scores[n][:, column] for n in own], axis=0),
                    np.mean([scores[n][:, column] for n in other], axis=0),
                    runs=500,
                ),
            }
        report[name] = entry
    return report


def closest(adjacency: csr_array, k: int = PARTNERS) -> list[set[int]]:
    """Each gene's k strongest partners in either direction."""
    strength = abs(adjacency).maximum(abs(adjacency).T).tocsr()
    partners = []
    for row in range(strength.shape[0]):
        start, end = strength.indptr[row], strength.indptr[row + 1]
        cols, weights = strength.indices[start:end], strength.data[start:end]
        partners.append(set(cols[np.argsort(-weights)[:k]].tolist()))
    return partners


def rewiring_scores(healthy: csr_array, tumour: csr_array) -> np.ndarray:
    """1 - Jaccard overlap of each gene's closest partners in the two networks."""
    before, after = closest(healthy), closest(tumour)
    return np.array(
        [1 - len(a & b) / len(a | b) if a | b else 0.0 for a, b in zip(before, after, strict=True)]
    )


def rewiring(context: Context, lists: dict[str, set[str]]) -> dict:
    report = {}
    pairs = [
        (h, t)
        for t in context.plain
        if context.stats[t]["state"] == "tumour"
        for h in context.plain
        if context.stats[h]["state"] == "healthy"
        and context.stats[h]["tissue"] == context.stats[t]["tissue"]
    ]
    for healthy, tumour in pairs:
        change = rewiring_scores(context.adjacency[healthy], context.adjacency[tumour])
        quantile = {
            n: context.expression[f"{n}_mean"].argsort().argsort() / len(context.genes)
            for n in (healthy, tumour)
        }
        degree = [
            np.log1p(np.asarray(context.links[n].sum(axis=1)).ravel()) for n in (healthy, tumour)
        ]
        blocks = {
            "baseline": np.column_stack(
                [context.expression[f"{n}_{s}"] for n in (healthy, tumour) for s in ("mean", "sd")]
                + [quantile[tumour] - quantile[healthy], *degree]
            ),
            "rewiring": change[:, None],
        }
        entry = {"median_rewiring": float(np.median(change))}
        for name, members in lists.items():
            y = np.array([g in members for g in context.genes])
            base = cross_validated(blocks, ("baseline",), y, context.groups, [])
            full = cross_validated(blocks, ("baseline", "rewiring"), y, context.groups, [])
            entry[name] = {
                "rewiring_alone_auc": float(roc_auc_score(y, change)),
                "baseline_auc": float(roc_auc_score(y, base)),
                "with_rewiring_auc": float(roc_auc_score(y, full)),
                "gain": paired_gain(y, full, base, runs=500),
            }
        report[f"{healthy}->{tumour}"] = entry
    return report


def cleanings(adjacency: csr_array) -> dict[str, csr_array]:
    """Undirected link sets after each cleaning of one signed network."""
    present = (adjacency != 0).astype(np.float64)
    positive = (adjacency > 0).astype(np.float64)
    cleaned = adjacency
    for clean in CLEAN_2016:
        cleaned = clean(cleaned)
    rows = [(i, j) for i, partners in enumerate(closest(adjacency)) for j in partners]
    strongest = csr_array((np.ones(len(rows)), tuple(np.array(rows).T)), shape=adjacency.shape)
    return {
        "all": undirected(adjacency),
        "two-way": undirected(present.multiply(present.T)),
        "two-way-positive": undirected(positive.multiply(positive.T)),
        f"top-{PARTNERS}": undirected(strongest),
        "2016": undirected(cleaned),
    }


def cleaning(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
    """Essential genes: neighbour share and personalised PageRank on each cleaning."""
    context = load(root, gene_table)
    essential = read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols)
    y = np.array([g in essential for g in context.genes])
    labels = y[:, None]
    scorers = {"neighbours": neighbour_scores, "pagerank": personalised_pagerank}
    report = {}
    for name in context.links:
        links = cleanings(context.adjacency[name])
        scores = {
            (variant, scorer): held_out_scores(m, labels, context.groups, scorer=f)[:, 0]
            for variant, m in links.items()
            for scorer, f in scorers.items()
        }
        baseline = scores["all", "neighbours"]
        report[name] = {
            f"{variant}/{scorer}": {
                "links": int(links[variant].nnz // 2),
                "auc": float(roc_auc_score(y, s)),
                "vs_all_neighbours": paired_gain(y, s, baseline, runs=300),
            }
            for (variant, scorer), s in scores.items()
        }
    output.mkdir(parents=True, exist_ok=True)
    (output / "cleaning.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def run(root: Path, gene_table: Path, references: Path, output: Path) -> dict:
    context = load(root, gene_table)
    reactome = read_gmt(references / "reactome/ReactomePathways.gmt.zip", context.symbols)
    depmap, models = depmap_selective(
        references / "depmap-24q2/CRISPRGeneEffect.csv",
        references / "depmap-24q2/Model.csv",
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv",
        context.symbols,
        LINEAGES,
    )
    gwas = gwas_traits(
        references / "gwas-catalog/gwas-catalog-associations_ontology-annotated-full.zip",
        context.symbols,
        set(GWAS_TRAITS),
    )
    lists = {
        "cancer": read_gene_list(
            references / "gene-lists/oncokb-cancerGeneList.txt", context.symbols
        ),
        "essential": read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols),
    }
    report = {
        "function": function(context, {"hgnc": hgnc_groups(context.table), "reactome": reactome}),
        "uncharacterised": uncharacterised(context, reactome),
        "dependencies": {
            "models": models,
            "sets": matched(context, depmap, {name: name for name in depmap}),
        },
        "gwas": matched(context, gwas, GWAS_TRAITS),
        "rewiring": rewiring(context, lists),
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "applications.json").write_text(json.dumps(report, indent=2) + "\n")
    with (output / "uncharacterised.tsv").open("w") as stream:
        stream.write("symbol\tgene\tpathway\tscore\texpected_precision\n")
        for row in report["uncharacterised"]["genes"]:
            stream.write(
                f"{row['symbol']}\t{row['gene']}\t{row['pathway']}\t{row['score']:.4f}\t"
                f"{row['expected_precision']:.3f}\n"
            )
    return report
