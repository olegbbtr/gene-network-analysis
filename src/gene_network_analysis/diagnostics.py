"""Checks that explain the TF rankings: how senders were chosen, how TFs connect, which
genes top each ranking, and how the rankings change without one co-expressed TF family.
"""

import datetime
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from scipy.stats import fisher_exact

from .evaluation import normalize_gene_id, ranking_metrics, read_reference, resolve_labels
from .graph import read_graph
from .neighbours import neighbour_vote
from .pipelines import PIPELINES, Runner
from .references import read_gene_table
from .report import compare, revision, sha256

PROTEIN_CODING = "protein-coding gene"
EXCLUDED_GROUP = "KRAB domain containing"
FAMILY = "Zinc fingers C2H2-type"
COMPARED = ("out-degree", "repression-order")


def sender_selection(senders: np.ndarray, labels: np.ndarray, coding: np.ndarray) -> dict:
    """Is a protein-coding gene in the graph more likely to send edges if it is a TF?"""
    table = [
        [int((senders & coding & labels).sum()), int((senders & coding & ~labels).sum())],
        [int((~senders & coding & labels).sum()), int((~senders & coding & ~labels).sum())],
    ]
    odds, p = fisher_exact(table)
    return {
        "protein_coding_senders": table[0][0] + table[0][1],
        "tf_share_among_senders": table[0][0] / max(1, table[0][0] + table[0][1]),
        "protein_coding_target_only": table[1][0] + table[1][1],
        "tf_share_among_target_only": table[1][0] / max(1, table[1][0] + table[1][1]),
        "odds_ratio": float(odds) if np.isfinite(odds) else None,
        "fisher_p": float(p),
        "table": {"rows": ["sender", "target only"], "columns": ["TF", "other"], "counts": table},
    }


def tf_partner_share(adjacency: csr_array, universe: np.ndarray, labels: np.ndarray) -> dict:
    """Share of outgoing weight to TFs, for TF senders against other senders."""
    out = neighbour_vote(adjacency, labels, "out")[universe]
    tf = labels[universe]
    return {
        "tf_senders": float(out[tf].mean()),
        "other_senders": float(out[~tf].mean()),
        "tf_share_of_all_targets": float(labels.mean()),
    }


def chromosome(location: str) -> str:
    match = re.match(r"\d+|X|Y|MT", location)
    return match.group() if match else "unknown"


def top_genes(scores: np.ndarray, labels: np.ndarray, nodes: list, genes: dict, k: int) -> dict:
    """The k first-ranked genes (ties by node ID, as in ranking_metrics), annotated."""
    order = np.lexsort((np.asarray(nodes), -scores))[:k]
    listed = []
    for i in order:
        info = genes.get(normalize_gene_id(nodes[i]), {})
        listed.append(
            {
                "gene": nodes[i],
                "symbol": info.get("symbol"),
                "location": info.get("location"),
                "tf": bool(labels[i]),
                "family": FAMILY in info.get("gene_group", []),
                "gene_groups": info.get("gene_group", []),
            }
        )
    return {
        "reference_genes": sum(g["tf"] for g in listed),
        "family": FAMILY,
        "family_genes": sum(g["family"] for g in listed),
        "family_reference_genes": sum(g["family"] and g["tf"] for g in listed),
        "chromosomes": dict(
            Counter(chromosome(g["location"] or "") for g in listed).most_common(5)
        ),
        "genes": listed,
    }


def evaluate(names, scores, labels, nodes, genes, records, top_k) -> dict:
    """Comparison rows and top genes for one gene set."""
    records = [dict(record, metrics={}) for record in records]
    for record in records:
        record["metrics"]["set"] = ranking_metrics(scores[record["name"]], labels, nodes, top_k)
    return {
        "comparison": compare(records, scores, "set", labels, top_k),
        "top_genes": {name: top_genes(scores[name], labels, nodes, genes, top_k) for name in names},
    }


def diagnose(
    graph_path: Path,
    references: dict[str, Path],
    gene_table: Path,
    output: Path,
    top_k: int = 100,
) -> dict:
    graph = read_graph(graph_path)
    raw = graph.adjacency
    nodes = list(graph.nodes)
    _, known = resolve_labels(graph.nodes, set())
    senders = np.diff(raw.indptr) > 0
    universe = np.flatnonzero(senders & known)
    universe_nodes = [nodes[i] for i in universe]
    genes = read_gene_table(gene_table)
    info = [genes.get(normalize_gene_id(node), {}) for node in nodes]
    coding = np.array([g.get("locus_group") == PROTEIN_CODING for g in info])
    excluded = np.array([EXCLUDED_GROUP in g.get("gene_group", []) for g in info])[universe]
    kept = ~excluded
    kept_nodes = [node for node, keep in zip(universe_nodes, kept, strict=True) if keep]
    selected = [p for p in PIPELINES if p.name in COMPARED or (p.uses_labels and p.reported)]
    runner = Runner(raw)
    records = [runner.record(p) for p in selected]
    names = [p.name for p in selected]
    results = {}
    for name, path in references.items():
        labels = resolve_labels(graph.nodes, read_reference(path))[0]
        scores = {p.name: runner.score(p, universe, labels) for p in selected}
        results[name] = {
            "path": path.name,
            "sha256": sha256(path),
            "sender_selection": sender_selection(senders, labels, coding),
            "tf_partner_share": tf_partner_share(raw, universe, labels),
            "all_senders": evaluate(
                names, scores, labels[universe], universe_nodes, genes, records, top_k
            ),
            "excluded_genes": int(excluded.sum()),
            "excluded_reference_genes": int((excluded & labels[universe]).sum()),
            "without_excluded": evaluate(
                names,
                {key: score[kept] for key, score in scores.items()},
                labels[universe][kept],
                kept_nodes,
                genes,
                records,
                top_k,
            ),
        }
    report = {
        "date": datetime.date.today().isoformat(),
        "implementation": revision(),
        "input": {
            "graph": {"path": graph_path.name, "sha256": sha256(graph_path)},
            "gene_table": {"path": gene_table.name, "sha256": sha256(gene_table)},
        },
        "selection": {
            "genes": len(universe),
            "top_k": top_k,
            "protein_coding": f"HGNC locus group '{PROTEIN_CODING}'",
            "excluded": f"HGNC gene group '{EXCLUDED_GROUP}', a TF family expressed together",
            "family": f"HGNC gene group '{FAMILY}', counted in top genes",
        },
        "references": results,
        "notes": [
            "Pipelines that use labels score each gene from other genes' labels, never its own.",
            "The best baseline is out-degree; repression-order is the best label-free pipeline.",
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "diagnostics.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report
