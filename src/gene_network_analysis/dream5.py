"""Inputs for the pipelines from a DREAM5 network inference dataset (Marbach et al. 2012).

The graph is inferred from the expression data the same way as in `simulate.py`, with
every gene as a candidate sender, so no TF knowledge enters it. DREAM5 names genes G1, G2,
...; they are prefixed with ENSGD5 and the network number because the report ranks only
Ensembl-like IDs.
"""

import csv
from pathlib import Path

import numpy as np

from .graph import write_graph
from .simulate import infer


def gene_name(network: str, gene: str) -> str:
    return f"ENSGD5N{network}{gene}"


def read_expression(path: Path) -> tuple[list[str], np.ndarray]:
    """Samples x genes, each gene standardised; constant genes are dropped."""
    with path.open() as stream:
        genes = stream.readline().split()
        data = np.loadtxt(stream, delimiter="\t", ndmin=2)
    if data.shape[1] != len(genes):
        raise ValueError(f"{path}: {len(genes)} gene names but {data.shape[1]} columns")
    spread = data.std(axis=0)
    keep = spread > 0
    data = (data[:, keep] - data[:, keep].mean(axis=0)) / spread[keep]
    return [g for g, k in zip(genes, keep, strict=True) if k], data


def read_gold(path: Path, signed: Path | None) -> list[tuple[str, str, int]]:
    """True edges (weight 1 rows); signs from the signed file where it has one, else 0."""
    signs = {}
    if signed:
        with signed.open() as stream:
            for source, target, sign in csv.reader(stream, delimiter="\t"):
                signs[source, target] = {"+": 1, "-": -1}.get(sign, 0)
    with path.open() as stream:
        return [
            (source, target, signs.get((source, target), 0))
            for source, target, weight in csv.reader(stream, delimiter="\t")
            if weight == "1"
        ]


def prepare(
    network: str,
    expression: Path,
    tfs: Path,
    gold: Path,
    output: Path,
    *,
    signed: Path | None = None,
    inference: str = "correlation",
    threshold: float = 2.35,
) -> dict:
    genes, data = read_expression(expression)
    graph = infer(data, np.arange(len(genes)), threshold, inference)
    names = [gene_name(network, g) for g in genes]
    true_edges = read_gold(gold, signed)
    tf_list = tfs.read_text().split()
    regulators = sorted({source for source, _, _ in true_edges})
    output.mkdir(parents=True, exist_ok=True)
    write_graph(output / "graph.txt", graph, names)
    (output / "tfs.txt").write_text("".join(f"{gene_name(network, g)}\n" for g in tf_list))
    (output / "regulators.txt").write_text(
        "".join(f"{gene_name(network, g)}\n" for g in regulators)
    )
    with (output / "true-edges.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["source_gene", "target_gene", "sign"])
        for source, target, sign in true_edges:
            writer.writerow([gene_name(network, source), gene_name(network, target), sign])
    return {
        "genes": len(genes),
        "samples": int(data.shape[0]),
        "tfs": len(tf_list),
        "regulators_with_true_targets": len(regulators),
        "true_edges": len(true_edges),
        "inferred_edges": int(graph.nnz),
        "inferred_negative": int((graph.data < 0).sum()),
    }
