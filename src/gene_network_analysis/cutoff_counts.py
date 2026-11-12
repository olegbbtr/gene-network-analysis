"""Link counts at cutoffs taken from each network's own |z| distribution.

The public networks keep links with |z| > 2.35, a cutoff inherited from the 2016 graph. Here
each network's full z-score matrix is recomputed from its expression data exactly as in
`simulate.infer`, without a cutoff, and every gene's outgoing and incoming links are counted
above the 25th, 50th and 75th percentiles of |z| over all gene pairs of that network. The
counts (log1p) form an extra block of the local graph group in `label_free`.
"""

import json
from pathlib import Path

import numpy as np

from .public import DATASETS, read_dataset, read_expression, rows, standardise
from .references import read_gene_table
from .simulate import association

PERCENTILES = (25, 50, 75)
CACHE = "cutoff-counts.npz"


def zscores(data: np.ndarray) -> np.ndarray:
    """Every gene a sender: the z matrix of `simulate.infer` before its cutoff (float32)."""
    correlation = association(data, "correlation", 0.0).astype(np.float32)
    fisher = np.arctanh(np.clip(correlation, -0.999999, 0.999999))
    del correlation
    fisher -= fisher.mean(axis=0)
    fisher /= fisher.std(axis=0)
    np.fill_diagonal(fisher, 0)
    return fisher


def counts(z: np.ndarray) -> tuple[np.ndarray, list[float]]:
    """Out- and in-link counts per gene above each percentile of off-diagonal |z|."""
    magnitude = np.abs(z)
    off = ~np.eye(len(z), dtype=bool)
    sample = magnitude[off]
    cutoffs = [float(c) for c in np.percentile(sample, PERCENTILES)]
    del sample
    columns = []
    for cutoff in cutoffs:
        above = magnitude > cutoff  # the diagonal is 0 and never above a positive cutoff
        columns += [above.sum(axis=1), above.sum(axis=0)]
    return np.log1p(np.column_stack(columns).astype(np.float64)), cutoffs


def build(references: Path, networks: Path, gene_table: Path) -> dict:
    genes, _ = read_expression(networks / "expression.tsv")
    symbols = {row["symbol"]: gene for gene, row in read_gene_table(gene_table).items()}
    blocks, cutoffs = {}, {}
    for dataset in DATASETS:
        values = rows(read_dataset(dataset, references, symbols), genes)
        z = zscores(standardise(values))
        blocks[dataset.name], cutoffs[dataset.name] = counts(z)
        del z
        print(dataset.name, cutoffs[dataset.name], flush=True)
    np.savez(
        networks / CACHE,
        genes=np.array(genes),
        cutoffs=json.dumps(cutoffs),
        **blocks,
    )
    return {"percentiles": PERCENTILES, "cutoffs": cutoffs}


def load(networks: Path, genes: list[str], names: list[str]) -> np.ndarray:
    saved = np.load(networks / CACHE)
    if saved["genes"].tolist() != genes:
        raise ValueError(f"{networks / CACHE} was built for a different gene order")
    return np.column_stack([saved[name] for name in names])


def cutoffs(networks: Path) -> dict:
    return json.loads(str(np.load(networks / CACHE)["cutoffs"]))
