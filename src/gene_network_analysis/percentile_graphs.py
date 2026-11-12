"""Public networks that keep the strongest SHARE of gene pairs, per network.

Each network's full z matrix is recomputed as in `cutoff_counts` (every gene a sender, as in
`simulate.infer`) and links are kept where |z| is in the top `share` of off-diagonal pairs of
that network. Expression summaries and the cutoff-free link counts are shared with the
original networks.
"""

import json
import shutil
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from .cutoff_counts import CACHE, zscores
from .graph import write_graph
from .public import DATASETS, read_dataset, read_expression, rows, standardise
from .references import read_gene_table


def build(references: Path, networks: Path, gene_table: Path, shares: list[float]) -> dict:
    genes, _ = read_expression(networks / "expression.tsv")
    symbols = {row["symbol"]: gene for gene, row in read_gene_table(gene_table).items()}
    roots = {s: networks.parent / f"{networks.name}-p{round(s * 100):02d}" for s in shares}
    stats = {s: {} for s in shares}
    for root in roots.values():
        root.mkdir(parents=True, exist_ok=True)
        shutil.copy(networks / "expression.tsv", root / "expression.tsv")
        target = root / CACHE
        if not target.exists():
            target.symlink_to((networks / CACHE).resolve())
    original = json.loads((networks / "stats.json").read_text())
    for dataset in DATASETS:
        z = zscores(standardise(rows(read_dataset(dataset, references, symbols), genes)))
        magnitude = np.abs(z)
        off = ~np.eye(len(z), dtype=bool)
        for share, root in roots.items():
            cutoff = float(np.quantile(magnitude[off], 1 - share))
            r, c = np.nonzero(magnitude > cutoff)
            graph = csr_array((z[r, c].astype(np.float64), (r, c)), shape=z.shape)
            (root / dataset.name).mkdir(exist_ok=True)
            write_graph(root / dataset.name / "graph.txt", graph, genes)
            stats[share][dataset.name] = original[dataset.name] | {
                "edges": int(graph.nnz),
                "negative_edges": int((graph.data < 0).sum()),
                "cutoff": cutoff,
            }
            print(dataset.name, share, cutoff, graph.nnz, flush=True)
        del z, magnitude
    for share, root in roots.items():
        stats[share]["settings"] = original["settings"] | {"threshold": None, "share": share}
        (root / "stats.json").write_text(json.dumps(stats[share], indent=2) + "\n")
    return {str(s): str(r) for s, r in roots.items()}
