"""Do network descriptors predict genes that people cannot lose, beyond expression?

Essential-gene lists come from cancer cell lines. gnomAD v2.1.1 measures something else:
how many loss-of-function variants each gene carries in 141,456 people against the number
expected from its sequence. A gene is called constrained when the upper bound of that
ratio (LOEUF) is below LOEUF_CUTOFF, gnomAD's own threshold. Only genes with at least
MIN_EXPECTED expected variants are used, as gnomAD advises, because shorter genes cannot
reach the threshold.

The models are those of `robustness`, on held-out chromosomes. The baseline is the gene's
expression plus the log of its expected variant count, so gene length cannot pass as a
network effect. The comparison is repeated without the cell-essential genes (CEGv2 and
DepMap common essentials) to see whether the network only finds those again.

The threshold, the filter and the baseline were fixed before the first run.
"""

import csv
import gzip
import json
from pathlib import Path

import numpy as np

from .applications import load
from .graph import Graph
from .prediction import network_features, read_gene_list
from .robustness import blocks, compare, held_out
from .structure import ARM, location_codes

LOEUF_CUTOFF = 0.35
MIN_EXPECTED = 10.0


def read_constraint(path: Path) -> dict[str, tuple[float, float]]:
    """LOEUF and expected loss-of-function count per Ensembl gene ID, where both are given."""
    table = {}
    with gzip.open(path, "rt", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            try:
                table[row["gene_id"]] = (float(row["oe_lof_upper"]), float(row["exp_lof"]))
            except ValueError:
                continue
    return table


def labels(
    table: dict[str, tuple[float, float]], genes: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per gene: constrained, usable (enough expected variants), log expected count."""
    loeuf = np.array([table.get(g, (np.nan, np.nan))[0] for g in genes])
    expected = np.array([table.get(g, (np.nan, np.nan))[1] for g in genes])
    usable = expected >= MIN_EXPECTED
    return loeuf < LOEUF_CUTOFF, usable, np.log(np.where(usable, expected, 1.0))


def run(root: Path, gene_table: Path, references: Path, output: Path, runs: int = 500) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    arms = location_codes(tuple(genes), context.table, ARM)
    names = context.plain
    features = {
        n: network_features(Graph(tuple(genes), context.adjacency[n]), genes, arms) for n in names
    }
    constrained, usable, expected = labels(
        read_constraint(references / "gnomad-v2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"),
        genes,
    )
    essential = read_gene_list(references / "gene-lists/CEGv2.txt", context.symbols)
    essential |= read_gene_list(
        references / "depmap-24q2/CRISPRInferredCommonEssentials.csv", context.symbols
    )
    kept = {
        "all_genes": usable,
        "without_cell_essential": usable & np.array([g not in essential for g in genes]),
    }
    report: dict = {
        "loeuf_cutoff": LOEUF_CUTOFF,
        "min_expected": MIN_EXPECTED,
        "genes": {name: int(keep.sum()) for name, keep in kept.items()},
        "constrained_among_cell_essential": float(
            constrained[usable & ~kept["without_cell_essential"]].mean()
        ),
        "results": {name: {} for name in kept},
    }
    for set_name, members in ({name: [name] for name in names} | {"all": names}).items():
        parts = blocks(context, features, members)
        parts["expression"] = np.column_stack([parts["expression"], expected])
        links = [features[n].links for n in members]
        for name, keep in kept.items():
            predictions = held_out(parts, links, constrained, context.groups, keep)
            report["results"][name][set_name] = compare(constrained[keep], predictions, runs)
        print(
            set_name,
            json.dumps({name: report["results"][name][set_name] for name in kept}),
            flush=True,
        )
    output.mkdir(parents=True, exist_ok=True)
    (output / "constraint.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
