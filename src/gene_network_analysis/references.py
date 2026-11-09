"""Curated TF->target references, converted to Ensembl gene IDs through HGNC."""

import csv
from collections.abc import Iterable, Iterator
from pathlib import Path


def _values(field: str) -> list[str]:
    return [value for value in field.strip('"').split("|") if value]


class GeneMap:
    def __init__(self, hgnc: Path):
        self.uniprot: dict[str, set[str]] = {}
        self.symbols: tuple[dict[str, set[str]], ...] = ({}, {}, {})
        with hgnc.open(newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                gene = row["ensembl_gene_id"]
                if not gene:
                    continue
                for protein in _values(row["uniprot_ids"]):
                    self.uniprot.setdefault(protein, set()).add(gene)
                for table, field in zip(self.symbols, ("symbol", "prev_symbol", "alias_symbol")):
                    for symbol in _values(row[field]):
                        table.setdefault(symbol.upper(), set()).add(gene)

    def by_uniprot(self, protein: str) -> str | None:
        genes = self.uniprot.get(protein, set())
        return next(iter(genes)) if len(genes) == 1 else None

    def by_symbol(self, symbol: str) -> str | None:
        """Approved symbol first, then a previous symbol, then an alias, if unambiguous."""
        for table in self.symbols:
            genes = table.get(symbol.upper(), set())
            if genes:
                return next(iter(genes)) if len(genes) == 1 else None
        return None


def _sign(stimulation: bool, inhibition: bool) -> int:
    return 1 if stimulation and not inhibition else -1 if inhibition and not stimulation else 0


def omnipath_pairs(path: Path, genes: GeneMap, levels: set[str] | None = None) -> Iterator:
    """OmniPath interactions; complexes count as each member regulating the target."""
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if levels is not None and not levels & set(row["dorothea_level"].split(";")):
                continue
            source = row["source"]
            members = (
                source.split(":", 1)[1].split("_") if source.startswith("COMPLEX:") else [source]
            )
            sign = _sign(
                row["consensus_stimulation"] == "True", row["consensus_inhibition"] == "True"
            )
            target = genes.by_uniprot(row["target"])
            for member in members:
                yield genes.by_uniprot(member), target, sign


def trrust_pairs(path: Path, genes: GeneMap) -> Iterator:
    modes = {"Activation": 1, "Repression": -1, "Unknown": 0}
    with path.open() as stream:
        for line in stream:
            source, target, mode, _ = line.rstrip("\n").split("\t")
            yield genes.by_symbol(source), genes.by_symbol(target), modes[mode]


def combine(pairs: Iterable) -> tuple[dict[tuple[str, str], int], dict]:
    """One sign per pair: the known sign if all evidence agrees, else 0."""
    signs: dict[tuple[str, str], set[int]] = {}
    total = unmapped = 0
    for source, target, sign in pairs:
        total += 1
        if source is None or target is None or source == target:
            unmapped += 1
            continue
        signs.setdefault((source, target), set()).add(sign)
    result = {}
    for pair, values in signs.items():
        known = values - {0}
        result[pair] = known.pop() if len(known) == 1 else 0
    return result, {"input_pairs": total, "dropped": unmapped, "pairs": len(result)}


def write_edges(path: Path, pairs: dict[tuple[str, str], int]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["source_gene", "target_gene", "sign"])
        writer.writerows([*pair, sign] for pair, sign in sorted(pairs.items()))


def read_edges(path: Path) -> list[tuple[str, str, int]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return [(row["source_gene"], row["target_gene"], int(row["sign"])) for row in reader]


GENE_COLUMNS = ("ensembl_gene_id", "symbol", "location", "locus_group", "gene_group")


def write_gene_table(hgnc: Path, output: Path) -> int:
    """HGNC genes with an Ensembl ID, reduced to the columns the diagnostics use."""
    with hgnc.open(newline="") as stream:
        rows = sorted(
            (row for row in csv.DictReader(stream, delimiter="\t") if row["ensembl_gene_id"]),
            key=lambda row: row["ensembl_gene_id"],
        )
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(GENE_COLUMNS)
        writer.writerows([row[column] for column in GENE_COLUMNS] for row in rows)
    return len(rows)


def read_gene_table(path: Path) -> dict[str, dict]:
    with path.open(newline="") as stream:
        return {
            row["ensembl_gene_id"]: {
                **row,
                "gene_group": [g for g in row["gene_group"].split("|") if g],
            }
            for row in csv.DictReader(stream, delimiter="\t")
        }
