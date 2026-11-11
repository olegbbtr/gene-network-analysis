"""Graphs inferred from public human expression data: TCGA tumours (UCSC Xena) and GTEx
healthy tissues, paired by tissue.

Each graph is inferred as in `simulate.py`, with every gene a candidate sender, so no TF
knowledge enters it. All graphs use one gene set: the protein-coding genes present in
every dataset with the highest mean expression quantile across them.
"""

import gzip
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .graph import write_graph
from .references import read_gene_table
from .simulate import infer


@dataclass(frozen=True)
class Matrix:
    genes: list[str]
    samples: list[str]
    values: np.ndarray  # genes x samples


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_xena(path: Path) -> Matrix:
    """A UCSC Xena genomic matrix: gene symbol rows, sample columns."""
    with _open(path) as stream:
        samples = stream.readline().rstrip("\n").split("\t")[1:]
        genes, rows = [], []
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[0])
            rows.append(np.array(fields[1:], dtype=float))
    return Matrix(genes, samples, np.vstack(rows))


def read_gct(path: Path) -> Matrix:
    """A GTEx GCT file of TPMs, as log2(TPM + 1), rows keyed by unversioned Ensembl ID."""
    with _open(path) as stream:
        stream.readline()
        stream.readline()
        header = stream.readline().rstrip("\n").split("\t")
        first = header.index("Description") + 1
        genes, rows = [], []
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            gene = fields[first - 2]
            if gene.endswith("_PAR_Y"):
                continue
            genes.append(gene.split(".", 1)[0])
            rows.append(np.array(fields[first:], dtype=float))
    return Matrix(genes, header[first:], np.log2(np.vstack(rows) + 1))


def _of_type(matrix: Matrix, code: str) -> Matrix:
    keep = [i for i, s in enumerate(matrix.samples) if s.split("-")[3:4] == [code]]
    return Matrix(matrix.genes, [matrix.samples[i] for i in keep], matrix.values[:, keep])


def primary_tumours(matrix: Matrix) -> Matrix:
    """TCGA primary tumour samples (sample type 01)."""
    return _of_type(matrix, "01")


def adjacent_normals(matrix: Matrix) -> Matrix:
    """TCGA samples of normal tissue next to the tumour (sample type 11)."""
    return _of_type(matrix, "11")


def to_ensembl(matrix: Matrix, symbols: dict[str, str]) -> Matrix:
    """Rename symbol rows to Ensembl IDs; unknown and duplicated symbols are dropped."""
    ids = [symbols.get(g) for g in matrix.genes]
    counts: dict[str, int] = {}
    for gene in ids:
        if gene:
            counts[gene] = counts.get(gene, 0) + 1
    keep = [i for i, gene in enumerate(ids) if gene and counts[gene] == 1]
    return Matrix([ids[i] for i in keep], matrix.samples, matrix.values[keep])


def select_genes(matrices: list[Matrix], allowed: set[str], size: int) -> list[str]:
    """The `size` allowed genes, varying in every dataset, with the highest mean quantile."""
    common = sorted(allowed.intersection(*(m.genes for m in matrices)))
    varying = np.ones(len(common), dtype=bool)
    quantiles = np.zeros(len(common))
    for matrix in matrices:
        values = rows(matrix, common)
        varying &= values.std(axis=1) > 0
        quantiles += values.mean(axis=1).argsort().argsort() / len(common)
    common = [g for g, keep in zip(common, varying, strict=True) if keep]
    quantiles = quantiles[varying]
    order = np.lexsort((np.array(common), -quantiles))
    return sorted(common[i] for i in order[:size])


def rows(matrix: Matrix, genes: list[str]) -> np.ndarray:
    index = {g: i for i, g in enumerate(matrix.genes)}
    return matrix.values[[index[g] for g in genes]]


def standardise(values: np.ndarray) -> np.ndarray:
    """Genes x samples to samples x genes, each gene with mean 0 and SD 1."""
    if (values.std(axis=1) == 0).any():
        raise ValueError("a selected gene is constant across samples")
    return zscores(values)


def zscores(values: np.ndarray) -> np.ndarray:
    """Genes x samples to samples x genes, mean 0 and SD 1; constant genes become 0."""
    data = values.T
    spread = data.std(axis=0)
    centred = data - data.mean(axis=0)
    return np.divide(centred, spread, out=np.zeros_like(centred), where=spread > 0)


def remove_copy_number(expression: np.ndarray, copy_number: np.ndarray) -> np.ndarray:
    """Residuals of each gene's expression after a linear fit on its own copy number.

    Both are samples x genes; genes without copy-number variation are left unchanged.
    """
    centred = copy_number - copy_number.mean(axis=0)
    variance = (centred**2).sum(axis=0)
    slope = np.divide(
        (centred * expression).sum(axis=0),
        variance,
        out=np.zeros(len(variance)),
        where=variance > 0,
    )
    return expression - centred * slope


@dataclass(frozen=True)
class Dataset:
    name: str
    tissue: str
    state: str  # "tumour" or "healthy"
    path: str  # relative to the references directory
    copy_number: str | None = None


def _tcga(code: str, tissue: str, copy_number: bool = False) -> Dataset:
    folder = f"tcga-{code}"
    return Dataset(
        folder,
        tissue,
        "tumour",
        f"{folder}/HiSeqV2.gz",
        f"{folder}/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz" if copy_number else None,
    )


def _gtex(tissue: str, slug: str) -> Dataset:
    return Dataset(
        f"gtex-{tissue}",
        tissue,
        "healthy",
        f"gtex-v8-{tissue}/gene_tpm_2017-06-05_v8_{slug}.gct.gz",
    )


DATASETS = (
    _tcga("luad", "lung", copy_number=True),
    _tcga("lusc", "lung"),
    _tcga("brca", "breast"),
    _tcga("coad", "colon"),
    _tcga("thca", "thyroid"),
    _tcga("prad", "prostate"),
    _tcga("stad", "stomach"),
    _tcga("lihc", "liver"),
    _gtex("lung", "lung"),
    _gtex("breast", "breast_mammary_tissue"),
    _gtex("colon", "colon_transverse"),
    _gtex("thyroid", "thyroid"),
    _gtex("prostate", "prostate"),
    _gtex("stomach", "stomach"),
    _gtex("liver", "liver"),
)


def read_dataset(dataset: Dataset, root: Path, symbols: dict[str, str]) -> Matrix:
    if dataset.state == "healthy":
        return read_gct(root / dataset.path)
    return to_ensembl(primary_tumours(read_xena(root / dataset.path)), symbols)


def read_samples(
    dataset: Dataset, root: Path, symbols: dict[str, str]
) -> tuple[Matrix, Matrix | None]:
    """Expression and, if the dataset has it, copy number on the samples both measure."""
    matrix = read_dataset(dataset, root, symbols)
    if not dataset.copy_number:
        return matrix, None
    copy_number = to_ensembl(primary_tumours(read_xena(root / dataset.copy_number)), symbols)
    return same_samples(matrix, copy_number)


def same_samples(first: Matrix, second: Matrix) -> tuple[Matrix, Matrix]:
    shared = sorted(set(first.samples) & set(second.samples))
    return tuple(
        Matrix(m.genes, shared, m.values[:, [m.samples.index(s) for s in shared]])
        for m in (first, second)
    )


def write_dataset(
    output: Path,
    genes: list[str],
    data: np.ndarray,
    threshold: float,
    inference: str,
    constant_genes: bool = False,
) -> dict:
    """With `constant_genes`, a gene that does not vary is allowed and gets no links."""
    samples = zscores(data) if constant_genes else standardise(data)
    with np.errstate(invalid="ignore" if constant_genes else "warn"):
        graph = infer(samples, np.arange(len(genes)), threshold, inference)
    output.mkdir(parents=True, exist_ok=True)
    write_graph(output / "graph.txt", graph, genes)
    return {
        "genes": len(genes),
        "samples": int(data.shape[1]),
        "edges": int(graph.nnz),
        "negative_edges": int((graph.data < 0).sum()),
    }


def build(
    root: Path,
    gene_table: Path,
    output: Path,
    *,
    datasets: tuple[Dataset, ...] = DATASETS,
    size: int = 12000,
    threshold: float = 2.35,
    inference: str = "correlation",
) -> dict:
    """One graph per dataset on one shared gene set, plus copy-number-removed graphs.

    Also writes `expression.tsv`: each gene's mean and SD in every dataset.
    """
    table = read_gene_table(gene_table)
    symbols = {row["symbol"]: gene for gene, row in table.items()}
    allowed = {g for g, row in table.items() if row["locus_group"] == "protein-coding gene"}
    matrices, copy_numbers = {}, {}
    for dataset in datasets:
        matrix, copy_number = read_samples(dataset, root, symbols)
        if copy_number is not None:
            copy_numbers[dataset.name] = copy_number
            allowed &= set(copy_number.genes)
        matrices[dataset.name] = matrix
    genes = select_genes(list(matrices.values()), allowed, size)
    stats: dict = {}
    summary = {}
    for dataset in datasets:
        values = rows(matrices[dataset.name], genes)
        summary[f"{dataset.name}_mean"] = values.mean(axis=1)
        summary[f"{dataset.name}_sd"] = values.std(axis=1)
        stats[dataset.name] = {"tissue": dataset.tissue, "state": dataset.state} | write_dataset(
            output / dataset.name, genes, values, threshold, inference
        )
        if dataset.name in copy_numbers:
            adjusted = remove_copy_number(values.T, rows(copy_numbers[dataset.name], genes).T).T
            name = f"{dataset.name}-copy-number-removed"
            stats[name] = {"tissue": dataset.tissue, "state": dataset.state} | write_dataset(
                output / name, genes, adjusted, threshold, inference
            )
    write_expression(output / "expression.tsv", genes, summary)
    stats["settings"] = {"genes": size, "threshold": threshold, "inference": inference}
    (output / "stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    return stats


def build_normals(
    root: Path,
    gene_table: Path,
    networks: Path,
    output: Path,
    *,
    datasets: tuple[Dataset, ...] = DATASETS,
) -> dict:
    """Graphs from the adjacent normal tissue of every tumour dataset.

    They use the genes and settings of the existing `networks` and are written to a new
    folder that also links those networks, so all of them can be compared. TCGA processed
    these samples together with its tumours, unlike the GTEx healthy tissues.
    """
    symbols = {row["symbol"]: gene for gene, row in read_gene_table(gene_table).items()}
    genes, summary = read_expression(networks / "expression.tsv")
    stats = json.loads((networks / "stats.json").read_text())
    settings = stats.pop("settings")
    output.mkdir(parents=True, exist_ok=True)
    for name in stats:
        if not (output / name).exists():
            (output / name).symlink_to((networks / name).resolve())
    for dataset in datasets:
        if dataset.state != "tumour":
            continue
        values = rows(to_ensembl(adjacent_normals(read_xena(root / dataset.path)), symbols), genes)
        name = f"{dataset.name}-normal"
        summary[f"{name}_mean"] = values.mean(axis=1)
        summary[f"{name}_sd"] = values.std(axis=1)
        stats[name] = {"tissue": dataset.tissue, "state": "normal"} | write_dataset(
            output / name,
            genes,
            values,
            settings["threshold"],
            settings["inference"],
            constant_genes=True,
        )
    write_expression(output / "expression.tsv", genes, summary)
    stats["settings"] = settings
    (output / "stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    return stats


def write_expression(path: Path, genes: list[str], columns: dict[str, np.ndarray]) -> None:
    with path.open("w") as stream:
        stream.write("\t".join(["gene", *columns]) + "\n")
        for i, gene in enumerate(genes):
            stream.write("\t".join([gene, *(f"{c[i]:.6f}" for c in columns.values())]) + "\n")


def read_expression(path: Path) -> tuple[list[str], dict[str, np.ndarray]]:
    with path.open() as stream:
        header = stream.readline().rstrip("\n").split("\t")[1:]
        genes, values = [], []
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[0])
            values.append([float(v) for v in fields[1:]])
    array = np.array(values)
    return genes, {name: array[:, i] for i, name in enumerate(header)}
