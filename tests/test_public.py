import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from gene_network_analysis.graph import read_graph
from gene_network_analysis.public import (
    Dataset,
    Matrix,
    adjacent_normals,
    build,
    build_normals,
    primary_tumours,
    read_expression,
    read_gct,
    read_xena,
    remove_copy_number,
    select_genes,
    to_ensembl,
)


def write_inputs(directory: Path, rng: np.random.Generator) -> dict[str, Path]:
    genes = [f"G{i}" for i in range(30)]
    tumours = [f"TCGA-AA-{i:04d}-01" for i in range(60)]
    copy_number = rng.normal(size=(30, 60))
    expression = rng.normal(size=(30, 60)) + 3 * copy_number + np.arange(30)[:, None]
    expression[1] += 2 * expression[0]
    paths = {
        "genes": directory / "genes.tsv",
        "tcga": directory / "tcga.gz",
        "cn": directory / "cn.gz",
        "gtex": directory / "gtex.gct.gz",
    }
    rows = ["ensembl_gene_id\tsymbol\tlocation\tlocus_group\tgene_group"]
    rows += [f"ENSG{i:011d}\tG{i}\t1p36\tprotein-coding gene\t" for i in range(30)]
    paths["genes"].write_text("\n".join(rows) + "\n")
    normals = [f"TCGA-AA-{i:04d}-11" for i in range(20)]
    normal = rng.normal(size=(30, 20))
    normal[2] = normal[3] + 0.01 * rng.normal(size=20)
    normal[5] = 7.0  # constant in normal tissue
    for key, values in (("tcga", expression), ("cn", copy_number)):
        lines = ["sample\t" + "\t".join(tumours + normals)]
        lines += [
            "\t".join([g, *(f"{v:.6f}" for v in row), *(f"{v:.6f}" for v in extra)])
            for g, row, extra in zip(genes, values, normal, strict=True)
        ]
        with gzip.open(paths[key], "wt") as stream:
            stream.write("\n".join(lines) + "\n")
    gtex = np.exp2(rng.normal(size=(30, 40)) + 3)
    lines = ["#1.2", "30\t40", "id\tName\tDescription\t" + "\t".join(f"S{i}" for i in range(40))]
    lines += [
        "\t".join([str(i), f"ENSG{i:011d}.3", f"G{i}", *(f"{v:.6f}" for v in row)])
        for i, row in enumerate(gtex)
    ]
    with gzip.open(paths["gtex"], "wt") as stream:
        stream.write("\n".join(lines) + "\n")
    return paths


class PublicTests(unittest.TestCase):
    def test_readers_keep_primary_tumours_and_log_transform_tpm(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = write_inputs(Path(directory), np.random.default_rng(0))
            tumour = primary_tumours(read_xena(paths["tcga"]))
            self.assertEqual(len(tumour.samples), 60)
            self.assertEqual(tumour.values.shape, (30, 60))
            gtex = read_gct(paths["gtex"])
            self.assertEqual(gtex.genes[0], "ENSG00000000000")
            self.assertTrue((gtex.values >= 0).all())

    def test_symbols_map_to_ensembl_and_duplicates_are_dropped(self):
        matrix = Matrix(["A", "B", "B", "X"], ["s"], np.arange(4.0)[:, None])
        mapped = to_ensembl(matrix, {"A": "ENSG1", "B": "ENSG2"})
        self.assertEqual(mapped.genes, ["ENSG1"])

    def test_gene_selection_ranks_by_mean_expression_in_every_dataset(self):
        first = Matrix(
            ["a", "b", "c", "d"], ["s", "t"], np.array([[4, 5], [3, 4], [2, 3], [1, 2.0]])
        )
        second = Matrix(["a", "b", "c"], ["s", "t"], np.array([[1, 2], [3, 4], [2, 3.0]]))
        self.assertEqual(select_genes([first, second], {"a", "b", "c", "d"}, 2), ["a", "b"])

    def test_genes_constant_in_any_dataset_are_not_selected(self):
        first = Matrix(["a", "b"], ["s", "t"], np.array([[9.0, 9.0], [1.0, 2.0]]))
        second = Matrix(["a", "b"], ["s", "t"], np.array([[5.0, 6.0], [1.0, 2.0]]))
        self.assertEqual(select_genes([first, second], {"a", "b"}, 2), ["b"])

    def test_copy_number_removal_leaves_no_linear_dependence(self):
        rng = np.random.default_rng(1)
        copy_number = rng.normal(size=(200, 3))
        expression = 2 * copy_number + rng.normal(size=(200, 3))
        expression[:, 2] = rng.normal(size=200)
        copy_number[:, 2] = 0
        residual = remove_copy_number(expression, copy_number)
        for gene in range(2):
            self.assertAlmostEqual(np.corrcoef(residual[:, gene], copy_number[:, gene])[0, 1], 0)
        np.testing.assert_allclose(residual[:, 2], expression[:, 2])

    def test_build_writes_three_graphs_on_one_gene_set(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = write_inputs(Path(directory), np.random.default_rng(2))
            output = Path(directory) / "out"
            datasets = (
                Dataset("tcga-x", "lung", "tumour", "tcga.gz", "cn.gz"),
                Dataset("gtex-x", "lung", "healthy", "gtex.gct.gz"),
            )
            stats = build(Path(directory), paths["genes"], output, datasets=datasets, size=25)
            names = ("tcga-x", "tcga-x-copy-number-removed", "gtex-x")
            nodes = {name: set(read_graph(output / name / "graph.txt").nodes) for name in names}
            self.assertEqual(stats["tcga-x"]["samples"], 60)
            genes, expression = read_expression(output / "expression.tsv")
            self.assertEqual(len(genes), 25)
            self.assertEqual(
                set(expression), {"tcga-x_mean", "tcga-x_sd", "gtex-x_mean", "gtex-x_sd"}
            )
            self.assertTrue(all(len(n) <= 25 for n in nodes.values()))
            self.assertTrue(set().union(*nodes.values()) <= {f"ENSG{i:011d}" for i in range(30)})
            self.assertEqual(json.loads((output / "stats.json").read_text())["gtex-x"]["genes"], 25)

    def test_normal_tissue_graphs_join_the_existing_networks(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = write_inputs(Path(directory), np.random.default_rng(3))
            first, second = Path(directory) / "out", Path(directory) / "with-normals"
            datasets = (
                Dataset("tcga-x", "lung", "tumour", "tcga.gz", "cn.gz"),
                Dataset("gtex-x", "lung", "healthy", "gtex.gct.gz"),
            )
            build(Path(directory), paths["genes"], first, datasets=datasets, size=30)
            self.assertEqual(len(adjacent_normals(read_xena(paths["tcga"])).samples), 20)
            stats = build_normals(Path(directory), paths["genes"], first, second, datasets=datasets)
            self.assertEqual(stats["tcga-x-normal"]["samples"], 20)
            self.assertEqual(stats["tcga-x-normal"]["state"], "normal")
            self.assertEqual(stats["settings"]["threshold"], 2.35)
            genes, expression = read_expression(second / "expression.tsv")
            self.assertEqual(len(genes), 30)
            self.assertIn("tcga-x_mean", expression)
            self.assertEqual(expression["tcga-x-normal_sd"][5], 0.0)
            graph = read_graph(second / "tcga-x-normal" / "graph.txt")
            weights = graph.adjacency.toarray()
            index = {gene: i for i, gene in enumerate(graph.nodes)}
            linked = weights[index["ENSG00000000002"], index["ENSG00000000003"]]
            self.assertGreater(linked, 2.35)  # the planted pair of normal tissue
            self.assertNotIn("ENSG00000000005", index)  # the constant gene has no links
            old = read_graph(first / "gtex-x" / "graph.txt")
            self.assertEqual(read_graph(second / "gtex-x" / "graph.txt").nodes, old.nodes)


if __name__ == "__main__":
    unittest.main()
