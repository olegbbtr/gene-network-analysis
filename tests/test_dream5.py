import tempfile
import unittest
from pathlib import Path

import numpy as np

from gene_network_analysis.dream5 import prepare, read_expression
from gene_network_analysis.evaluation import read_reference, resolve_labels
from gene_network_analysis.graph import read_graph


def write_dataset(directory: Path) -> dict[str, Path]:
    rng = np.random.default_rng(0)
    data = rng.normal(size=(300, 30))
    data[:, 1] += 2 * data[:, 0]
    data[:, 2] -= 2 * data[:, 0]
    data[:, 29] = 5.0
    paths = {name: directory / name for name in ("expression", "tfs", "gold", "signed")}
    lines = ["\t".join(f"G{i + 1}" for i in range(30))]
    lines += ["\t".join(f"{v:.6f}" for v in row) for row in data]
    paths["expression"].write_text("\n".join(lines) + "\n")
    paths["tfs"].write_text("G1\nG4\n")
    paths["gold"].write_text("G1\tG2\t1\nG1\tG3\t1\nG4\tG5\t0\n")
    paths["signed"].write_text("G1\tG2\t+\nG1\tG3\t-\n")
    return paths


class Dream5Tests(unittest.TestCase):
    def test_expression_is_standardised_and_constant_genes_dropped(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = write_dataset(Path(directory))
            genes, data = read_expression(paths["expression"])
            self.assertEqual(len(genes), 29)
            self.assertNotIn("G30", genes)
            np.testing.assert_allclose(data.std(axis=0), 1.0)

    def test_convert_writes_inputs_the_pipelines_read(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = write_dataset(Path(directory))
            output = Path(directory) / "out"
            stats = prepare(
                "1",
                paths["expression"],
                paths["tfs"],
                paths["gold"],
                output,
                signed=paths["signed"],
            )
            graph = read_graph(output / "graph.txt")
            self.assertEqual(graph.adjacency.nnz, stats["inferred_edges"])
            index = {node: i for i, node in enumerate(graph.nodes)}
            source, up, down = (index[f"ENSGD5N1G{i}"] for i in (1, 2, 3))
            self.assertGreater(graph.adjacency[source, up], 0)
            self.assertLess(graph.adjacency[source, down], 0)
            labels, known = resolve_labels(graph.nodes, read_reference(output / "tfs.txt"))
            self.assertTrue(known.all())
            self.assertEqual(labels[source], True)
            self.assertEqual(stats["regulators_with_true_targets"], 1)
            truth = (output / "true-edges.tsv").read_text().splitlines()
            self.assertEqual(truth[1:], ["ENSGD5N1G1\tENSGD5N1G2\t1", "ENSGD5N1G1\tENSGD5N1G3\t-1"])


if __name__ == "__main__":
    unittest.main()
