import tempfile
import unittest
from pathlib import Path

import numpy as np

from gene_network_analysis.graph import component_summary, read_graph
from gene_network_analysis.simulate import (
    Settings,
    association,
    expression,
    infer,
    simulate,
    true_network,
)


class SimulationTests(unittest.TestCase):
    def test_true_network_is_a_dag_from_regulators(self):
        settings = Settings(genes=300, regulators=30, samples=50)
        network, layer = true_network(settings, np.random.default_rng(0))
        self.assertTrue(np.all(np.diff(network.indptr)[settings.regulators :] == 0))
        self.assertEqual(component_summary(network)["strong_components"], settings.genes)
        rows, columns = network.nonzero()
        regulator_targets = columns < settings.regulators
        self.assertTrue(np.all(layer[rows[regulator_targets]] < layer[columns[regulator_targets]]))

    def test_expression_reflects_a_direct_edge(self):
        settings = Settings(genes=300, regulators=30, samples=2000)
        rng = np.random.default_rng(1)
        network, _ = true_network(settings, rng)
        data = expression(network, settings, rng)
        source, target = network.nonzero()[0][0], network.nonzero()[1][0]
        correlation = np.corrcoef(data[:, source], data[:, target])[0, 1]
        self.assertEqual(np.sign(correlation), np.sign(network[source, target]))

    def test_inference_is_directional_and_thresholded(self):
        rng = np.random.default_rng(2)
        data = rng.normal(size=(200, 40))
        data[:, 1] += data[:, 0]
        graph = infer((data - data.mean(0)) / data.std(0), np.arange(20), 2.0)
        self.assertTrue(np.all(np.abs(graph.data) > 2.0))
        self.assertEqual(np.diff(graph.indptr)[20:].sum(), 0)
        self.assertNotEqual(graph[0, 1], graph[1, 0])

    def test_partial_correlation_removes_an_indirect_link(self):
        rng = np.random.default_rng(4)
        a = rng.normal(size=5000)
        b = a + rng.normal(size=5000)
        c = b + rng.normal(size=5000)
        data = np.column_stack([a, b, c])
        data = (data - data.mean(0)) / data.std(0)
        direct = association(data, "correlation", 0.0)
        partial = association(data, "partial", 1e-6)
        self.assertGreater(abs(direct[0, 2]), 0.4)
        self.assertLess(abs(partial[0, 2]), 0.05)
        self.assertGreater(partial[1, 2], 0.4)

    def test_simulate_writes_inputs_the_pipelines_read(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            stats = simulate(Settings(genes=300, regulators=30, samples=100), output)
            graph = read_graph(output / "graph.txt")
            self.assertEqual(graph.adjacency.nnz, stats["inferred_edges"])
            self.assertEqual(len((output / "regulators.txt").read_text().split()), 30)
            self.assertTrue((output / "true-edges.tsv").read_text().startswith("source_gene"))


if __name__ == "__main__":
    unittest.main()
