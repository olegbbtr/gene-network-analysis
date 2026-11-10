import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.graph import Graph
from gene_network_analysis.structure import (
    MODULES,
    balance,
    locality,
    modules,
    signed_pairs,
    undirected_signs,
)


def graph_from(edges: list[tuple[int, int, float]], size: int) -> csr_array:
    rows, columns, weights = zip(*edges, strict=True)
    return csr_array((weights, (rows, columns)), shape=(size, size))


class BalanceTests(unittest.TestCase):
    def test_triangle_balance_counts_negative_edges(self):
        positive = undirected_signs(graph_from([(0, 1, 1), (1, 2, 1), (2, 0, 1)], 3))
        self.assertEqual(balance(positive)["triangles"], 1)
        self.assertEqual(balance(positive)["balanced_share"], 1.0)
        one_negative = undirected_signs(graph_from([(0, 1, -1), (1, 2, 1), (0, 2, 1)], 3))
        self.assertEqual(balance(one_negative)["balanced_share"], 0.0)
        two_negative = undirected_signs(graph_from([(0, 1, -1), (1, 2, -1), (0, 2, 1)], 3))
        self.assertEqual(balance(two_negative)["balanced_share"], 1.0)

    def test_pairs_in_both_directions_count_once_with_the_sign_of_their_sum(self):
        signs = undirected_signs(graph_from([(0, 1, 3), (1, 0, -1), (1, 2, -2)], 3))
        self.assertEqual(signs[0, 1], 1)
        self.assertEqual(signs[1, 0], 1)
        self.assertEqual(balance(signs)["edges"], 2)

    def test_curated_pairs_with_unknown_sign_are_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "edges.tsv"
            path.write_text("source_gene\ttarget_gene\tsign\nA\tB\t1\nB\tC\t-1\nC\tA\t0\n")
            self.assertEqual(balance(signed_pairs(path))["edges"], 2)


class LocalityTests(unittest.TestCase):
    def test_edges_within_an_arm_are_enriched_over_shuffled_locations(self):
        locations = ["1p36", "1p35", "1p34", "2q11", "2q12", "2q13", "3p21", "3p22"]
        nodes = tuple(f"ENSG{i}" for i in range(8))
        genes = {node: {"location": loc} for node, loc in zip(nodes, locations, strict=True)}
        edges = [(0, 1, 3.0), (1, 2, 3.0), (3, 4, 3.0), (4, 5, 3.0), (6, 7, -3.0), (0, 6, -3.0)]
        result = locality(Graph(nodes, graph_from(edges, 8)), genes, permutations=50)
        self.assertEqual(result["arm"]["positive"]["same_share"], 1.0)
        self.assertGreater(result["arm"]["all"]["fold"], 1.0)
        self.assertEqual(result["band"]["all"]["same_share"], 0.0)


class ModuleTests(unittest.TestCase):
    def test_hubs_sending_into_a_module_are_assigned_to_it(self):
        markers = MODULES["cell_cycle"][:3]
        symbols = markers + ["HUB"] + [f"OTHER{i}" for i in range(40)]
        nodes = tuple(f"ENSG{i}" for i in range(len(symbols)))
        genes = {node: {"symbol": s} for node, s in zip(nodes, symbols, strict=True)}
        edges = [(3, target, 2.5) for target in range(3)] + [(0, 1, 2.5), (1, 2, 2.5)]
        result = modules(Graph(nodes, graph_from(edges, len(nodes))), genes, top_k=2)
        self.assertEqual(result["top_2_out_degree"]["cell_cycle"], 2)
        self.assertEqual(result["top_20"][0]["symbol"], "HUB")
        self.assertGreater(result["modules"]["cell_cycle"]["within_density_fold"], 1)
        np.testing.assert_equal(result["modules"]["sex"]["present"], [])


if __name__ == "__main__":
    unittest.main()
