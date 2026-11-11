import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from sklearn.metrics import roc_auc_score

from gene_network_analysis.graph import Graph
from gene_network_analysis.prediction import (
    aligned,
    cross_validated,
    descriptors,
    evaluate,
    labelled_feature,
    neighbour_share,
    network_features,
    pagerank_features,
    personalised,
    read_gene_list,
    undirected,
)


def planted(rng: np.random.Generator, genes: int = 600):
    """Listed genes form a dense positive module; expression carries no signal."""
    labels = np.zeros(genes, dtype=bool)
    labels[rng.choice(genes, 60, replace=False)] = True
    members = np.flatnonzero(labels)
    rows, columns = [], []
    for source in members:
        for target in rng.choice(members, 15, replace=False):
            if source != target:
                rows.append(source)
                columns.append(target)
    for _ in range(3000):
        rows.append(rng.integers(genes))
        columns.append(rng.integers(genes))
    matrix = csr_array((np.full(len(rows), 3.0), (rows, columns)), shape=(genes, genes))
    matrix.setdiag(0)
    matrix.eliminate_zeros()
    names = [f"ENSG{i:011d}" for i in range(genes)]
    return names, Graph(tuple(names), matrix), labels


class PredictionTests(unittest.TestCase):
    def test_aligned_reorders_and_gives_missing_genes_no_edges(self):
        graph = Graph(("b", "a"), csr_array(np.array([[0, 2.0], [0, 0]])))
        matrix = aligned(graph, ["a", "b", "c"])
        self.assertEqual(matrix[1, 0], 2.0)
        self.assertEqual(matrix.nnz, 1)

    def test_descriptors_count_signs_clustering_and_arm_share(self):
        matrix = csr_array(np.array([[0, 3.0, -3.0], [0, 0, 3.0], [0, 0, 0]]))
        features = descriptors(matrix, np.array([0, 0, 1]))
        self.assertAlmostEqual(features["log_positive_out"][0], np.log1p(1))
        self.assertAlmostEqual(features["log_negative_out"][0], np.log1p(1))
        self.assertEqual(features["clustering"][0], 1.0)
        self.assertEqual(features["arm_share"][0], 0.5)
        self.assertEqual(features["net_flow"][0], 6.0)

    def test_neighbour_share_ignores_labels_hidden_from_training(self):
        links = undirected(csr_array(np.array([[0, 1.0, 1.0], [0, 0, 0], [0, 0, 0]])))
        share = neighbour_share([links], np.array([False, True, False]))
        self.assertEqual(share[0, 0], 0.5)
        self.assertEqual(share[1, 0], 0.0)

    def test_neighbour_feature_does_not_see_test_labels(self):
        rng = np.random.default_rng(0)
        labels = rng.random(400) < 0.2
        links = undirected(csr_array(np.diag(np.ones(399), 1)))
        groups = np.arange(400) % 5
        blocks = {"expression": rng.normal(size=(400, 1))}
        # Every gene's only neighbours are in other folds, and labels are random:
        # a leak of test labels would show up as AUC far from 0.5.
        predictions = cross_validated(blocks, ("neighbours",), labels, groups, [links])
        self.assertLess(abs(roc_auc_score(labels, predictions) - 0.5), 0.1)

    def test_personalised_feature_does_not_see_test_labels(self):
        rng = np.random.default_rng(2)
        labels = rng.random(400) < 0.2
        links = undirected(csr_array(np.diag(np.ones(399), 1)))
        groups = np.arange(400) % 5
        predictions = cross_validated({}, ("personalised",), labels, groups, [links])
        self.assertLess(abs(roc_auc_score(labels, predictions) - 0.5), 0.1)

    def test_a_training_genes_personalised_score_ignores_its_own_label(self):
        links = undirected(csr_array(np.diag(np.ones(39), 1)))
        groups = np.arange(40) % 8
        train = np.flatnonzero(groups != 0)
        labels = np.zeros(40, dtype=bool)
        labels[[3, 10, 20]] = True
        flipped = labels.copy()
        flipped[3] = False
        scores = [
            labelled_feature("personalised", personalised, [links], y, groups, train, 5)[3, 0]
            for y in (labels, flipped)
        ]
        self.assertEqual(scores[0], scores[1])

    def test_pagerank_features_favour_the_end_of_a_chain_and_reverse_its_start(self):
        chain = csr_array(np.diag(np.ones(3), 1))  # 0 -> 1 -> 2 -> 3
        features = pagerank_features(chain)
        self.assertEqual(features[:, 0].argmax(), 3)
        self.assertEqual(features[:, 1].argmax(), 0)

    def test_a_planted_module_is_found_by_network_but_not_expression(self):
        rng = np.random.default_rng(1)
        names, graph, labels = planted(rng)
        arms = np.zeros(len(names), dtype=int)
        features = network_features(graph, names, arms)
        groups = np.arange(len(names)) % 10
        result = evaluate(rng.normal(size=(len(names), 2)), [features], labels, groups, runs=50)
        self.assertLess(result["auc"]["expression"], 0.65)
        self.assertGreater(result["auc"]["shape"], 0.9)
        self.assertGreater(result["gains"]["network_beyond_expression"]["ci95"][0], 0)

    def test_gene_lists_map_symbols_and_skip_the_header(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "list.txt"
            path.write_text("GENE\tHGNC_ID\nAARS\tHGNC:20\nNOPE\t-\n")
            self.assertEqual(read_gene_list(path, {"AARS": "ENSG1", "GENE": "X"}), {"ENSG1"})


if __name__ == "__main__":
    unittest.main()
