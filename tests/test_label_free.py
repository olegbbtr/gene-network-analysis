import unittest

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.label_free import (
    ALL,
    COMBINATIONS,
    GROUPS,
    LABELS,
    MODELS,
    WITHOUT,
    coherence,
    consistent_partners,
    leave_one_out,
    removal,
    shapley,
)
from gene_network_analysis.prediction import LABELLED


class LabelFreeTests(unittest.TestCase):
    def test_no_model_uses_a_feature_computed_from_labels(self):
        for name, parts in MODELS.items():
            self.assertFalse(set(parts) & set(LABELLED), name)

    def test_shown_and_leave_one_out_models_are_fitted(self):
        self.assertEqual(set(MODELS), set(LABELS) | set(WITHOUT.values()))
        self.assertEqual(set(MODELS[ALL]), {p for parts in GROUPS.values() for p in parts})
        self.assertEqual(WITHOUT["graph"], "gene+biological")
        self.assertNotIn("cutoff_counts", MODELS[WITHOUT["graph"]])

    def test_removal_is_the_auc_lost_without_a_block(self):
        labels = np.array([0, 0, 1, 1] * 25)
        perfect = labels.astype(float)
        noise = np.random.default_rng(0).random(len(labels))
        predictions = {ALL: perfect} | {m: perfect for m in WITHOUT.values()}
        predictions[WITHOUT["graph"]] = noise
        losses = removal(labels, predictions)
        self.assertAlmostEqual(losses["gene"]["auc_loss"], 0.0)
        self.assertAlmostEqual(losses["biological"]["auc_loss"], 0.0)
        self.assertGreater(losses["graph"]["auc_loss"], 0.3)

    def test_shapley_splits_shared_information_and_sums_to_the_gain(self):
        # gene and local carry the same 0.2; global adds 0.1 on its own; biological nothing
        def auc(name):
            groups = set(name.split("+"))
            return 0.5 + 0.2 * bool(groups & {"gene", "local"}) + 0.1 * ("global" in groups)

        values = shapley({name: auc(name) for name in COMBINATIONS})
        self.assertAlmostEqual(sum(values.values()), auc(ALL) - 0.5)
        self.assertAlmostEqual(values["gene"], 0.1)
        self.assertAlmostEqual(values["local"], 0.1)
        self.assertAlmostEqual(values["global"], 0.1)
        self.assertAlmostEqual(values["biological"], 0.0)

    def test_leave_one_out_misses_information_two_groups_share(self):
        def auc(name):
            groups = set(name.split("+"))
            return 0.5 + 0.2 * bool(groups & {"gene", "local"}) + 0.1 * ("global" in groups)

        losses = leave_one_out({name: auc(name) for name in COMBINATIONS})
        expected_losses = {"gene": 0.0, "local": 0.0, "global": 0.1, "biological": 0.0}
        for group, expected in expected_losses.items():
            self.assertAlmostEqual(losses[group], expected)

    def test_consistent_partners_need_a_positive_link_in_every_network(self):
        first = csr_array(np.array([[0, 1.0, 1.0], [0, 0, 0], [0, 0, 0]]))
        second = csr_array(np.array([[0, 2.0, -1.0], [0, 0, 0], [0, 0, 0]]))
        partners = consistent_partners([first, second]).toarray()
        np.testing.assert_array_equal(partners, [[0, 1, 0], [1, 0, 0], [0, 0, 0]])

    def test_coherence_counts_neighbours_that_share_a_pathway(self):
        adjacency = csr_array(np.array([[0, 1.0, 1.0], [0, 0, 0], [0, 0, 0]]))
        links = csr_array(np.array([[0, 1.0, 1.0], [1.0, 0, 0], [1.0, 0, 0]]))
        shared = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 1]], dtype=bool)
        values = coherence([links], [adjacency], shared, np.array([1.0, 1.0, 0.0]))
        np.testing.assert_allclose(values[:, 0], [0.5, 1.0, 0.0])
        np.testing.assert_allclose(values[:, 1], np.log1p([1, 1, 0]))
        np.testing.assert_allclose(values[:, 2], [1.0, 1.0, 0.0])


if __name__ == "__main__":
    unittest.main()
