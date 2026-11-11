import unittest

import numpy as np
from scipy.sparse import csr_array
from sklearn.metrics import roc_auc_score

from gene_network_analysis.graph_features import rank_positions
from gene_network_analysis.pagerank_discovery import JOINT_MODELS, predict_variants
from gene_network_analysis.robustness import held_out


class PageRankDiscoveryTests(unittest.TestCase):
    def test_joint_encoding_can_capture_a_pattern_missed_by_either_linear_feature(self):
        size = 400
        position = (np.arange(size) + 0.5) / size
        known = (position > 0.3) & (position < 0.7)
        parts = {name: np.zeros((size, 1)) for name in ("expression", "degree", "shape")}
        parts["pagerank"] = np.column_stack([position**3, np.zeros(size)])
        parts["pagerank_rank"] = rank_positions(parts["pagerank"])
        scores = predict_variants(
            parts, [csr_array((size, size))], known, np.arange(size) % 8, JOINT_MODELS
        )
        joint = roc_auc_score(known, scores["expression+network+pagerank-both"])
        self.assertGreater(joint, 0.95)
        for key in ("expression+network+pagerank", "expression+network+pagerank-rank"):
            self.assertLess(roc_auc_score(known, scores[key]), 0.65)

    def test_pagerank_can_add_signal_without_changing_the_original_model(self):
        rng = np.random.default_rng(7)
        size = 400
        known = rng.random(size) < 0.3
        groups = np.arange(size) % 8
        parts = {name: rng.normal(size=(size, 2)) for name in ("expression", "degree", "shape")}
        parts["pagerank"] = np.column_stack(
            [known + rng.normal(scale=0.1, size=size), rng.normal(size=size)]
        )
        links = [csr_array(np.diag(np.ones(size - 1), 1))]
        scores = predict_variants(parts, links, known, groups)
        original = held_out(parts, links, known, groups)
        for name, expected in original.items():
            np.testing.assert_array_equal(scores[name], expected)
        self.assertGreater(roc_auc_score(known, scores["expression+network+pagerank"]), 0.95)
        self.assertLess(roc_auc_score(known, scores["expression+network"]), 0.65)

    def test_personalised_variant_has_no_signal_when_labels_are_random(self):
        rng = np.random.default_rng(2)
        size = 400
        known = rng.random(size) < 0.2
        groups = np.arange(size) % 8
        parts = {
            name: rng.normal(size=(size, 2))
            for name in ("expression", "degree", "shape", "pagerank")
        }
        links = [csr_array(np.diag(np.ones(size - 1), 1))]
        scores = predict_variants(parts, links, known, groups)
        self.assertLess(
            abs(roc_auc_score(known, scores["expression+network+personalised"]) - 0.5), 0.1
        )


if __name__ == "__main__":
    unittest.main()
