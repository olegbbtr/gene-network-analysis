import unittest

import numpy as np
from scipy.sparse import csr_array
from sklearn.metrics import roc_auc_score

from gene_network_analysis.pagerank_discovery import predict_variants
from gene_network_analysis.robustness import held_out


class PageRankDiscoveryTests(unittest.TestCase):
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
