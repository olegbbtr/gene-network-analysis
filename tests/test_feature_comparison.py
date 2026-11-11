import unittest

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.feature_comparison import MODELS
from gene_network_analysis.pagerank_discovery import predict_variants
from gene_network_analysis.robustness import held_out


class FeatureComparisonTests(unittest.TestCase):
    def test_network_only_ignores_expression_and_original_models_are_preserved(self):
        rng = np.random.default_rng(8)
        size = 400
        labels = rng.random(size) < 0.2
        groups = np.arange(size) % 8
        parts = {
            key: rng.normal(size=(size, 2))
            for key in (
                "expression",
                "degree",
                "shape",
                "pagerank",
                "pagerank_rank",
                "graph_theory",
                "weights",
            )
        }
        links = [csr_array(np.diag(np.ones(size - 1), 1))]
        models = {key: MODELS[key] for key in ("expression", "network", "expression+network")}
        first = predict_variants(parts, links, labels, groups, models)
        original = held_out(parts, links, labels, groups)
        for key, values in original.items():
            np.testing.assert_array_equal(first[key], values)
        changed = parts | {"expression": rng.normal(size=(size, 2))}
        second = predict_variants(changed, links, labels, groups, models)
        np.testing.assert_array_equal(first["network"], second["network"])
        self.assertGreater(np.max(np.abs(first["expression"] - second["expression"])), 0.01)


if __name__ == "__main__":
    unittest.main()
