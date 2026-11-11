import unittest

import numpy as np

from gene_network_analysis.hubs import SIGNED_POWER, UNSIGNED_POWER, connectivity, scored


class HubTests(unittest.TestCase):
    def test_connectivity_follows_the_wgcna_definitions(self):
        rng = np.random.default_rng(0)
        base = rng.normal(size=500)
        data = np.column_stack([base, base, -base, rng.normal(size=500)])
        data = (data - data.mean(axis=0)) / data.std(axis=0)
        found = connectivity(data)
        noise = np.corrcoef(data.T)[0, 3]
        expected_unsigned = 1 + 1 + abs(noise) ** UNSIGNED_POWER
        expected_signed = 1 + 0 + ((1 + noise) / 2) ** SIGNED_POWER
        self.assertAlmostEqual(found["wgcna_unsigned"][0], expected_unsigned, places=4)
        self.assertAlmostEqual(found["wgcna_signed"][0], expected_signed, places=4)
        self.assertLess(found["wgcna_unsigned"][3], 0.01)

    def test_kinds_are_scored_only_among_the_genes_they_are_defined_for(self):
        scores = {"wgcna_unsigned": np.array([5.0, 4.0, 3.0, 2.0, 1.0, 9.0])}
        label = np.array([True, True, False, False, False, True])
        known = np.array([True, True, True, True, True, False])
        result = scored(scores, {"kind": (label, known)})["kind"]
        self.assertEqual(result["auc"]["wgcna_unsigned"], 1.0)
        self.assertEqual(result["among_top_hubs"], 2)
        self.assertEqual(result["expected_among_top_hubs"], 40.0)


if __name__ == "__main__":
    unittest.main()
