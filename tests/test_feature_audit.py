import unittest

import numpy as np

from gene_network_analysis.feature_audit import shuffled_within_groups


class FeatureAuditTests(unittest.TestCase):
    def test_shuffle_preserves_chromosome_counts_and_has_a_reproducible_seed(self):
        groups = np.repeat(["1", "2", "3"], 30)
        labels = np.arange(90) % 3 == 0
        shuffled = shuffled_within_groups(labels, groups, 0)
        np.testing.assert_array_equal(shuffled, shuffled_within_groups(labels, groups, 0))
        self.assertFalse(np.array_equal(labels, shuffled))
        self.assertFalse(np.array_equal(shuffled, shuffled_within_groups(labels, groups, 1)))
        for group in np.unique(groups):
            self.assertEqual(labels[groups == group].sum(), shuffled[groups == group].sum())


if __name__ == "__main__":
    unittest.main()
