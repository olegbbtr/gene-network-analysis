import unittest

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.regulons import (
    MIN_TARGETS,
    elevated,
    enrichment,
    regulons,
    sign_agreement,
    tissue_match,
)


def network(size=200):
    """TF 0 is linked to its own targets (1-20); TF 100 is linked to genes that are not its."""
    genes = [f"g{i}" for i in range(size)]
    pairs = [("g0", f"g{i}", 1 if i % 2 else -1) for i in range(1, 21)]
    pairs += [("g100", f"g{i}", 0) for i in range(101, 121)]
    pairs += [("g0", "g0", 1), ("g50", "g51", 1), ("missing", "g1", 1)]
    dense = np.zeros((size, size))
    dense[0, 1:21] = 1.0
    dense[0, 2:21:2] = -1.0
    dense[100, 150:170] = 1.0
    return genes, pairs, csr_array(dense)


class RegulonTests(unittest.TestCase):
    def test_only_tfs_with_enough_targets_are_kept(self):
        genes, pairs, _ = network()
        tfs, member, signs = regulons(pairs, genes)
        self.assertEqual(list(tfs), [0, 100])  # g50 has one target; g0 -> g0 is dropped
        self.assertEqual(member.sum(axis=1).tolist(), [20, 20])
        self.assertGreaterEqual(20, MIN_TARGETS)
        self.assertEqual(signs[0, 1], 1)
        self.assertEqual(signs[0, 2], -1)

    def test_a_tf_linked_to_its_targets_is_enriched_and_one_that_is_not_is_not(self):
        genes, pairs, adjacency = network()
        tfs, member, _ = regulons(pairs, genes)
        links = ((adjacency + adjacency.T) != 0).astype(float).tocsr()
        found = enrichment(links, tfs, member)
        self.assertEqual(found["shared"].tolist(), [20, 0])
        self.assertAlmostEqual(found["fold"][0], 199 / 20)
        self.assertLess(found["p"][0], 1e-20)
        self.assertEqual(found["p"][1], 1.0)
        self.assertEqual(found["own_rank"].tolist(), [1.0, 0.5])  # TF 100: 0 for both, a tie

    def test_a_tf_is_elevated_against_its_median_elsewhere(self):
        means = np.array([[5.0, 1.0, 1.0, 4.5], [2.0, 2.0, 2.0, 2.9], [3.0, 0.0, 3.0, 0.0]])
        self.assertEqual(
            elevated(means).tolist(),
            [[True, False, False, True], [False, False, False, False], [True, False, True, False]],
        )

    def test_tissue_match_compares_each_tf_with_itself_elsewhere(self):
        tissue = np.array([0, 0, 1, 2])  # four networks in three tissues
        high = np.array([[True, False, False], [False, True, True], [False, False, False]])
        strength = np.array([[3.0, 0.0, 9.0], [1.0, 0.0, 9.0], [0.0, 2.0, 9.0], [0.0, 4.0, 9.0]])
        significant = strength > 1.5
        found = tissue_match(strength, significant, tissue, high, runs=200)
        self.assertEqual(found["tfs"], 2)  # the third TF is elevated nowhere
        np.testing.assert_allclose(found["gain"], [2.0, 3.0])
        self.assertEqual(found["mean_difference"], 2.5)
        self.assertEqual(found["tfs_higher_in_own_tissues"], 2)
        self.assertEqual(found["significant_in_own_tissues"], 0.75)
        self.assertEqual(found["significant_elsewhere"], 0.0)

    def test_link_signs_are_compared_with_curated_signs(self):
        genes, pairs, adjacency = network()
        tfs, _, signs = regulons(pairs, genes)
        result = sign_agreement(adjacency, tfs, signs)
        self.assertEqual(result["pairs"], 20)  # TF 100's curated signs are unknown
        self.assertEqual(result["agreement"], 1.0)
        self.assertAlmostEqual(result["chance"], 0.5)


if __name__ == "__main__":
    unittest.main()
