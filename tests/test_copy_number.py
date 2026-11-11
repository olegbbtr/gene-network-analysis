import unittest

import numpy as np

from gene_network_analysis.copy_number import arm_means, cross, curated, direction, outward


def driven(rng: np.random.Generator, samples: int = 400):
    """Gene 0's copy number moves its own expression and gene 1's; gene 1 drives nothing."""
    copy_number = rng.normal(size=(samples, 4))
    expression = copy_number + 0.5 * rng.normal(size=(samples, 4))
    expression[:, 1] += 2 * expression[:, 0]
    standard = (expression - expression.mean(0)) / expression.std(0)
    cn = (copy_number - copy_number.mean(0)) / copy_number.std(0)
    return cn, standard


class CopyNumberTests(unittest.TestCase):
    def test_direction_points_from_the_driver(self):
        cn, expression = driven(np.random.default_rng(0))
        effects = cross(cn, expression)
        self.assertGreater(direction(effects, np.array([0]), np.array([1]))[0], 0.3)
        self.assertLess(direction(effects, np.array([1]), np.array([0]))[0], -0.3)

    def test_curated_pairs_on_one_chromosome_are_skipped(self):
        cn, expression = driven(np.random.default_rng(1))
        effects = cross(cn, expression)
        eligible = np.ones(4, dtype=bool)
        chromosomes = np.array(["1", "2", "3", "4"])
        result = curated(effects, np.array([[0, 1]]), eligible, chromosomes, {(0, 1)})
        self.assertEqual(result["tf_to_target"]["share"], 1.0)
        self.assertEqual(result["linked_pairs"], 1)
        same = curated(effects, np.array([[0, 1]]), eligible, np.array(["1", "1", "3", "4"]), set())
        self.assertEqual(same["pairs"], 0)

    def test_arm_means_give_genes_on_one_arm_the_same_copy_number(self):
        cn = np.random.default_rng(2).normal(size=(50, 3))
        regional = arm_means(cn, np.array([0, 0, 1]))
        np.testing.assert_allclose(regional[:, 0], regional[:, 1])
        np.testing.assert_allclose(regional[:, 2], (cn[:, 2] - cn[:, 2].mean()) / cn[:, 2].std())

    def test_outward_counts_links_pointing_away_on_other_chromosomes(self):
        effects = np.zeros((4, 4), dtype=np.float32)
        effects[0, 1:] = 0.5  # gene 0's copy number moves every other gene
        links = np.array([[0, 1], [2, 0], [0, 3], [1, 2]])
        chromosomes = np.array(["1", "2", "3", "1"])  # the 0-3 link stays on chromosome 1
        away, total = outward(effects, links, np.ones(4, dtype=bool), chromosomes)
        np.testing.assert_array_equal(total, [2, 2, 2, 0])
        self.assertEqual(away[0], 2)
        self.assertEqual(away[1], 0)


if __name__ == "__main__":
    unittest.main()
