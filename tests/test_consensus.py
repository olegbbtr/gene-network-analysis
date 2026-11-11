import unittest

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.consensus import link_counts, listed_pairs, tally, together


class ConsensusTests(unittest.TestCase):
    def test_links_are_counted_once_per_network_with_their_sign(self):
        first = csr_array(np.array([[0, 2.0, 0], [3.0, 0, -1.0], [0, 0, 0]]))
        second = csr_array(np.array([[0, -2.0, 0], [0.5, 0, 0], [0, -4.0, 0]]))
        counts = link_counts([first, second])
        self.assertEqual(counts["any"][0, 1], 2)
        self.assertEqual(counts["positive"][0, 1], 1)  # second network: -2 + 0.5 is negative
        self.assertEqual(counts["negative"][1, 0], 1)
        self.assertEqual(counts["negative"][1, 2], 2)
        np.testing.assert_array_equal(counts["any"], counts["any"].T)

    def test_sets_above_the_size_limit_do_not_join_genes(self):
        genes = [f"g{i}" for i in range(300)]
        joined = together(genes, {"small": {"g0", "g1"}, "huge": set(genes), "one": {"g5"}})
        self.assertTrue(joined[0, 1] and joined[1, 0])
        self.assertFalse(joined[0, 2] or joined[5, 6])

    def test_tally_gives_the_share_of_known_pairs_per_count(self):
        count = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], dtype=np.int8)
        known = {"listed": listed_pairs(["a", "b", "c"], [("b", "a"), ("x", "a")])}
        result = tally(count, known, np.triu(np.ones((3, 3), dtype=bool), 1), 3)
        self.assertEqual(result["base"], {"pairs": 3, "listed": 1 / 3})
        self.assertEqual(result["by_networks"][0], {"networks": 0, "pairs": 1, "listed": 0.0})
        self.assertEqual(result["by_networks"][1], {"networks": 1, "pairs": 2, "listed": 0.5})
        self.assertIsNone(result["by_networks"][2]["listed"])


if __name__ == "__main__":
    unittest.main()
