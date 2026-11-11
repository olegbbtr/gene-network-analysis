import unittest

import numpy as np

from gene_network_analysis.partners import (
    MIN_PARTNERS,
    partner_means,
    specificity,
    summarise,
    within_ranks,
)


def screen(specific: bool, knockdowns=60, genes=400, seed=0):
    """Each knockdown has its own 20 partners. They fall under their own knockdown only
    (`specific`), or one shared set of genes falls under every knockdown."""
    rng = np.random.default_rng(seed)
    effects = rng.normal(size=(knockdowns, genes))
    partners = np.zeros((knockdowns, genes), dtype=bool)
    for k in range(knockdowns):
        chosen = rng.choice(genes, 20, replace=False) if specific else np.arange(20)
        partners[k, chosen] = True
        effects[k, chosen] -= 3
    return effects, partners


class PartnerTests(unittest.TestCase):
    def test_ranks_are_shares_of_the_row_with_ties_averaged(self):
        ranks = within_ranks(np.array([[3.0, 1.0, 1.0, 9.0]]))
        np.testing.assert_allclose(ranks, [[0.625, 0.25, 0.25, 0.875]])

    def test_partners_that_fall_only_under_their_own_knockdown_are_specific(self):
        effects, partners = screen(specific=True)
        means, sizes = partner_means(within_ranks(-effects), partners)
        own, share = specificity(means)
        result = summarise(own, share, np.ones(len(own), dtype=bool))
        self.assertTrue((sizes == 20).all())
        self.assertGreater(result["median_partner_rank"], 0.9)
        self.assertEqual(result["partners_respond"], 1.0)
        self.assertGreater(result["most_specific"], 0.95)

    def test_a_programme_that_falls_under_every_knockdown_is_not_specific(self):
        effects, partners = screen(specific=False)
        means, _ = partner_means(within_ranks(-effects), partners)
        own, share = specificity(means)
        result = summarise(own, share, np.ones(len(own), dtype=bool))
        self.assertGreater(result["median_partner_rank"], 0.9)  # the partners do fall
        self.assertLess(abs(result["median_specificity"] - 0.5), 0.1)  # but under any knockdown
        self.assertLess(result["most_specific"], 0.15)

    def test_knockdowns_with_few_partners_are_left_out(self):
        effects, partners = screen(specific=True)
        partners[0] = False
        partners[0, : MIN_PARTNERS - 1] = True
        means, _ = partner_means(within_ranks(-effects), partners)
        own, share = specificity(means)
        self.assertTrue(np.isnan(own[0]))
        result = summarise(own, share, np.ones(len(own), dtype=bool))
        self.assertEqual(result["knockdowns"], len(own) - 1)
        self.assertEqual(summarise(own, share, np.arange(len(own)) < 5)["knockdowns"], 4)


if __name__ == "__main__":
    unittest.main()
