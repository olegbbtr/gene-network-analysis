import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from gene_network_analysis.perturbation import (
    AFFECTED,
    MIN_CELLS,
    by_links,
    direction,
    one_per_gene,
    read_screen,
    tf_odds,
)


class PerturbationTests(unittest.TestCase):
    def test_the_knockdown_that_lowers_its_target_most_is_kept(self):
        names = [
            "1_A_P1_ENSG01",
            "2_A_P2_ENSG01",
            "3_B_P1_ENSG02",
            "4_B_P2_ENSG02",
            "5_non-targeting_non-targeting_non-targeting",
            "6_C_P1_nan",
            "7_D_P1_ENSG04",
        ]
        knockdown = np.array([-0.5, -0.9, np.nan, np.nan, 0.0, -0.8, -0.9])
        cells = np.array([200.0, 60.0, 30.0, 90.0, 500.0, 100.0, MIN_CELLS - 1])
        self.assertEqual(one_per_gene(names, knockdown, cells), {"ENSG01": 1, "ENSG02": 3})

    def test_a_screen_file_is_read_as_one_row_per_target(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "screen.h5ad"
            with h5py.File(path, "w") as store:
                store["X"] = np.arange(12, dtype=np.float32).reshape(4, 3)
                names = [
                    b"1_A_P1_ENSG01",
                    b"2_A_P2_ENSG01",
                    b"3_x_x_non-targeting",
                    b"4_B_P1_ENSG02",
                ]
                store["obs/gene_transcript"] = names
                store["obs/pct_expr"] = np.array([-0.2, -0.7, 0.0, np.nan])
                store["obs/num_cells_filtered"] = np.array([100.0, 100.0, 100.0, 50.0])
                store["obs/control_expr"] = np.array([1.0, 1.0, np.nan, 3.0])
                store["obs/anderson_darling_counts"] = np.array([5, 40, 0, 7])
                store["var/gene_id"] = [b"ENSG01", b"ENSG02", b"ENSG03"]
            screen = read_screen(path)
        self.assertEqual(screen.targets, ["ENSG01", "ENSG02"])
        self.assertEqual(screen.measured, ["ENSG01", "ENSG02", "ENSG03"])
        np.testing.assert_array_equal(screen.effects, [[3, 4, 5], [9, 10, 11]])
        np.testing.assert_array_equal(screen.changed, [40, 7])
        np.testing.assert_allclose(screen.covariates[1], [-0.2, np.log(50), 3.0])  # median fill

    def test_tf_odds_compare_tfs_with_other_knockdowns(self):
        broad = np.array([True, True, False, False, True, False, False, False, True])
        tf = np.array([True, True, True, False, False, False, False, False, True])
        keep = np.array([True] * 8 + [False])
        result = tf_odds(broad, tf, keep)
        self.assertEqual((result["tfs"], result["broad_tfs"], result["broad_others"]), (3, 2, 1))
        self.assertAlmostEqual(result["odds_ratio"], 8.0)

    def test_linked_pairs_with_strong_effects_stand_out(self):
        rng = np.random.default_rng(0)
        effects = rng.normal(size=(200, 300))
        count = np.zeros(effects.shape, dtype=np.int8)
        count[:20, :30] = 2  # 600 pairs linked in two networks
        effects[:20, :30] -= 6  # knocking these down lowers their linked genes
        valid = np.ones(effects.shape, dtype=bool)
        valid[:, 299] = False
        result = by_links(effects, count, valid, 3)
        linked, unlinked = result["by_networks"][2], result["by_networks"][0]
        self.assertEqual((linked["pairs"], result["by_networks"][1]["pairs"]), (600, 0))
        self.assertGreater(linked["affected"], 0.9)
        self.assertLess(unlinked["affected"], 0.001)
        self.assertEqual(linked["falls_among_affected"], 1.0)
        share = sum(r["pairs"] * (r["affected"] or 0) for r in result["by_networks"])
        self.assertAlmostEqual(share / valid.sum(), AFFECTED, places=3)

    def test_direction_counts_one_way_edges_with_a_strong_effect(self):
        edges = np.array([[0, 1, 1, 0], [0, 0, 0, 0], [1, 0, 0, 1], [0, 0, 0, 0]], dtype=bool)
        effects = np.array(
            [[0, 5.0, 9.0, 0], [1.0, 0, 0, 0], [9.0, 0, 0, 0.5], [0, 0, 4.0, 0]], dtype=float
        )
        valid = ~np.eye(4, dtype=bool)
        result = direction(edges, effects, valid, 3.0)
        # 0->1 agrees (5 against 1); 2->3 does not (0.5 against 4); 0<->2 is two-way
        self.assertEqual((result["pairs"], result["agreement"]), (2, 0.5))


if __name__ == "__main__":
    unittest.main()
