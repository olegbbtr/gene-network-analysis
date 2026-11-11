import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.robustness import (
    TOP,
    discovery,
    held_out,
    near_misses,
    read_effects,
    transfer,
    transfer_summary,
)


def synthetic(size=400, seed=0):
    rng = np.random.default_rng(seed)
    labels = rng.random(size) < 0.3
    groups = np.array([f"chr{i % 8}" for i in range(size)])
    signal = labels + rng.normal(scale=0.8, size=size)
    parts = {
        "expression": rng.normal(size=(size, 2)),
        "degree": np.column_stack([signal, rng.normal(size=size)]),
        "shape": rng.normal(size=(size, 3)),
    }
    dense = np.triu(rng.random((size, size)) < 0.03, 1)
    links = csr_array((dense | dense.T).astype(float))
    return parts, links, labels, groups


class RobustnessTests(unittest.TestCase):
    def test_discovery_scores_only_genes_outside_the_training_list(self):
        known = np.array([True, True, False, False, False, False])
        target = np.array([True, True, True, False, True, False])
        scores = np.array([9.0, 8.0, 3.0, 2.0, 1.0, 0.0])
        result = discovery(known, target, {"expression": scores, "expression+network": scores}, 50)
        self.assertEqual(result["listed"], 2)
        self.assertEqual(result["base_rate"], 0.5)
        self.assertEqual(result["auc"]["expression"], 0.75)
        self.assertEqual(result["top_hits"]["expression"], 2)

    def test_removed_genes_do_not_influence_the_kept_ones(self):
        parts, links, labels, groups = synthetic()
        keep = np.arange(len(labels)) % 5 != 0
        first = held_out(parts, [links], labels, groups, keep)
        changed = {name: np.where(keep[:, None], block, 100.0) for name, block in parts.items()}
        second = held_out(changed, [links], labels | ~keep, groups, keep)
        for name, scores in first.items():
            self.assertEqual(len(scores), keep.sum())
            np.testing.assert_allclose(scores, second[name])

    def test_effects_are_placed_by_gene_and_unscreened_genes_stay_empty(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "effects.csv"
            path.write_text(",B (2),X (9),A (1)\nACH-1,-1.5,0.1,\nACH-2,-0.2,0.3,0.4\n")
            values = read_effects(path, {"A": "ga", "B": "gb", "C": "gc"}, ["ga", "gb", "gc"])
        np.testing.assert_array_equal(values, [[np.nan, -1.5, np.nan], [0.4, -0.2, np.nan]])

    def test_near_misses_are_the_top_predictions_outside_both_lists(self):
        size = TOP + 50
        known = np.zeros(size, dtype=bool)
        known[0] = True  # scored highest, but it is in the training list
        target = np.zeros(size, dtype=bool)
        target[1:TOP] = True  # so the top new predictions are genes 1..TOP; one is unlisted
        scores = -np.arange(size, dtype=float)
        effects = np.zeros((4, size))
        effects[:2, TOP] = -1.0  # the missed gene is needed by half of the cell lines
        effects[:, TOP + 1] = np.nan
        result = near_misses(scores, known, target, effects)
        self.assertEqual(
            result["missed"],
            {
                "genes": 1,
                "screened": 1,
                "median_share_of_lines": 0.5,
                "needed_by_a_quarter_of_lines": 1.0,
            },
        )
        self.assertEqual((result["others"]["genes"], result["others"]["screened"]), (49, 48))
        self.assertEqual(result["others"]["needed_by_a_quarter_of_lines"], 0.0)

    def test_transfer_losses_are_grouped_by_the_kinds_of_network(self):
        matrix = {
            "t1": {"t1": 0.9, "t2": 0.88, "h": 0.7},
            "t2": {"t1": 0.86, "t2": 0.9, "h": 0.6},
            "h": {"t1": 0.8, "t2": 0.75, "h": 0.85},
        }
        state = {"t1": "tumour", "t2": "tumour", "h": "healthy"}
        summary = transfer_summary(matrix, {"t1": 0.87, "t2": 0.7, "h": 0.65}, state)
        by_state = summary["median_loss_by_state"]
        self.assertEqual(
            set(by_state), {"tumour to tumour", "tumour to healthy", "healthy to tumour"}
        )
        self.assertAlmostEqual(by_state["tumour to tumour"], 0.03)  # losses 0.02 and 0.04
        self.assertAlmostEqual(by_state["tumour to healthy"], 0.2)  # 0.15 and 0.25
        self.assertAlmostEqual(by_state["healthy to tumour"], 0.125)  # 0.1 and 0.15
        self.assertEqual(summary["pairs_above_own_expression"], 3)  # t1>t2, t1>h, h>t2

    def test_transfer_ignores_each_dataset_own_scale(self):
        parts, links, labels, groups = synthetic()
        rescaled = {name: block * 10 + 3 for name, block in parts.items()}
        matrix = transfer({"a": parts, "b": rescaled}, {"a": links, "b": links}, labels, groups)
        self.assertAlmostEqual(matrix["a"]["b"], matrix["b"]["b"], places=6)
        self.assertAlmostEqual(matrix["a"]["a"], matrix["b"]["b"], places=6)
        self.assertGreater(matrix["a"]["b"], 0.65)


if __name__ == "__main__":
    unittest.main()
