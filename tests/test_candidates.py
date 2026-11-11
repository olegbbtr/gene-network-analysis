import unittest

import numpy as np

from gene_network_analysis.candidates import summarise


class CandidateTests(unittest.TestCase):
    def test_excludes_both_lists_and_keeps_rank_and_unscreened_denominator(self):
        rows = [
            {"gene": "known", "depmap": False},
            {"gene": "common", "depmap": True},
            {"gene": "candidate", "depmap": False},
            {"gene": "unscreened", "depmap": False},
        ]
        effects = np.array(
            [
                [-1, -1, -0.5, np.nan],
                [-1, -1, -0.8, np.nan],
                [-1, -1, 0, np.nan],
                [-1, -1, np.nan, np.nan],
            ]
        )
        first, second = summarise(rows, {"known"}, {"common"}, effects)
        self.assertEqual(first["rank_outside_cegv2"], 3)
        self.assertEqual(first["screened_lines"], 3)
        self.assertEqual(first["dependent_lines"], 2)
        self.assertAlmostEqual(first["dependent_share"], 2 / 3)
        self.assertEqual(second["screened_lines"], 0)
        self.assertIsNone(second["dependent_share"])

    def test_reference_drift_is_reported_instead_of_changing_candidate_identity(self):
        with self.assertRaisesRegex(ValueError, "membership differs"):
            summarise([{"gene": "a", "depmap": False}], set(), {"a"}, np.zeros((3, 1)))


if __name__ == "__main__":
    unittest.main()
