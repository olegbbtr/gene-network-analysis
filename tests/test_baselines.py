import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

from gene_network_analysis.baselines import mean_correlation, read_string


class BaselineTests(unittest.TestCase):
    def test_mean_correlation_excludes_the_gene_itself(self):
        rng = np.random.default_rng(0)
        base = rng.normal(size=(200, 1))
        samples = np.hstack([base, base, rng.normal(size=(200, 1))])
        samples = (samples - samples.mean(0)) / samples.std(0)
        known = np.array([True, False, False])
        scores = mean_correlation(samples, known)
        self.assertAlmostEqual(scores[1], 1.0)
        self.assertEqual(scores[0], 0.0)  # its only known gene is itself
        self.assertLess(abs(scores[2]), 0.3)

    def test_string_keeps_confident_links_between_listed_genes(self):
        with tempfile.TemporaryDirectory() as directory:
            info, links = Path(directory) / "info.gz", Path(directory) / "links.gz"
            with gzip.open(info, "wt") as stream:
                stream.write("#id\tpreferred_name\tsize\tannotation\n")
                stream.write("P.A\tA\t1\tx\nP.B\tB\t1\tx\nP.C\tC\t1\tx\n")
            with gzip.open(links, "wt") as stream:
                stream.write("protein1 protein2 combined_score\n")
                stream.write("P.A P.B 900\nP.B P.A 900\nP.A P.C 400\n")
            matrix = read_string(
                info, links, ["g1", "g2", "g3"], {"A": "g1", "B": "g2", "C": "g3"}.get
            )
        self.assertEqual(matrix.nnz, 2)
        self.assertEqual(matrix[0, 1], 1.0)


if __name__ == "__main__":
    unittest.main()
