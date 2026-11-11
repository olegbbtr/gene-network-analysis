import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

from gene_network_analysis.constraint import labels, read_constraint


class ConstraintTests(unittest.TestCase):
    def test_rows_without_numbers_are_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "constraint.txt.bgz"
            with gzip.open(path, "wt") as stream:
                stream.write("gene\toe_lof_upper\texp_lof\tgene_id\n")
                stream.write("A\t0.2\t40.5\tENSG1\nB\tNA\t3\tENSG2\nC\t1.1\t12\tENSG3\n")
            table = read_constraint(path)
        self.assertEqual(table, {"ENSG1": (0.2, 40.5), "ENSG3": (1.1, 12.0)})

    def test_short_genes_are_not_usable_and_missing_genes_are_not_constrained(self):
        table = {"a": (0.2, 40.0), "b": (0.2, 5.0), "c": (0.9, 20.0)}
        constrained, usable, expected = labels(table, ["a", "b", "c", "d"])
        self.assertEqual(list(constrained), [True, True, False, False])
        self.assertEqual(list(usable), [True, False, True, False])
        np.testing.assert_allclose(expected, np.log([40.0, 1.0, 20.0, 1.0]))


if __name__ == "__main__":
    unittest.main()
