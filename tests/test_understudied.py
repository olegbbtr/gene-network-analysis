import unittest

import numpy as np

from gene_network_analysis.understudied import by_group, terciles


class UnderstudiedTests(unittest.TestCase):
    def test_equal_values_share_a_group(self):
        group = terciles(np.array([0, 0, 0, 0, 1, 2, 5, 9, 9]))
        self.assertEqual(list(group), [0, 0, 0, 0, 1, 1, 2, 2, 2])

    def test_gain_is_measured_inside_each_group(self):
        rng = np.random.default_rng(0)
        group = np.repeat([0, 1, 2], 200)
        labels = rng.random(600) < 0.3
        noise = rng.normal(size=600)
        informed = np.where(group == 0, labels + 0.1 * noise, noise)
        result = by_group(labels, {"string": noise, "string+network": informed}, group, 100)
        self.assertGreater(result["fewest"]["gain"]["ci95"][0], 0.3)
        self.assertEqual(result["most"]["gain"]["gain"], 0.0)
        self.assertEqual(result["middle"]["genes"], 200)


if __name__ == "__main__":
    unittest.main()
