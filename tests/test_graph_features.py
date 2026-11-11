import unittest

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.graph_features import (
    core_depth,
    eigenvector_score,
    rank_positions,
    weight_summaries,
)


class GraphFeatureTests(unittest.TestCase):
    def test_rank_positions_average_ties_and_keep_highest_score_first(self):
        scores = np.array([[9.0, 1.0], [9.0, 2.0], [3.0, 3.0], [0.0, 4.0]])
        expected = np.array([[1 / 6, 1], [1 / 6, 2 / 3], [2 / 3, 1 / 3], [1, 0]])
        np.testing.assert_allclose(rank_positions(scores), expected)
        np.testing.assert_array_equal(rank_positions(scores), rank_positions(np.exp(scores)))

    def test_clique_survives_at_depth_three_but_its_leaf_does_not(self):
        matrix = np.zeros((6, 6))
        matrix[:4, :4] = 1
        np.fill_diagonal(matrix, 0)
        matrix[0, 4] = matrix[4, 0] = 1
        np.testing.assert_array_equal(core_depth(csr_array(matrix)), [3, 3, 3, 3, 1, 0])

    def test_star_has_no_dense_core_even_at_its_central_node(self):
        matrix = np.zeros((8, 8))
        matrix[0, 1:] = matrix[1:, 0] = 1
        np.testing.assert_array_equal(core_depth(csr_array(matrix)), np.ones(8))
        permutation = np.array([7, 3, 1, 0, 6, 2, 5, 4])
        np.testing.assert_array_equal(core_depth(csr_array(matrix[permutation][:, permutation])), 1)

    def test_eigenvector_score_favours_the_star_centre_and_is_normalised(self):
        matrix = np.zeros((5, 5))
        matrix[0, 1:] = matrix[1:, 0] = 1
        scores = eigenvector_score(csr_array(matrix))
        self.assertEqual(scores.argmax(), 0)
        self.assertAlmostEqual(scores.sum(), 1)
        np.testing.assert_allclose(scores[1:], scores[1])
        np.testing.assert_allclose(eigenvector_score(csr_array((3, 3))), np.full(3, 1 / 3))

    def test_weight_distributions_do_not_depend_on_neighbour_identity(self):
        first = np.zeros((5, 5))
        second = first.copy()
        first[0, 1:3] = [2, -4]
        second[0, 3:5] = [-4, 2]
        a, b = (weight_summaries(csr_array(matrix)) for matrix in (first, second))
        np.testing.assert_allclose(a[0, :4], [-1, 3, 3, 0.5])
        np.testing.assert_array_equal(a[0, :4], b[0, :4])
        np.testing.assert_array_equal(a[3], np.zeros(8))
        np.testing.assert_allclose(a[1, 4:], [2, 0, 2, 1])


if __name__ == "__main__":
    unittest.main()
