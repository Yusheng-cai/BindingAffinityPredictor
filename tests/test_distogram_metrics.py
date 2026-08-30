import unittest

import numpy as np

from affinity_benchmark.metrics.distogram import (
    average_precision,
    binary_auroc,
    distance_bin_indices,
    distogram_centers,
    distogram_negative_log_likelihood,
    expected_distances,
    precision_recall_f1,
)


class DistogramMetricTests(unittest.TestCase):
    def test_nesso_centers_include_open_end_bins(self):
        centers = distogram_centers()
        self.assertEqual(len(centers), 64)
        self.assertEqual(centers[0], 1.5)
        self.assertEqual(centers[-1], 24.5)

    def test_distance_bins_include_exact_first_boundary_in_first_bin(self):
        bins = distance_bin_indices(np.array([1.0, 2.0, 2.01, 22.0, 30.0]))
        self.assertEqual(bins[0], 0)
        self.assertEqual(bins[1], 0)
        self.assertEqual(bins[2], 1)
        self.assertEqual(bins[-2], 62)
        self.assertEqual(bins[-1], 63)

    def test_nll_and_expected_distance_for_certain_prediction(self):
        logits = np.full((2, 64), -20.0)
        logits[0, 0] = 20.0
        logits[1, 63] = 20.0
        nll = distogram_negative_log_likelihood(logits, np.array([1.0, 30.0]))
        np.testing.assert_allclose(nll, 0.0, atol=1e-12)
        np.testing.assert_allclose(expected_distances(logits), [1.5, 24.5], atol=1e-12)

    def test_binary_ranking_metrics(self):
        labels = np.array([1, 0, 1, 0])
        perfect = np.array([0.9, 0.2, 0.8, 0.1])
        self.assertEqual(average_precision(labels, perfect), 1.0)
        self.assertEqual(binary_auroc(labels, perfect), 1.0)
        tied = np.ones(4)
        self.assertEqual(binary_auroc(labels, tied), 0.5)

    def test_precision_recall_f1(self):
        precision, recall, f1 = precision_recall_f1(
            np.array([1, 1, 0, 0]), np.array([1, 0, 1, 0])
        )
        self.assertEqual((precision, recall, f1), (0.5, 0.5, 0.5))


if __name__ == "__main__":
    unittest.main()
