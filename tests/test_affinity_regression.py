"""Focused tests for affinity benchmark metric definitions."""

import unittest

import numpy as np

from affinity_benchmark.metrics.affinity_regression import (
    assay_metrics,
    paired_bootstrap_weighted_metric_difference,
    compound_weighted_average,
    log10_micromolar_to_kcal_per_mol,
    pairwise_mae,
)


class AffinityRegressionTests(unittest.TestCase):
    def test_log10_micromolar_conversion(self) -> None:
        values = log10_micromolar_to_kcal_per_mol([0.0, -3.0])
        self.assertAlmostEqual(values[0], -8.185, places=3)
        self.assertAlmostEqual(values[1], -12.278, places=3)

    def test_pairwise_mae_is_offset_invariant(self) -> None:
        observed = [1.0, 2.0, 4.0]
        predicted = [6.0, 7.0, 9.0]
        self.assertAlmostEqual(pairwise_mae(observed, predicted), 0.0)

    def test_centered_mae_removes_constant_offset(self) -> None:
        result = assay_metrics([1.0, 2.0, 4.0], [6.0, 7.0, 9.0])
        self.assertAlmostEqual(result["centered_mae_kcal_mol"], 0.0)
        self.assertAlmostEqual(result["pearson_r"], 1.0)

    def test_compound_weighted_average(self) -> None:
        metrics = {"a": {"n": 1, "pearson_r": 0.0}, "b": {"n": 3, "pearson_r": 1.0}}
        self.assertAlmostEqual(compound_weighted_average(metrics, "pearson_r"), 0.75)

    def test_paired_bootstrap_preserves_model_pairing(self) -> None:
        observed = np.arange(10, dtype=float)
        result = paired_bootstrap_weighted_metric_difference(
            {"target": (observed, observed, -observed)},
            "pearson_r",
            iterations=100,
            seed=7,
        )
        self.assertAlmostEqual(result["lower_95"], 2.0)
        self.assertAlmostEqual(result["median"], 2.0)
        self.assertAlmostEqual(result["upper_95"], 2.0)


if __name__ == "__main__":
    unittest.main()
