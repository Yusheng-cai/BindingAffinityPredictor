"""Focused checks for hotspot-ablation geometry helpers."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from analyze_rfdiffusion_hotspot_ablation import (  # noqa: E402
    contact_set_and_distances,
    rounded_summary,
)


class HotspotAblationMetricTests(unittest.TestCase):
    def test_contact_set_uses_target_residue_numbers(self) -> None:
        target = {
            10: np.asarray([0.0, 0.0, 0.0]),
            25: np.asarray([20.0, 0.0, 0.0]),
            90: np.asarray([40.0, 0.0, 0.0]),
        }
        binder = np.asarray([[2.0, 0.0, 0.0], [37.0, 0.0, 0.0]])
        contacts, distances = contact_set_and_distances(target, binder, threshold=5.0)
        self.assertEqual(contacts, {10, 90})
        self.assertEqual(distances.shape, (2, 3))

    def test_contact_threshold_is_inclusive(self) -> None:
        target = {7: np.asarray([0.0, 0.0, 0.0])}
        binder = np.asarray([[0.0, 0.0, 10.0]])
        contacts, _ = contact_set_and_distances(target, binder, threshold=10.0)
        self.assertEqual(contacts, {7})

    def test_rounded_summary_reports_prespecified_statistics(self) -> None:
        self.assertEqual(
            rounded_summary([1.0, 2.0, 9.0]),
            {"minimum": 1.0, "median": 2.0, "mean": 4.0, "maximum": 9.0},
        )


if __name__ == "__main__":
    unittest.main()
