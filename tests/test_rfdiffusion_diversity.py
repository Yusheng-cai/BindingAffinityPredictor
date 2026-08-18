"""Focused checks for the RFdiffusion structural-diversity metrics."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from analyze_rfdiffusion_diversity import (  # noqa: E402
    apply_transform,
    fit_transform,
    intrinsic_trace_rmsd,
    jaccard,
    resample_trace,
    rmsd,
    symmetric_chamfer,
)


class StructuralDiversityMetricTests(unittest.TestCase):
    def setUp(self) -> None:
        self.trace = np.asarray(
            [[0.0, 0.0, 0.0], [1.0, 0.3, 0.2], [2.2, 1.0, -0.4], [3.0, 1.8, 0.5], [4.3, 2.0, 1.1]]
        )
        angle = math.radians(37.0)
        self.rotation = np.asarray(
            [[math.cos(angle), -math.sin(angle), 0.0], [math.sin(angle), math.cos(angle), 0.0], [0.0, 0.0, 1.0]]
        )
        self.transformed = self.trace @ self.rotation + np.asarray([8.0, -3.0, 2.5])

    def test_kabsch_fit_recovers_a_rigid_transform(self) -> None:
        rotation, translation = fit_transform(self.transformed, self.trace)
        recovered = apply_transform(self.transformed, rotation, translation)
        self.assertLess(rmsd(recovered, self.trace), 1e-10)

    def test_intrinsic_trace_rmsd_removes_global_pose(self) -> None:
        self.assertLess(intrinsic_trace_rmsd(self.trace, self.transformed), 1e-10)

    def test_resampling_preserves_trace_endpoints(self) -> None:
        sampled = resample_trace(self.trace, count=64)
        np.testing.assert_allclose(sampled[0], self.trace[0])
        np.testing.assert_allclose(sampled[-1], self.trace[-1])
        self.assertEqual(sampled.shape, (64, 3))

    def test_chamfer_is_symmetric(self) -> None:
        shifted = self.trace + np.asarray([0.5, 0.0, 0.0])
        self.assertAlmostEqual(
            symmetric_chamfer(self.trace, shifted),
            symmetric_chamfer(shifted, self.trace),
            places=12,
        )

    def test_contact_jaccard(self) -> None:
        self.assertAlmostEqual(jaccard({1, 2, 3}, {2, 3, 4}), 0.5)
        self.assertEqual(jaccard(set(), set()), 1.0)


if __name__ == "__main__":
    unittest.main()
