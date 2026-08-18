"""Focused tests for experimental-reference geometry helpers."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

from analyze_insr_experimental_reference import (  # noqa: E402
    component_contact_set,
    component_coverage,
    minimum_cross_distance,
)


class ExperimentalReferenceMetricTests(unittest.TestCase):
    def test_component_contact_set_preserves_target_numbering(self) -> None:
        target = {
            4: np.asarray([0.0, 0.0, 0.0]),
            19: np.asarray([12.0, 0.0, 0.0]),
            51: np.asarray([30.0, 0.0, 0.0]),
        }
        component = np.asarray([[2.0, 0.0, 0.0], [16.0, 0.0, 0.0]])
        self.assertEqual(component_contact_set(target, component, 5.0), {4, 19})

    def test_component_coverage_is_directional_and_inclusive(self) -> None:
        component = np.asarray([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
        binder = np.asarray([[2.0, 0.0, 0.0]])
        self.assertEqual(component_coverage(component, binder, 8.0), 1.0)
        self.assertEqual(component_coverage(component, binder, 2.0), 0.5)

    def test_minimum_cross_distance(self) -> None:
        first = np.asarray([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
        second = np.asarray([[4.0, 0.0, 0.0]])
        self.assertEqual(minimum_cross_distance(first, second), 4.0)


if __name__ == "__main__":
    unittest.main()
