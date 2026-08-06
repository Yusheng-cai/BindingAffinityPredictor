"""Tests for explicit affinity-unit and logarithmic-scale handling."""

import math
import unittest

from affinity_benchmark.data.affinity import (
    concentration_to_log10_micromolar,
    concentration_to_molar,
    concentration_to_px,
    log10_micromolar_to_px,
    px_to_log10_micromolar,
)


class AffinityUnitTests(unittest.TestCase):
    def test_reference_ki_conversions(self) -> None:
        self.assertAlmostEqual(concentration_to_molar(4.8, "nM"), 4.8e-9)
        log_uM = concentration_to_log10_micromolar(0.0048, "uM")
        px = concentration_to_px(0.0048, "uM")
        self.assertAlmostEqual(log_uM, -2.318758762624413)
        self.assertAlmostEqual(px, 8.318758762624412)
        self.assertAlmostEqual(log10_micromolar_to_px(log_uM), px)
        self.assertAlmostEqual(px_to_log10_micromolar(px), log_uM)

    def test_equivalent_units_have_identical_logs(self) -> None:
        expected = concentration_to_log10_micromolar(4.8, "nM")
        self.assertAlmostEqual(expected, concentration_to_log10_micromolar(0.0048, "uM"))

    def test_rejects_invalid_values_and_units(self) -> None:
        for value in (0, -1, math.inf, math.nan):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    concentration_to_molar(value, "uM")
        with self.assertRaises(TypeError):
            concentration_to_molar(True, "uM")
        with self.assertRaises(ValueError):
            concentration_to_molar(1, "mol/L")


if __name__ == "__main__":
    unittest.main()
