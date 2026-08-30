import importlib.util
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "compare_rnp_paired_models.py"
SPEC = importlib.util.spec_from_file_location("compare_rnp_paired_models", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TestRnpPairedComparison(unittest.TestCase):
    def test_improvement_direction_respects_metric_orientation(self):
        nesso = np.asarray([3.0, 2.0, 1.0])
        boltz = np.asarray([2.0, 2.0, 2.0])
        lower = MODULE.paired_comparison(
            nesso,
            boltz,
            higher_is_better=False,
            rng=np.random.default_rng(1),
            iterations=20,
        )
        higher = MODULE.paired_comparison(
            nesso,
            boltz,
            higher_is_better=True,
            rng=np.random.default_rng(1),
            iterations=20,
        )
        self.assertEqual(lower["boltz2_better"], 1)
        self.assertEqual(lower["nesso1_better"], 1)
        self.assertEqual(higher["boltz2_better"], 1)
        self.assertEqual(higher["nesso1_better"], 1)
        self.assertEqual(lower["ties"], 1)

    def test_binned_familiarity_summary_uses_fixed_bins_and_includes_100(self):
        records = [
            {
                "familiarity_score_0_to_100": score,
                "metrics": {"metric": value},
            }
            for score, value in [(0, 1.0), (19.9, 3.0), (20, 5.0), (100, 9.0)]
        ]
        result = MODULE.binned_familiarity_summary(
            records,
            "metric",
            rng=np.random.default_rng(2),
            iterations=20,
        )
        self.assertEqual([entry["n"] for entry in result], [2, 1, 0, 0, 1])
        self.assertEqual(result[0]["median"], 2.0)
        self.assertEqual(result[-1]["median"], 9.0)


if __name__ == "__main__":
    unittest.main()
