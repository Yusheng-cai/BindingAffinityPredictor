"""Focused tests for FlashBind ID mapping and native-output validation."""

import json
import tempfile
import unittest
from pathlib import Path

from affinity_benchmark.adapters.affinity_models import (
    flashbind_record_id,
    load_flashbind_predictions,
)


class FlashBindAdapterTests(unittest.TestCase):
    def test_record_id_matches_upstream_fep4_convention(self) -> None:
        sample = {"target_id": "tyk2", "ligand": {"name": "lig_ejm_31"}}
        self.assertEqual(flashbind_record_id(sample), "tyk2_lig-ejm-31")

    def test_successful_prediction_preserves_native_values(self) -> None:
        document = {
            "tyk2_lig-ejm-31": {
                "status": "success",
                "pred_value": -1.25,
                "pred_value_raw": -0.75,
                "mw": 412.3,
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "predictions.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            self.assertEqual(load_flashbind_predictions(path), document)

    def test_nonfinite_native_value_is_rejected(self) -> None:
        document = {
            "tyk2_lig-ejm-31": {
                "status": "success",
                "pred_value": float("nan"),
                "pred_value_raw": -0.75,
                "mw": 412.3,
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "predictions.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "invalid pred_value"):
                load_flashbind_predictions(path)


if __name__ == "__main__":
    unittest.main()
