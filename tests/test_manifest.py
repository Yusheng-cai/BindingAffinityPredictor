"""Tests for the frozen exp001 manifest and validation failures."""

import copy
import json
import unittest
from pathlib import Path

from affinity_benchmark.data.manifest import load_manifest, validate_manifest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPOSITORY_ROOT / "data" / "manifests" / "exp001_tyk2_4gih.json"
UNLABELED_MANIFEST_PATH = (
    REPOSITORY_ROOT / "data" / "manifests" / "exp002_gal3_9d63_galactose.json"
)


class ManifestTests(unittest.TestCase):
    def test_exp001_manifest_is_valid(self) -> None:
        manifest = load_manifest(MANIFEST_PATH)
        sample = manifest["samples"][0]
        self.assertEqual(sample["protein"]["sequence_length"], 302)
        self.assertEqual(sample["measurement"]["type"], "ki")
        self.assertFalse(sample["structure_reference"]["use_as_model_input"])

    def test_unlabeled_structural_manifest_is_valid(self) -> None:
        manifest = load_manifest(UNLABELED_MANIFEST_PATH)
        sample = manifest["samples"][0]
        self.assertIsNone(sample["measurement"])
        self.assertEqual(sample["ligand"]["ccd_id"], "GAL")
        self.assertFalse(sample["structure_reference"]["use_as_model_input"])

    def test_sequence_length_mismatch_is_rejected(self) -> None:
        manifest = self._raw_manifest()
        manifest["samples"][0]["protein"]["sequence_length"] = 301
        with self.assertRaisesRegex(ValueError, "sequence has 302 residues"):
            validate_manifest(manifest)

    def test_inconsistent_log_value_is_rejected(self) -> None:
        manifest = self._raw_manifest()
        manifest["samples"][0]["measurement"]["pX_M"] = 2.0
        with self.assertRaisesRegex(ValueError, "expected 8.318"):
            validate_manifest(manifest)

    def test_reference_cannot_be_enabled_as_input(self) -> None:
        manifest = self._raw_manifest()
        manifest["samples"][0]["structure_reference"]["use_as_model_input"] = True
        with self.assertRaisesRegex(ValueError, "must be false"):
            validate_manifest(manifest)

    def test_multichain_protein_is_valid(self) -> None:
        manifest = self._raw_manifest()
        manifest["samples"][0]["protein"] = {
            "name": "two-chain test construct",
            "chains": [
                {"id": "A", "sequence": "ACDE", "sequence_length": 4},
                {"id": "B", "sequence": "FGHI", "sequence_length": 4},
            ],
        }
        validate_manifest(manifest)

    def test_duplicate_protein_chain_id_is_rejected(self) -> None:
        manifest = self._raw_manifest()
        manifest["samples"][0]["protein"] = {
            "chains": [
                {"id": "A", "sequence": "ACDE", "sequence_length": 4},
                {"id": "A", "sequence": "FGHI", "sequence_length": 4},
            ]
        }
        with self.assertRaisesRegex(ValueError, "duplicates 'A'"):
            validate_manifest(manifest)

    @staticmethod
    def _raw_manifest() -> dict:
        with MANIFEST_PATH.open(encoding="utf-8") as handle:
            return copy.deepcopy(json.load(handle))


if __name__ == "__main__":
    unittest.main()
