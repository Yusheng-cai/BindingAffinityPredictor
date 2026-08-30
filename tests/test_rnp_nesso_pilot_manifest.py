"""Focused checks for the frozen Runs N' Poses Nesso-1 pilot."""

import json
import unittest
from collections import Counter
from pathlib import Path

from affinity_benchmark.data.manifest import validate_manifest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPOSITORY_ROOT / "data/manifests/rnp_nesso1_pilot100.json"


class RunsNPosesPilotManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with MANIFEST_PATH.open(encoding="utf-8") as handle:
            cls.manifest = json.load(handle)

    def test_generic_manifest_schema_is_valid(self) -> None:
        validate_manifest(self.manifest)

    def test_selection_is_blinded_and_has_one_hundred_samples(self) -> None:
        self.assertEqual(len(self.manifest["samples"]), 100)
        self.assertEqual(
            self.manifest["selection"]["status"], "frozen_before_nesso_inference"
        )
        self.assertTrue(
            all(
                sample["structure_reference"]["use_as_model_input"] is False
                for sample in self.manifest["samples"]
            )
        )

    def test_bin_counts_match_preregistered_quotas(self) -> None:
        observed = Counter(
            sample["runs_n_poses"]["similarity_bin"]
            for sample in self.manifest["samples"]
        )
        self.assertEqual(dict(observed), self.manifest["selection"]["bin_quotas"])

    def test_hard_diversity_identifiers_are_unique(self) -> None:
        getters = (
            lambda sample: sample["runs_n_poses"]["cluster"],
            lambda sample: sample["structure_reference"]["pdb_id"],
            lambda sample: sample["ligand"]["input_smiles"],
        )
        for getter in getters:
            values = [getter(sample) for sample in self.manifest["samples"]]
            self.assertEqual(len(values), len(set(values)))


if __name__ == "__main__":
    unittest.main()
