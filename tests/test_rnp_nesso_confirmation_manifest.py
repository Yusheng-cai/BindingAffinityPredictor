"""Checks for the frozen 50-system Nesso Runs N' Poses confirmation cohort."""

import json
import unittest
from collections import Counter
from pathlib import Path

from affinity_benchmark.data.manifest import validate_manifest


ROOT = Path(__file__).resolve().parents[1]
DISCOVERY = ROOT / "data/manifests/rnp_nesso1_pilot100.json"
CONFIRMATION = ROOT / "data/manifests/rnp_nesso1_confirmation50.json"


class RunsNPosesConfirmationManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.discovery = json.loads(DISCOVERY.read_text())
        cls.confirmation = json.loads(CONFIRMATION.read_text())

    def test_generic_schema_and_count(self) -> None:
        validate_manifest(self.confirmation)
        self.assertEqual(len(self.confirmation["samples"]), 50)

    def test_no_exact_sample_overlap(self) -> None:
        discovery_ids = {sample["sample_id"] for sample in self.discovery["samples"]}
        confirmation_ids = {
            sample["sample_id"] for sample in self.confirmation["samples"]
        }
        self.assertFalse(discovery_ids & confirmation_ids)

    def test_expected_bin_counts(self) -> None:
        observed = Counter(
            sample["runs_n_poses"]["similarity_bin"]
            for sample in self.confirmation["samples"]
        )
        self.assertEqual(dict(observed), self.confirmation["selection"]["bin_quotas"])

    def test_within_cohort_diversity(self) -> None:
        getters = (
            lambda sample: sample["runs_n_poses"]["cluster"],
            lambda sample: sample["structure_reference"]["pdb_id"],
            lambda sample: sample["ligand"]["input_smiles"],
        )
        for getter in getters:
            values = [getter(sample) for sample in self.confirmation["samples"]]
            self.assertEqual(len(values), len(set(values)))


if __name__ == "__main__":
    unittest.main()
