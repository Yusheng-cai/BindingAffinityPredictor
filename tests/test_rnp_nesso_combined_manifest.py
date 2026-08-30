"""Checks for the pooled 150-system Nesso Runs N' Poses analysis."""

import json
import unittest
from collections import Counter
from pathlib import Path

from affinity_benchmark.data.manifest import validate_manifest


ROOT = Path(__file__).resolve().parents[1]
COMBINED = ROOT / "data/manifests/rnp_nesso1_combined150.json"
COMPONENTS = (
    ROOT / "data/manifests/rnp_nesso1_pilot100.json",
    ROOT / "data/manifests/rnp_nesso1_confirmation50.json",
)
RESULTS = ROOT / "reports/exp014_nesso1_rnp_combined150/results"


class RunsNPosesCombinedManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.combined = json.loads(COMBINED.read_text())
        cls.components = [json.loads(path.read_text()) for path in COMPONENTS]
        cls.summary = json.loads((RESULTS / "summary.json").read_text())

    def test_manifest_is_valid_union(self) -> None:
        validate_manifest(self.combined)
        combined_ids = {sample["sample_id"] for sample in self.combined["samples"]}
        component_ids = {
            sample["sample_id"]
            for manifest in self.components
            for sample in manifest["samples"]
        }
        self.assertEqual(combined_ids, component_ids)
        self.assertEqual(len(combined_ids), 150)

    def test_combined_bin_counts(self) -> None:
        observed = Counter(
            sample["runs_n_poses"]["similarity_bin"]
            for sample in self.combined["samples"]
        )
        self.assertEqual(dict(sorted(observed.items())), self.combined["selection"]["bin_counts"])

    def test_all_records_scored(self) -> None:
        aggregate = self.summary["aggregate"]
        self.assertEqual(aggregate["systems_requested"], 150)
        self.assertEqual(aggregate["systems_scored"], 150)
        self.assertEqual(aggregate["systems_failed"], 0)
        self.assertIn("all records in a sampled cluster retained", aggregate["bootstrap"]["unit"])


if __name__ == "__main__":
    unittest.main()
