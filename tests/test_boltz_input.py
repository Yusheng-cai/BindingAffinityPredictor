"""Guard the committed Boltz input against drift from the frozen manifest."""

import re
import unittest
from pathlib import Path

from affinity_benchmark.data.manifest import load_manifest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPOSITORY_ROOT / "data" / "manifests" / "exp001_tyk2_4gih.json"
INPUT_PATH = (
    REPOSITORY_ROOT
    / "experiments"
    / "exp001_boltz2_4gih_smoke"
    / "inputs"
    / "4gih_0x5.yaml"
)


class BoltzInputTests(unittest.TestCase):
    def test_sequence_and_smiles_match_manifest(self) -> None:
        sample = load_manifest(MANIFEST_PATH)["samples"][0]
        input_text = INPUT_PATH.read_text(encoding="utf-8")
        self.assertIn(f"sequence: {sample['protein']['sequence']}", input_text)
        self.assertIn(f"smiles: '{sample['ligand']['input_smiles']}'", input_text)

    def test_holo_information_is_absent(self) -> None:
        input_text = INPUT_PATH.read_text(encoding="utf-8")
        for forbidden_key in ("templates", "constraints", "pocket", "contact", "msa"):
            with self.subTest(forbidden_key=forbidden_key):
                self.assertIsNone(
                    re.search(rf"^\s*{forbidden_key}:\s*", input_text, flags=re.MULTILINE)
                )
        self.assertNotIn("4GIH", input_text.upper())


if __name__ == "__main__":
    unittest.main()
