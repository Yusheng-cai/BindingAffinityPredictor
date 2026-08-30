import importlib.util
import tempfile
import unittest
from pathlib import Path

import yaml


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_rnp_paired_inputs.py"
SPEC = importlib.util.spec_from_file_location("prepare_rnp_paired_inputs", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TestRnpPairedInputPreparation(unittest.TestCase):
    def setUp(self):
        self.sample = {
            "sample_id": "test__sample",
            "target_id": "test__target",
            "protein": {"chains": [{"id": "1.A", "sequence": "ACDE"}]},
            "ligand": {"input_smiles": "C[C@H](O)C(=O)O"},
        }

    def test_nesso_has_sequence_and_smiles_but_no_coordinates(self):
        document = MODULE.document_for_sample(self.sample, "nesso1", None)
        self.assertEqual(document["sequences"][0]["protein"]["sequence"], "ACDE")
        self.assertEqual(
            document["sequences"][1]["ligand"]["smiles"],
            "C[C@H](O)C(=O)O",
        )
        self.assertNotIn("msa", document["sequences"][0]["protein"])
        serialized = yaml.safe_dump(document)
        self.assertNotIn("pdb", serialized.lower())
        self.assertNotIn("template", serialized.lower())

    def test_boltz_requires_and_resolves_existing_msa(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            msa = root / "test__target" / "test__target_1.A.csv"
            msa.parent.mkdir(parents=True)
            msa.write_text("key,sequence\n0,ACDE\n", encoding="utf-8")
            document = MODULE.document_for_sample(self.sample, "boltz2", root)
            protein = document["sequences"][0]["protein"]
            self.assertEqual(document["version"], 1)
            self.assertEqual(protein["msa"], str(msa.resolve()))

    def test_boltz_rejects_missing_msa(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                MODULE.document_for_sample(self.sample, "boltz2", Path(tmp))


if __name__ == "__main__":
    unittest.main()
