import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_boltz_msas.py"
SPEC = importlib.util.spec_from_file_location("validate_boltz_msas", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TestValidateBoltzMsas(unittest.TestCase):
    def test_query_sequence_and_depth(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "msa.csv"
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["key", "sequence"])
                writer.writerow([-1, "ACDE"])
                writer.writerow([-1, "AC-E"])
            record = MODULE.inspect_msa(path, "ACDE")
            self.assertEqual(record["rows"], 2)
            self.assertTrue(record["query_matches_manifest"])

    def test_wrong_query_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "msa.csv"
            path.write_text("key,sequence\n-1,AAAA\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                MODULE.inspect_msa(path, "ACDE")


if __name__ == "__main__":
    unittest.main()
