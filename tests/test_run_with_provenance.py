import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_with_provenance.py"
SPEC = importlib.util.spec_from_file_location("run_with_provenance", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TestRunWithProvenance(unittest.TestCase):
    def test_file_record_is_absolute_and_content_addressed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            payload = b'{"schema_version": 1}\n'
            path.write_bytes(payload)
            record = MODULE.file_record(path)
            self.assertEqual(record["path"], str(path.resolve()))
            self.assertEqual(record["size_bytes"], len(payload))
            self.assertEqual(record["sha256"], hashlib.sha256(payload).hexdigest())


if __name__ == "__main__":
    unittest.main()
