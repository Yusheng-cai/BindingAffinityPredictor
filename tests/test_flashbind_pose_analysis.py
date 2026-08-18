import importlib.util
import sys
import unittest
from pathlib import Path


HAS_MODEL_SPECIFIC_DEPENDENCIES = all(
    importlib.util.find_spec(module) is not None
    for module in ("gemmi", "lmdb", "rdkit")
)


@unittest.skipUnless(
    HAS_MODEL_SPECIFIC_DEPENDENCIES,
    "requires the isolated FlashBind pose-analysis dependencies",
)
class FlashBindPoseAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[1]
        scripts = root / "scripts"
        sys.path.insert(0, str(root / "src"))
        sys.path.insert(0, str(scripts))
        import analyze_flashbind_crystal_pose as analysis

        cls.analysis = analysis

    def test_blank_chain_id_fragment_is_joined_without_coordinate_changes(self):
        import gemmi

        model = gemmi.Model("1")
        for chain_id, start, stop in (("A", 1, 25), ("", 26, 50)):
            chain = gemmi.Chain(chain_id)
            for residue_number in range(start, stop + 1):
                residue = gemmi.Residue()
                residue.name = "ALA"
                residue.seqid = gemmi.SeqId(residue_number, " ")
                atom = gemmi.Atom()
                atom.name = "CA"
                atom.element = gemmi.Element("C")
                atom.pos = gemmi.Position(float(residue_number), 0.0, 0.0)
                residue.add_atom(atom)
                chain.add_residue(residue)
            model.add_chain(chain)

        chain, selection = self.analysis.choose_predicted_chain(model, "A" * 50)
        sequence, residues = self.analysis.observed_chain(chain)

        self.assertEqual(sequence, "A" * 50)
        self.assertEqual([residue.seqid.num for residue in residues], list(range(1, 51)))
        self.assertIsNotNone(selection["chain_id_repair"])
        self.assertEqual(selection["expected_sequence_coverage"], 1.0)


if __name__ == "__main__":
    unittest.main()
