import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


HAS_DEPS = all(
    importlib.util.find_spec(module) is not None
    for module in ("gemmi", "rdkit", "Bio", "safetensors")
)


@unittest.skipUnless(HAS_DEPS, "requires the Nesso/Boltz analysis dependencies")
class BoltzRunsNPosesAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[1]
        sys.path.insert(0, str(root / "src"))
        sys.path.insert(0, str(root / "scripts"))
        import analyze_boltz2_rnp_structures as analysis

        cls.analysis = analysis

    def test_chain_coordinates_preserve_sequence_indices(self):
        import gemmi

        chain = gemmi.Chain("A")
        for residue_number, residue_name in ((1, "ALA"), (3, "GLY")):
            residue = gemmi.Residue()
            residue.name = residue_name
            residue.seqid = gemmi.SeqId(residue_number, " ")
            for atom_name, xyz in (
                ("CA", (float(residue_number), 0.0, 0.0)),
                ("CB", (float(residue_number), 1.0, 0.0)),
            ):
                if residue_name == "GLY" and atom_name == "CB":
                    continue
                atom = gemmi.Atom()
                atom.name = atom_name
                atom.element = gemmi.Element("C")
                atom.pos = gemmi.Position(*xyz)
                residue.add_atom(atom)
            chain.add_residue(residue)

        token, ca, heavy = self.analysis.chain_coordinates(chain, "ACG")
        self.assertTrue(np.allclose(token[0], [1.0, 1.0, 0.0]))
        self.assertTrue(np.all(np.isnan(token[1])))
        self.assertTrue(np.allclose(token[2], [3.0, 0.0, 0.0]))
        self.assertTrue(np.allclose(ca[2], [3.0, 0.0, 0.0]))
        self.assertIsNone(heavy[1])

    def test_explicit_stereo_hydrogen_is_removed_before_coordinate_mapping(self):
        from rdkit import Chem

        molecule = Chem.MolFromSmiles("[H]/N=C(/N)c1ccccc1")
        cleaned = self.analysis.remove_nesso_tokenization_hydrogens(molecule)
        self.assertTrue(all(atom.GetAtomicNum() != 1 for atom in cleaned.GetAtoms()))


if __name__ == "__main__":
    unittest.main()
