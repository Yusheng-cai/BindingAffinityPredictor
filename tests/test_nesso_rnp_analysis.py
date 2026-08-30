import unittest

import numpy as np


try:
    from rdkit import Chem

    from scripts.analyze_nesso_rnp_distograms import (
        _cluster_bootstrap_interval,
        familiarity_metadata,
        remove_nesso_tokenization_hydrogens,
    )

    NESSO_ANALYSIS_DEPS_AVAILABLE = True
except ModuleNotFoundError:
    NESSO_ANALYSIS_DEPS_AVAILABLE = False


@unittest.skipUnless(
    NESSO_ANALYSIS_DEPS_AVAILABLE, "requires the separate Nesso analysis environment"
)
class NessoRunsNPosesAnalysisTests(unittest.TestCase):
    def test_nested_nesso_familiarity_uses_september_2021_score(self):
        sample = {
            "runs_n_poses": {
                "similarity_score_0_to_100": 88.0,
                "similarity_bin": "80-100",
                "nesso1_sep2021_familiarity": {
                    "similarity_score_0_to_100": 37.5
                },
            }
        }
        score, bin_name = familiarity_metadata(sample, "nesso1_sep2021")
        self.assertEqual(score, 37.5)
        self.assertEqual(bin_name, "30-40")

    def test_cluster_bootstrap_keeps_cluster_records_together(self):
        values = np.asarray([1.0, 2.0, 10.0])
        clusters = np.asarray(["shared", "shared", "single"], dtype=object)

        def assert_complete_blocks(sample):
            self.assertEqual(np.sum(sample == 1.0), np.sum(sample == 2.0))
            return float(np.mean(sample))

        result = _cluster_bootstrap_interval(
            values,
            clusters,
            assert_complete_blocks,
            iterations=100,
            rng=np.random.default_rng(12),
        )
        self.assertEqual(result["iterations_finite"], 100)

    def test_explicit_stereo_hydrogens_are_not_counted_as_nesso_tokens(self):
        smiles = "[H]/N=C(\\N)Nc1ccc(C(=O)Oc2ccc3cc(/C(N)=N\\[H])ccc3c2)cc1"
        molecule = Chem.MolFromSmiles(smiles)
        original_heavy_symbols = [
            atom.GetSymbol() for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1
        ]
        token_molecule = remove_nesso_tokenization_hydrogens(molecule)
        token_symbols = [atom.GetSymbol() for atom in token_molecule.GetAtoms()]

        self.assertEqual(molecule.GetNumAtoms(), 28)
        self.assertEqual(token_molecule.GetNumAtoms(), 26)
        self.assertEqual(token_symbols, original_heavy_symbols)


if __name__ == "__main__":
    unittest.main()
