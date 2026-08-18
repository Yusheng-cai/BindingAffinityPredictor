import unittest

import numpy as np

from scripts.compare_pose_models import describe, summarize_rmsd_components


class PoseComponentSummaryTests(unittest.TestCase):
    def test_describe_reports_expected_statistics(self) -> None:
        result = describe(np.array([1.0, 2.0, 5.0]))
        self.assertEqual(result["mean"], 8.0 / 3.0)
        self.assertEqual(result["median"], 2.0)
        self.assertEqual(result["minimum"], 1.0)
        self.assertEqual(result["maximum"], 5.0)

    def test_protein_and_ligand_components_remain_separate(self) -> None:
        flashbind = [
            {
                "protein_global_ca_rmsd_A": "0.2",
                "protein_pocket_ca_rmsd_A": "0.1",
                "ligand_rmsd_pocket_A": "0.6",
            },
            {
                "protein_global_ca_rmsd_A": "0.4",
                "protein_pocket_ca_rmsd_A": "0.3",
                "ligand_rmsd_pocket_A": "1.0",
            },
        ]
        boltz2 = [
            {
                "protein_global_ca_rmsd_A": "1.0",
                "protein_pocket_ca_rmsd_A": "0.5",
                "ligand_rmsd_pocket_A": "0.8",
            },
            {
                "protein_global_ca_rmsd_A": "2.0",
                "protein_pocket_ca_rmsd_A": "1.5",
                "ligand_rmsd_pocket_A": "1.2",
            },
        ]

        result = summarize_rmsd_components(flashbind, boltz2)

        self.assertEqual(set(result), {
            "protein_global_ca",
            "protein_pocket_ca",
            "ligand_after_pocket_alignment",
        })
        self.assertEqual(
            result["protein_global_ca"]["boltz2_msa1024_seed42"]["median"],
            1.5,
        )
        self.assertEqual(
            result["protein_pocket_ca"]["flashbind_released_FABind_plus"]["median"],
            0.2,
        )
        self.assertEqual(
            result["ligand_after_pocket_alignment"]["boltz2_msa1024_seed42"]["median"],
            1.0,
        )

    def test_describe_rejects_empty_values(self) -> None:
        with self.assertRaises(ValueError):
            describe(np.array([]))


if __name__ == "__main__":
    unittest.main()
