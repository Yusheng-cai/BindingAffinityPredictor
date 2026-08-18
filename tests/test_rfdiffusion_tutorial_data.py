"""Consistency checks for the tracked seed-42 browser data artifact."""

from __future__ import annotations

import json
import unittest
from collections import Counter, defaultdict
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = REPOSITORY_ROOT / "docs/rfdiffusion_binder_tutorial/real-run-data.js"
PREFIX = "window.RFDIFFUSION_REAL_RUN="


class RFdiffusionTutorialDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        source = DATA_PATH.read_text()
        assignment = next(line for line in source.splitlines() if line.startswith(PREFIX))
        cls.data = json.loads(assignment[len(PREFIX) :].removesuffix(";"))

    def test_selected_sequence_maps_one_to_one_to_binder_residues(self) -> None:
        binder = self.data["binder"]
        selected_sample = self.data["proteinMpnn"]["selectedSample"]
        selected = next(
            candidate
            for candidate in self.data["proteinMpnn"]["candidates"]
            if candidate["sample"] == selected_sample
        )
        exported_sequence = "".join(
            residue["oneLetter"] for residue in binder["assignedResidues"]
        )
        self.assertEqual(exported_sequence, selected["sequence"])
        self.assertEqual(len(exported_sequence), binder["length"])

    def test_each_binder_residue_has_exactly_four_backbone_atoms(self) -> None:
        atoms_by_residue: dict[int, list[str]] = defaultdict(list)
        for atom in self.data["binder"]["backboneAtoms"]:
            atoms_by_residue[atom["residue"]].append(atom["atom"])
        self.assertEqual(len(atoms_by_residue), 90)
        for residue, atom_names in atoms_by_residue.items():
            self.assertEqual(Counter(atom_names), Counter(["N", "CA", "C", "O"]), residue)

    def test_coordinate_completeness_is_not_overstated(self) -> None:
        binder = self.data["binder"]
        self.assertEqual(binder["coordinateCompleteness"], "backbone_N_CA_C_O_only")
        self.assertFalse(binder["sidechainsPacked"])


if __name__ == "__main__":
    unittest.main()
