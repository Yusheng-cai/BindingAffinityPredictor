#!/usr/bin/env python3
"""Compare one Boltz protein-ligand prediction with a withheld mmCIF reference."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
from rdkit import Chem


def load_structure(path: Path) -> gemmi.Structure:
    return gemmi.make_structure_from_block(gemmi.cif.read(str(path)).sole_block())


def select_protein_chain(model: gemmi.Model) -> gemmi.Chain:
    candidates = [chain for chain in model if sum("CA" in residue for residue in chain) >= 20]
    if len(candidates) != 1:
        raise ValueError(f"expected one protein chain; found {len(candidates)}")
    return candidates[0]


def select_ligand(
    model: gemmi.Model, residue_name: str, sequence_number: int | None = None
) -> tuple[gemmi.Chain, gemmi.Residue]:
    candidates = [
        (chain, residue)
        for chain in model
        for residue in chain
        if residue.name == residue_name
        and (sequence_number is None or residue.seqid.num == sequence_number)
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"expected one residue named {residue_name}; found {len(candidates)}"
        )
    return candidates[0]


def residue_index(residue: gemmi.Residue) -> int:
    """Use label_seq when present so reference numbering matches the input sequence."""

    return residue.label_seq if residue.label_seq is not None else residue.seqid.num


def atom_position(atom: gemmi.Atom) -> np.ndarray:
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)


def paired_ca(
    predicted: gemmi.Chain, reference: gemmi.Chain, indices: set[int] | None = None
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    pred = {
        residue_index(residue): atom_position(residue["CA"][0])
        for residue in predicted
        if "CA" in residue
    }
    ref = {
        residue_index(residue): atom_position(residue["CA"][0])
        for residue in reference
        if "CA" in residue
    }
    common = sorted(set(pred) & set(ref))
    if indices is not None:
        common = [index for index in common if index in indices]
    return np.array([pred[index] for index in common]), np.array([ref[index] for index in common]), common


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return rotation, translation, and fitted RMSD mapping mobile to target."""

    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (target - target_center)
    left, _, right_t = np.linalg.svd(covariance)
    rotation = right_t.T @ left.T
    if np.linalg.det(rotation) < 0:
        right_t[-1] *= -1
        rotation = right_t.T @ left.T
    translation = target_center - rotation @ mobile_center
    aligned = (rotation @ mobile.T).T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def transform(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return (rotation @ coords.T).T + translation


def heavy_atom_dict(residue: gemmi.Residue) -> dict[str, np.ndarray]:
    return {
        atom.name: atom_position(atom)
        for atom in residue
        if atom.element.name != "H"
    }


def contact_indices(
    protein: gemmi.Chain, ligand: gemmi.Residue, cutoff: float
) -> dict[int, dict[str, Any]]:
    ligand_coords = list(heavy_atom_dict(ligand).values())
    result: dict[int, dict[str, Any]] = {}
    for residue in protein:
        distances = [
            float(np.linalg.norm(atom_position(atom) - ligand_coord))
            for atom in residue
            if atom.element.name != "H"
            for ligand_coord in ligand_coords
        ]
        if distances and min(distances) <= cutoff:
            index = residue_index(residue)
            result[index] = {
                "residue_name": residue.name,
                "sequence_index": index,
                "minimum_heavy_atom_distance_angstrom": min(distances),
            }
    return result


def rmsd(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((first - second) ** 2, axis=1))))


def chiral_centers_from_coords(molecule: Chem.Mol, coords: np.ndarray) -> list[list[Any]]:
    """Assign CIP labels from coordinates on a copy of the molecular graph."""

    molecule = Chem.Mol(molecule)
    molecule.RemoveAllConformers()
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, xyz in enumerate(coords):
        conformer.SetAtomPosition(index, xyz)
    molecule.AddConformer(conformer)
    for atom in molecule.GetAtoms():
        if atom.HasProp("_CIPCode"):
            atom.ClearProp("_CIPCode")
    Chem.AssignAtomChiralTagsFromStructure(
        molecule, confId=0, replaceExistingTags=True
    )
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return [list(center) for center in Chem.FindMolChiralCenters(molecule, includeUnassigned=True)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--predicted-mol", type=Path, required=True)
    parser.add_argument("--reference-mol", type=Path, required=True)
    parser.add_argument("--predicted-ligand-name", default="LIG1")
    parser.add_argument("--reference-ligand-name", default="GAL")
    parser.add_argument("--reference-ligand-seqid", type=int)
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--pocket-cutoff", type=float, default=5.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    predicted_structure = load_structure(args.prediction)
    reference_structure = load_structure(args.reference)
    predicted_model = predicted_structure[0]
    reference_model = reference_structure[0]
    predicted_protein = select_protein_chain(predicted_model)
    reference_protein = select_protein_chain(reference_model)
    _, predicted_ligand = select_ligand(predicted_model, args.predicted_ligand_name)
    _, reference_ligand = select_ligand(
        reference_model, args.reference_ligand_name, args.reference_ligand_seqid
    )

    pred_ca, ref_ca, common_indices = paired_ca(predicted_protein, reference_protein)
    global_rotation, global_translation, global_ca_rmsd = kabsch(pred_ca, ref_ca)

    reference_contacts_4 = contact_indices(
        reference_protein, reference_ligand, args.contact_cutoff
    )
    predicted_contacts_4 = contact_indices(
        predicted_protein, predicted_ligand, args.contact_cutoff
    )
    reference_pocket = contact_indices(
        reference_protein, reference_ligand, args.pocket_cutoff
    )
    pocket_indices = set(reference_pocket)
    pred_pocket_ca, ref_pocket_ca, fitted_pocket_indices = paired_ca(
        predicted_protein, reference_protein, pocket_indices
    )
    pocket_rotation, pocket_translation, pocket_ca_rmsd = kabsch(
        pred_pocket_ca, ref_pocket_ca
    )

    with args.predicted_mol.open("rb") as handle:
        predicted_molecules = pickle.load(handle)  # noqa: S301 - trusted local Boltz artifact
    predicted_mol = predicted_molecules[args.predicted_ligand_name]
    with args.reference_mol.open("rb") as handle:
        reference_mol = pickle.load(handle)  # noqa: S301 - trusted checksummed Boltz cache

    matches = reference_mol.GetSubstructMatches(
        predicted_mol, useChirality=True, uniquify=False
    )
    if len(matches) != 1:
        raise ValueError(f"expected one chirality-aware atom mapping; found {len(matches)}")
    reference_indices = matches[0]
    predicted_coords_by_name = heavy_atom_dict(predicted_ligand)
    reference_coords_by_name = heavy_atom_dict(reference_ligand)
    predicted_coords = np.array(
        [
            predicted_coords_by_name[predicted_mol.GetAtomWithIdx(index).GetProp("name")]
            for index in range(predicted_mol.GetNumAtoms())
        ]
    )
    reference_coords = np.array(
        [
            reference_coords_by_name[reference_mol.GetAtomWithIdx(index).GetProp("name")]
            for index in reference_indices
        ]
    )
    atom_mapping = [
        {
            "predicted_atom": predicted_mol.GetAtomWithIdx(index).GetProp("name"),
            "reference_atom": reference_mol.GetAtomWithIdx(reference_index).GetProp("name"),
        }
        for index, reference_index in enumerate(reference_indices)
    ]

    global_ligand_rmsd = rmsd(
        transform(predicted_coords, global_rotation, global_translation), reference_coords
    )
    pocket_ligand_rmsd = rmsd(
        transform(predicted_coords, pocket_rotation, pocket_translation), reference_coords
    )
    _, _, internal_ligand_rmsd = kabsch(predicted_coords, reference_coords)

    predicted_bond_lengths = []
    reference_bond_lengths = []
    for bond in predicted_mol.GetBonds():
        first = bond.GetBeginAtomIdx()
        second = bond.GetEndAtomIdx()
        predicted_bond_lengths.append(
            float(np.linalg.norm(predicted_coords[first] - predicted_coords[second]))
        )
        reference_bond_lengths.append(
            float(np.linalg.norm(reference_coords[first] - reference_coords[second]))
        )
    bond_length_differences = np.array(predicted_bond_lengths) - np.array(
        reference_bond_lengths
    )
    predicted_chiral_centers = chiral_centers_from_coords(
        predicted_mol, predicted_coords
    )
    reference_chiral_centers = chiral_centers_from_coords(
        predicted_mol, reference_coords
    )

    reference_contact_set = set(reference_contacts_4)
    predicted_contact_set = set(predicted_contacts_4)
    intersection = reference_contact_set & predicted_contact_set
    union = reference_contact_set | predicted_contact_set

    report = {
        "protein_alignment": {
            "common_ca_atoms": len(common_indices),
            "global_ca_rmsd_angstrom": global_ca_rmsd,
            "reference_pocket_cutoff_angstrom": args.pocket_cutoff,
            "pocket_ca_atoms": len(fitted_pocket_indices),
            "pocket_sequence_indices": fitted_pocket_indices,
            "pocket_ca_rmsd_angstrom": pocket_ca_rmsd,
        },
        "ligand_pose": {
            "heavy_atoms": predicted_mol.GetNumAtoms(),
            "chirality_aware_graph_mappings": len(matches),
            "atom_mapping": atom_mapping,
            "rmsd_after_global_protein_ca_alignment_angstrom": global_ligand_rmsd,
            "rmsd_after_reference_pocket_ca_alignment_angstrom": pocket_ligand_rmsd,
            "internal_conformation_rmsd_after_ligand_only_alignment_angstrom": internal_ligand_rmsd,
        },
        "ligand_geometry": {
            "predicted_bond_length_range_angstrom": [
                min(predicted_bond_lengths),
                max(predicted_bond_lengths),
            ],
            "reference_bond_length_range_angstrom": [
                min(reference_bond_lengths),
                max(reference_bond_lengths),
            ],
            "bond_length_difference_rmse_angstrom": float(
                np.sqrt(np.mean(bond_length_differences**2))
            ),
            "bond_length_difference_max_abs_angstrom": float(
                np.max(np.abs(bond_length_differences))
            ),
            "predicted_coordinate_cip_centers": predicted_chiral_centers,
            "reference_coordinate_cip_centers": reference_chiral_centers,
            "all_coordinate_cip_centers_match": predicted_chiral_centers
            == reference_chiral_centers,
        },
        "contacts": {
            "cutoff_angstrom": args.contact_cutoff,
            "reference": list(reference_contacts_4.values()),
            "predicted": list(predicted_contacts_4.values()),
            "recovered_sequence_indices": sorted(intersection),
            "precision": len(intersection) / len(predicted_contact_set) if predicted_contact_set else None,
            "recall": len(intersection) / len(reference_contact_set) if reference_contact_set else None,
            "jaccard": len(intersection) / len(union) if union else None,
        },
        "method_notes": [
            "Protein residues are matched by mmCIF label_seq/input-sequence index.",
            "Ligand atoms are matched by one unique chirality-aware RDKit graph isomorphism.",
            "The global pose RMSD uses a Kabsch fit over every common protein C-alpha atom.",
            "The pocket pose RMSD uses a Kabsch fit over reference residues within the declared pocket cutoff.",
            "Hydrogens are excluded from ligand RMSD and contact calculations.",
        ],
    }

    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
