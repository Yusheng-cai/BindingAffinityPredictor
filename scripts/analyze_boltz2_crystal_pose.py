#!/usr/bin/env python3
"""Compare exp007 Boltz-2 poses with the 16 exact FEP+4 cocrystals."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pickle
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
from Bio.Align import PairwiseAligner
from rdkit import Chem

from affinity_benchmark.metrics.pose import kabsch, minimum_mapped_rmsd


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPOSITORY_ROOT / "data/manifests/fepplus4_crystal16.json"
DEFAULT_PARENT = REPOSITORY_ROOT / "data/manifests/fepplus4_87.json"
DEFAULT_RAW = REPOSITORY_ROOT / "data/raw/exp008_boltz2_crystal_pose"
DEFAULT_RUN = REPOSITORY_ROOT / "runs/exp008_boltz2_crystal_pose/seed42"
PREDICTION_ROOT = (
    REPOSITORY_ROOT
    / "runs/exp007_fepplus4_boltz2_nesso1/boltz2_msa1024/seed42"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atom_position(atom: gemmi.Atom) -> np.ndarray:
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)


def best_atoms(residue: gemmi.Residue, *, heavy_only: bool = False) -> dict[str, gemmi.Atom]:
    """Choose the highest-occupancy alternate location for each atom name."""

    selected: dict[str, gemmi.Atom] = {}
    for atom in residue:
        if heavy_only and atom.element.atomic_number == 1:
            continue
        name = atom.name.strip()
        current = selected.get(name)
        if current is None or float(atom.occ) > float(current.occ):
            selected[name] = atom
    return selected


def ca_atom(residue: gemmi.Residue) -> gemmi.Atom | None:
    return best_atoms(residue).get("CA")


def one_letter(residue: gemmi.Residue) -> str:
    info = gemmi.find_tabulated_residue(residue.name)
    code = info.one_letter_code
    return code if code and code != " " else "X"


def observed_chain(chain: gemmi.Chain) -> tuple[str, list[gemmi.Residue]]:
    residues = [residue for residue in chain if ca_atom(residue) is not None]
    return "".join(one_letter(residue) for residue in residues), residues


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    return aligner


def aligned_index_pairs(first: str, second: str) -> tuple[list[tuple[int, int]], float]:
    alignment = make_aligner().align(first, second)[0]
    pairs: list[tuple[int, int]] = []
    for first_block, second_block in zip(alignment.aligned[0], alignment.aligned[1]):
        first_start, first_end = map(int, first_block)
        second_start, second_end = map(int, second_block)
        if first_end - first_start != second_end - second_start:
            raise ValueError("unexpected unequal ungapped alignment blocks")
        pairs.extend(
            zip(range(first_start, first_end), range(second_start, second_end), strict=True)
        )
    identity = (
        sum(first[i] == second[j] for i, j in pairs) / len(pairs) if pairs else 0.0
    )
    return pairs, identity


def chain_match_score(chain: gemmi.Chain, expected_sequence: str) -> tuple[float, float, int]:
    sequence, _ = observed_chain(chain)
    pairs, identity = aligned_index_pairs(sequence, expected_sequence)
    matched = sum(sequence[i] == expected_sequence[j] for i, j in pairs)
    coverage = matched / len(expected_sequence)
    return coverage, identity, matched


def protein_chains(model: gemmi.Model) -> list[gemmi.Chain]:
    return [chain for chain in model if len(observed_chain(chain)[1]) >= 20]


def ligand_residues(model: gemmi.Model, residue_name: str) -> list[gemmi.Residue]:
    return [
        residue
        for chain in model
        for residue in chain
        if residue.name.strip().upper() == residue_name.upper()
    ]


def minimum_chain_ligand_distance(chain: gemmi.Chain, ligand: gemmi.Residue) -> float:
    ligand_xyz = [atom_position(atom) for atom in best_atoms(ligand, heavy_only=True).values()]
    return min(
        float(np.linalg.norm(atom_position(atom) - ligand_atom))
        for residue in chain
        for atom in best_atoms(residue, heavy_only=True).values()
        for ligand_atom in ligand_xyz
    )


def choose_predicted_pair(
    model: gemmi.Model, expected_sequence: str
) -> tuple[gemmi.Chain, gemmi.Residue, dict[str, Any]]:
    chains = protein_chains(model)
    ligands = ligand_residues(model, "LIG1")
    if not chains or len(ligands) != 1:
        raise ValueError(
            f"predicted model has {len(chains)} protein chains and "
            f"{len(ligands)} LIG1 residues"
        )
    scored = [(chain_match_score(chain, expected_sequence), chain) for chain in chains]
    score, chain = max(scored, key=lambda item: item[0])
    return chain, ligands[0], {
        "chain": chain.name,
        "expected_sequence_coverage": score[0],
        "aligned_identity": score[1],
    }


def choose_reference_pair(
    model: gemmi.Model, expected_sequence: str, ccd_id: str
) -> tuple[gemmi.Chain, gemmi.Residue, dict[str, Any]]:
    chains = protein_chains(model)
    ligands = ligand_residues(model, ccd_id)
    if not chains or not ligands:
        raise ValueError(
            f"reference has {len(chains)} protein chains and "
            f"{len(ligands)} {ccd_id} residues"
        )

    scored_chains = [(chain_match_score(chain, expected_sequence), chain) for chain in chains]
    best_coverage = max(score[0] for score, _ in scored_chains)
    eligible = [
        (score, chain)
        for score, chain in scored_chains
        if score[0] >= max(0.5, best_coverage - 0.02)
    ]
    candidates = [
        (minimum_chain_ligand_distance(chain, ligand), score, chain, ligand)
        for score, chain in eligible
        for ligand in ligands
    ]
    distance, score, chain, ligand = min(candidates, key=lambda item: item[0])
    if distance > 8.0:
        raise ValueError(f"nearest {ccd_id} is {distance:.2f} A from the target-like chain")
    return chain, ligand, {
        "chain": chain.name,
        "ligand_sequence_number": ligand.seqid.num,
        "expected_sequence_coverage": score[0],
        "aligned_identity": score[1],
        "minimum_chain_ligand_distance_angstrom": distance,
    }


def load_ccd_molecule(ccd_cif: Path, ccd_sdf: Path) -> Chem.Mol:
    block = gemmi.cif.read(str(ccd_cif)).sole_block()
    atom_ids = list(block.find_values("_chem_comp_atom.atom_id"))
    symbols = [symbol.upper() for symbol in block.find_values("_chem_comp_atom.type_symbol")]
    molecule = Chem.MolFromMolFile(str(ccd_sdf), removeHs=False, sanitize=True)
    if molecule is None:
        raise ValueError(f"RDKit could not parse {ccd_sdf}")
    if molecule.GetNumAtoms() != len(atom_ids):
        raise ValueError("CCD CIF and SDF atom counts differ")
    for index, (atom_id, symbol) in enumerate(zip(atom_ids, symbols, strict=True)):
        atom = molecule.GetAtomWithIdx(index)
        if atom.GetSymbol().upper() != symbol:
            raise ValueError(f"CCD atom-order mismatch at {atom_id}")
        atom.SetProp("name", atom_id)
    return Chem.RemoveHs(molecule)


def topology_copy(molecule: Chem.Mol) -> Chem.Mol:
    result = Chem.RWMol(Chem.Mol(molecule))
    for atom in result.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
    for bond in result.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    return result.GetMol()


def graph_mappings(predicted: Chem.Mol, reference: Chem.Mol) -> tuple[list[tuple[int, ...]], str]:
    if predicted.GetNumAtoms() != reference.GetNumAtoms():
        raise ValueError(
            "ligand heavy-atom counts differ: "
            f"predicted {predicted.GetNumAtoms()}, reference {reference.GetNumAtoms()}"
        )
    matches = reference.GetSubstructMatches(
        predicted, useChirality=True, uniquify=False, maxMatches=100000
    )
    mode = "bond_order_and_chirality"
    if not matches:
        matches = topology_copy(reference).GetSubstructMatches(
            topology_copy(predicted), useChirality=True, uniquify=False, maxMatches=100000
        )
        mode = "element_connectivity_and_chirality"
    full = [tuple(match) for match in matches if len(match) == predicted.GetNumAtoms()]
    if not full:
        raise ValueError("no complete chirality-compatible ligand graph mapping")
    return full, mode


def molecule_coordinates(molecule: Chem.Mol, residue: gemmi.Residue) -> np.ndarray:
    coordinates = best_atoms(residue, heavy_only=True)
    missing = [
        molecule.GetAtomWithIdx(index).GetProp("name")
        for index in range(molecule.GetNumAtoms())
        if molecule.GetAtomWithIdx(index).GetProp("name") not in coordinates
    ]
    if missing:
        raise ValueError(f"crystal residue is missing ligand atoms: {missing}")
    return np.array(
        [
            atom_position(coordinates[molecule.GetAtomWithIdx(index).GetProp("name")])
            for index in range(molecule.GetNumAtoms())
        ]
    )


def pocket_residue_ids(
    chain: gemmi.Chain, ligand: gemmi.Residue, cutoff: float
) -> set[int]:
    ligand_xyz = np.array(
        [atom_position(atom) for atom in best_atoms(ligand, heavy_only=True).values()]
    )
    result: set[int] = set()
    _, residues = observed_chain(chain)
    for index, residue in enumerate(residues):
        protein_xyz = np.array(
            [atom_position(atom) for atom in best_atoms(residue, heavy_only=True).values()]
        )
        distances = np.linalg.norm(
            protein_xyz[:, None] - ligand_xyz[None], axis=2
        )
        if protein_xyz.size and np.min(distances) <= cutoff:
            result.add(index)
    return result


def prediction_paths(target_id: str, sample_id: str) -> tuple[Path, Path, Path]:
    root = PREDICTION_ROOT / target_id / "raw" / f"boltz_results_{target_id}"
    prediction_dir = root / "predictions" / sample_id
    return (
        prediction_dir / f"{sample_id}_model_0.cif",
        root / "processed" / "mols" / f"{sample_id}.pkl",
        prediction_dir / f"confidence_{sample_id}_model_0.json",
    )


def analyze_sample(
    sample: dict[str, str], parent_sample: dict[str, Any], raw_root: Path
) -> dict[str, Any]:
    sample_id = sample["sample_id"]
    target_id = sample["target_id"]
    pdb_id = sample["pdb_id"].upper()
    ccd_id = sample["ccd_id"].upper()
    expected_sequence = parent_sample["protein"]["chains"][0]["sequence"]
    prediction_path, mol_path, confidence_path = prediction_paths(target_id, sample_id)
    reference_path = raw_root / "references" / f"{pdb_id.lower()}.cif"
    ccd_cif = raw_root / "ccd" / f"{ccd_id}.cif"
    ccd_sdf = raw_root / "ccd" / f"{ccd_id}_ideal.sdf"
    for path in (prediction_path, mol_path, reference_path, ccd_cif, ccd_sdf):
        if not path.exists():
            raise FileNotFoundError(path)

    predicted_structure = gemmi.read_structure(str(prediction_path))
    reference_structure = gemmi.read_structure(str(reference_path))
    predicted_chain, predicted_ligand, predicted_selection = choose_predicted_pair(
        predicted_structure[0], expected_sequence
    )
    reference_chain, reference_ligand, reference_selection = choose_reference_pair(
        reference_structure[0], expected_sequence, ccd_id
    )

    predicted_sequence, predicted_residues = observed_chain(predicted_chain)
    reference_sequence, reference_residues = observed_chain(reference_chain)
    residue_pairs, pair_identity = aligned_index_pairs(predicted_sequence, reference_sequence)
    identical_pairs = [
        (pred_index, ref_index)
        for pred_index, ref_index in residue_pairs
        if predicted_sequence[pred_index] == reference_sequence[ref_index]
    ]
    if len(identical_pairs) < 50:
        raise ValueError(f"only {len(identical_pairs)} identical CA residue pairs")

    predicted_ca = np.array(
        [atom_position(ca_atom(predicted_residues[i])) for i, _ in identical_pairs]
    )
    reference_ca = np.array(
        [atom_position(ca_atom(reference_residues[j])) for _, j in identical_pairs]
    )
    global_rotation, global_translation, global_ca_rmsd = kabsch(predicted_ca, reference_ca)

    reference_pocket = pocket_residue_ids(reference_chain, reference_ligand, cutoff=5.0)
    pocket_pairs = [(i, j) for i, j in identical_pairs if j in reference_pocket]
    if len(pocket_pairs) < 3:
        raise ValueError(f"only {len(pocket_pairs)} mapped pocket CA atoms")
    predicted_pocket_ca = np.array(
        [atom_position(ca_atom(predicted_residues[i])) for i, _ in pocket_pairs]
    )
    reference_pocket_ca = np.array(
        [atom_position(ca_atom(reference_residues[j])) for _, j in pocket_pairs]
    )
    pocket_rotation, pocket_translation, pocket_ca_rmsd = kabsch(
        predicted_pocket_ca, reference_pocket_ca
    )

    with mol_path.open("rb") as handle:
        predicted_molecules = pickle.load(handle)  # noqa: S301 - trusted local Boltz artifact
    predicted_mol = predicted_molecules["LIG1"]
    reference_mol = load_ccd_molecule(ccd_cif, ccd_sdf)
    mappings, mapping_mode = graph_mappings(predicted_mol, reference_mol)
    predicted_ligand_xyz = molecule_coordinates(predicted_mol, predicted_ligand)
    reference_ligand_xyz = molecule_coordinates(reference_mol, reference_ligand)

    global_ligand_rmsd, global_mapping = minimum_mapped_rmsd(
        predicted_ligand_xyz,
        reference_ligand_xyz,
        mappings,
        rotation=global_rotation,
        translation=global_translation,
    )
    pocket_ligand_rmsd, pocket_mapping = minimum_mapped_rmsd(
        predicted_ligand_xyz,
        reference_ligand_xyz,
        mappings,
        rotation=pocket_rotation,
        translation=pocket_translation,
    )
    internal_ligand_rmsd, _ = minimum_mapped_rmsd(
        predicted_ligand_xyz, reference_ligand_xyz, mappings, fit_mobile=True
    )

    reference_contacts = pocket_residue_ids(reference_chain, reference_ligand, cutoff=4.0)
    predicted_contacts = pocket_residue_ids(predicted_chain, predicted_ligand, cutoff=4.0)
    mapped_reference_to_predicted = {
        ref_index: pred_index for pred_index, ref_index in identical_pairs
    }
    comparable_reference = set(mapped_reference_to_predicted)
    reference_contact_mapped = reference_contacts & comparable_reference
    predicted_contact_as_reference = {
        ref_index
        for ref_index, pred_index in mapped_reference_to_predicted.items()
        if pred_index in predicted_contacts
    }
    intersection = reference_contact_mapped & predicted_contact_as_reference
    union = reference_contact_mapped | predicted_contact_as_reference

    confidence = json.loads(confidence_path.read_text()) if confidence_path.exists() else {}
    return {
        "sample_id": sample_id,
        "target_id": target_id,
        "pdb_id": pdb_id,
        "ccd_id": ccd_id,
        "status": "complete",
        "prediction_path": str(prediction_path.relative_to(REPOSITORY_ROOT)),
        "reference_path": str(reference_path.relative_to(REPOSITORY_ROOT)),
        "prediction_sha256": sha256(prediction_path),
        "reference_sha256": sha256(reference_path),
        "predicted_selection": predicted_selection,
        "reference_selection": reference_selection,
        "protein_alignment": {
            "sequence_mapped_identical_ca_atoms": len(identical_pairs),
            "pairwise_observed_sequence_identity": pair_identity,
            "global_ca_rmsd_angstrom": global_ca_rmsd,
            "reference_pocket_residues_5A": len(reference_pocket),
            "mapped_pocket_ca_atoms": len(pocket_pairs),
            "pocket_ca_rmsd_angstrom": pocket_ca_rmsd,
        },
        "ligand_pose": {
            "heavy_atoms": predicted_mol.GetNumAtoms(),
            "graph_mapping_mode": mapping_mode,
            "symmetry_compatible_mappings": len(mappings),
            "global_aligned_rmsd_angstrom": global_ligand_rmsd,
            "pocket_aligned_rmsd_angstrom": pocket_ligand_rmsd,
            "internal_conformation_rmsd_angstrom": internal_ligand_rmsd,
            "global_best_mapping": list(global_mapping),
            "pocket_best_mapping": list(pocket_mapping),
        },
        "contacts_4A": {
            "reference_count": len(reference_contact_mapped),
            "predicted_count": len(predicted_contact_as_reference),
            "recovered_count": len(intersection),
            "precision": len(intersection) / len(predicted_contact_as_reference)
            if predicted_contact_as_reference
            else None,
            "recall": len(intersection) / len(reference_contact_mapped)
            if reference_contact_mapped
            else None,
            "jaccard": len(intersection) / len(union) if union else None,
        },
        "boltz_confidence": confidence,
    }


def compact_row(result: dict[str, Any]) -> dict[str, Any]:
    if result["status"] != "complete":
        return {
            "sample_id": result["sample_id"],
            "target_id": result["target_id"],
            "pdb_id": result["pdb_id"],
            "status": "failed",
            "failure": result["failure"],
        }
    protein = result["protein_alignment"]
    ligand = result["ligand_pose"]
    contacts = result["contacts_4A"]
    confidence = result["boltz_confidence"]
    return {
        "sample_id": result["sample_id"],
        "target_id": result["target_id"],
        "pdb_id": result["pdb_id"],
        "ccd_id": result["ccd_id"],
        "status": "complete",
        "ligand_rmsd_pocket_A": ligand["pocket_aligned_rmsd_angstrom"],
        "ligand_rmsd_global_A": ligand["global_aligned_rmsd_angstrom"],
        "ligand_internal_rmsd_A": ligand["internal_conformation_rmsd_angstrom"],
        "protein_pocket_ca_rmsd_A": protein["pocket_ca_rmsd_angstrom"],
        "protein_global_ca_rmsd_A": protein["global_ca_rmsd_angstrom"],
        "contact_precision_4A": contacts["precision"],
        "contact_recall_4A": contacts["recall"],
        "contact_jaccard_4A": contacts["jaccard"],
        "boltz_confidence_score": confidence.get("confidence_score"),
        "boltz_complex_iplddt": confidence.get("complex_iplddt"),
    }


def summarize(results: list[dict[str, Any]]) -> dict[str, Any]:
    complete = [result for result in results if result["status"] == "complete"]
    pocket = np.array(
        [result["ligand_pose"]["pocket_aligned_rmsd_angstrom"] for result in complete]
    )
    global_values = np.array(
        [result["ligand_pose"]["global_aligned_rmsd_angstrom"] for result in complete]
    )
    by_target: dict[str, dict[str, Any]] = {}
    for target, target_results in sorted(
        (
            (target, [r for r in complete if r["target_id"] == target])
            for target in {r["target_id"] for r in complete}
        )
    ):
        values = np.array(
            [
                r["ligand_pose"]["pocket_aligned_rmsd_angstrom"]
                for r in target_results
            ]
        )
        by_target[target] = {
            "n": len(values),
            "median_pocket_aligned_ligand_rmsd_angstrom": float(np.median(values)),
            "success_below_2A": int(np.sum(values < 2.0)),
        }
    return {
        "requested": len(results),
        "complete": len(complete),
        "failed": len(results) - len(complete),
        "pocket_aligned_ligand_rmsd_angstrom": {
            "mean": float(np.mean(pocket)) if len(pocket) else None,
            "median": float(np.median(pocket)) if len(pocket) else None,
            "minimum": float(np.min(pocket)) if len(pocket) else None,
            "maximum": float(np.max(pocket)) if len(pocket) else None,
        },
        "global_aligned_ligand_rmsd_angstrom": {
            "mean": float(np.mean(global_values)) if len(global_values) else None,
            "median": float(np.median(global_values)) if len(global_values) else None,
        },
        "pose_success": {
            "threshold_angstrom": 2.0,
            "count": int(np.sum(pocket < 2.0)),
            "fraction_of_complete": float(np.mean(pocket < 2.0)) if len(pocket) else None,
        },
        "by_target": by_target,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--parent-manifest", type=Path, default=DEFAULT_PARENT)
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RUN)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    parent = json.loads(args.parent_manifest.read_text())
    parent_by_id = {sample["sample_id"]: sample for sample in parent["samples"]}
    results: list[dict[str, Any]] = []
    for sample in manifest["samples"]:
        print(f"Analyzing {sample['sample_id']} against {sample['pdb_id']}...")
        try:
            result = analyze_sample(sample, parent_by_id[sample["sample_id"]], args.raw_root)
        except Exception as error:  # preserve every per-pair failure in the run record
            result = {
                "sample_id": sample["sample_id"],
                "target_id": sample["target_id"],
                "pdb_id": sample["pdb_id"],
                "ccd_id": sample["ccd_id"],
                "status": "failed",
                "failure": f"{type(error).__name__}: {error}",
            }
        results.append(result)

    summary = summarize(results)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPOSITORY_ROOT, text=True
    ).strip()
    report = {
        "experiment_id": "exp008_boltz2_crystal_pose",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_git_revision": revision,
        "manifest_id": manifest["manifest_id"],
        "prediction_protocol": "boltz2_msa1024_seed42_one_structure_sample",
        "reference_coordinates_used_during_inference": False,
        "method": {
            "protein_mapping": (
                "global sequence alignment; identical observed residues with CA atoms"
            ),
            "global_fit": "all mapped kinase CA atoms",
            "pocket_fit": (
                "mapped kinase CA atoms whose reference residue has a heavy atom "
                "within 5 A of ligand"
            ),
            "ligand_metric": (
                "minimum heavy-atom RMSD over chirality-compatible graph mappings "
                "after fixed protein transform"
            ),
            "ligand_atoms_in_fit": False,
        },
        "summary": summary,
        "samples": results,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "pose_metrics.json"
    csv_path = args.output_dir / "pose_metrics.csv"
    json_path.write_text(json.dumps(report, indent=2) + "\n")
    rows = [compact_row(result) for result in results]
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2))
    print(f"Wrote {json_path} and {csv_path}")


if __name__ == "__main__":
    main()
