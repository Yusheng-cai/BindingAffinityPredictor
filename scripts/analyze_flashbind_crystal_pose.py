#!/usr/bin/env python3
"""Compare released FlashBind/FABind+ poses with 16 exact FEP+4 cocrystals."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi
import lmdb
import numpy as np
from rdkit import Chem

from affinity_benchmark.adapters.affinity_models import flashbind_record_id
from affinity_benchmark.metrics.pose import kabsch, minimum_mapped_rmsd
from analyze_boltz2_crystal_pose import (
    aligned_index_pairs,
    atom_position,
    best_atoms,
    ca_atom,
    chain_match_score,
    choose_reference_pair,
    graph_mappings,
    load_ccd_molecule,
    molecule_coordinates,
    observed_chain,
    pocket_residue_ids,
    protein_chains,
    sha256,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPOSITORY_ROOT / "data/manifests/fepplus4_crystal16.json"
DEFAULT_PARENT = REPOSITORY_ROOT / "data/manifests/fepplus4_87.json"
DEFAULT_RAW = REPOSITORY_ROOT / "data/raw/exp008_boltz2_crystal_pose"
DEFAULT_FLASHBIND_DATA = REPOSITORY_ROOT / "cache/flashbind/datasets/fep4/fep4"
DEFAULT_RUN = REPOSITORY_ROOT / "runs/exp010_flashbind_crystal_pose/released_poses"


def bytes_sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def choose_predicted_chain(
    model: gemmi.Model, expected_sequence: str
) -> tuple[gemmi.Chain, dict[str, Any]]:
    chains = protein_chains(model)
    if not chains:
        raise ValueError("released receptor has no protein chain with at least 20 CA atoms")
    candidates: list[tuple[gemmi.Chain, dict[str, Any]]] = [
        (chain, {"chain": chain.name, "chain_id_repair": None}) for chain in chains
    ]
    # The released p38 PDB loses its chain ID after residue 168. Gemmi correctly
    # parses this as two chains, although the residue numbering and target
    # sequence show that they are fragments of one kinase construct. Repair only
    # this narrow formatting case; do not concatenate ordinary biological chains.
    unnamed = [chain for chain in chains if chain.name == ""]
    named = [chain for chain in chains if chain.name != ""]
    if len(unnamed) == 1 and len(named) == 1:
        joined = gemmi.Chain(named[0].name)
        residues = sorted(
            (residue for chain in (named[0], unnamed[0]) for residue in chain),
            key=lambda residue: (residue.seqid.num, residue.seqid.icode),
        )
        for residue in residues:
            joined.add_residue(residue.clone())
        candidates.append(
            (
                joined,
                {
                    "chain": named[0].name,
                    "chain_id_repair": (
                        "joined named and blank-ID fragments by deposited residue number"
                    ),
                    "source_chain_ids": [named[0].name, unnamed[0].name],
                },
            )
        )
    scored = [
        (chain_match_score(chain, expected_sequence), chain, metadata)
        for chain, metadata in candidates
    ]
    score, chain, metadata = max(scored, key=lambda item: item[0])
    return chain, {
        **metadata,
        "expected_sequence_coverage": score[0],
        "aligned_identity": score[1],
    }


def molecule_from_sdf_bytes(content: bytes) -> Chem.Mol:
    molecule = Chem.MolFromMolBlock(
        content.decode("utf-8"), removeHs=False, sanitize=True
    )
    if molecule is None:
        raise ValueError("RDKit could not parse released FABind+ SDF record")
    molecule = Chem.RemoveHs(molecule)
    if molecule.GetNumConformers() != 1:
        raise ValueError(
            f"released FABind+ molecule has {molecule.GetNumConformers()} conformers"
        )
    return molecule


def conformer_coordinates(molecule: Chem.Mol) -> np.ndarray:
    return np.asarray(molecule.GetConformer().GetPositions(), dtype=float)


def coordinate_contacts(
    chain: gemmi.Chain, ligand_xyz: np.ndarray, cutoff: float
) -> set[int]:
    result: set[int] = set()
    _, residues = observed_chain(chain)
    for index, residue in enumerate(residues):
        protein_xyz = np.array(
            [atom_position(atom) for atom in best_atoms(residue, heavy_only=True).values()]
        )
        if protein_xyz.size == 0:
            continue
        distances = np.linalg.norm(protein_xyz[:, None] - ligand_xyz[None], axis=2)
        if float(np.min(distances)) <= cutoff:
            result.add(index)
    return result


def analyze_sample(
    sample: dict[str, str],
    parent_sample: dict[str, Any],
    raw_root: Path,
    flashbind_root: Path,
    ligand_environment: lmdb.Environment,
) -> dict[str, Any]:
    sample_id = sample["sample_id"]
    target_id = sample["target_id"]
    pdb_id = sample["pdb_id"].upper()
    ccd_id = sample["ccd_id"].upper()
    expected_sequence = parent_sample["protein"]["chains"][0]["sequence"]
    record_id = flashbind_record_id(parent_sample)
    receptor_path = flashbind_root / "pdb" / f"{target_id}.pdb"
    reference_path = raw_root / "references" / f"{pdb_id.lower()}.cif"
    ccd_cif = raw_root / "ccd" / f"{ccd_id}.cif"
    ccd_sdf = raw_root / "ccd" / f"{ccd_id}_ideal.sdf"
    for path in (receptor_path, reference_path, ccd_cif, ccd_sdf):
        if not path.exists():
            raise FileNotFoundError(path)
    with ligand_environment.begin() as transaction:
        sdf_content = transaction.get(record_id.encode("utf-8"))
    if sdf_content is None:
        raise KeyError(f"released ligand pose {record_id!r} is absent from LMDB")

    predicted_structure = gemmi.read_structure(str(receptor_path))
    reference_structure = gemmi.read_structure(str(reference_path))
    predicted_chain, predicted_selection = choose_predicted_chain(
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
    global_rotation, global_translation, global_ca_rmsd = kabsch(
        predicted_ca, reference_ca
    )

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

    predicted_mol = molecule_from_sdf_bytes(sdf_content)
    reference_mol = load_ccd_molecule(ccd_cif, ccd_sdf)
    mappings, mapping_mode = graph_mappings(predicted_mol, reference_mol)
    predicted_ligand_xyz = conformer_coordinates(predicted_mol)
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
    predicted_contacts = coordinate_contacts(predicted_chain, predicted_ligand_xyz, cutoff=4.0)
    mapped_reference_to_predicted = {
        ref_index: pred_index for pred_index, ref_index in identical_pairs
    }
    reference_contact_mapped = reference_contacts & set(mapped_reference_to_predicted)
    predicted_contact_as_reference = {
        ref_index
        for ref_index, pred_index in mapped_reference_to_predicted.items()
        if pred_index in predicted_contacts
    }
    intersection = reference_contact_mapped & predicted_contact_as_reference
    union = reference_contact_mapped | predicted_contact_as_reference

    return {
        "sample_id": sample_id,
        "target_id": target_id,
        "pdb_id": pdb_id,
        "ccd_id": ccd_id,
        "flashbind_record_id": record_id,
        "status": "complete",
        "receptor_path": str(receptor_path.relative_to(REPOSITORY_ROOT)),
        "reference_path": str(reference_path.relative_to(REPOSITORY_ROOT)),
        "receptor_sha256": sha256(receptor_path),
        "released_sdf_record_sha256": bytes_sha256(sdf_content),
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
    return {
        "sample_id": result["sample_id"],
        "target_id": result["target_id"],
        "pdb_id": result["pdb_id"],
        "ccd_id": result["ccd_id"],
        "flashbind_record_id": result["flashbind_record_id"],
        "status": "complete",
        "ligand_rmsd_pocket_A": ligand["pocket_aligned_rmsd_angstrom"],
        "ligand_rmsd_global_A": ligand["global_aligned_rmsd_angstrom"],
        "ligand_internal_rmsd_A": ligand["internal_conformation_rmsd_angstrom"],
        "protein_pocket_ca_rmsd_A": protein["pocket_ca_rmsd_angstrom"],
        "protein_global_ca_rmsd_A": protein["global_ca_rmsd_angstrom"],
        "contact_precision_4A": contacts["precision"],
        "contact_recall_4A": contacts["recall"],
        "contact_jaccard_4A": contacts["jaccard"],
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
    target_groups: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for result in complete:
        target_groups[result["target_id"]].append(result)
    for target, target_results in sorted(target_groups.items()):
        values = np.array(
            [r["ligand_pose"]["pocket_aligned_rmsd_angstrom"] for r in target_results]
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
    parser.add_argument("--flashbind-root", type=Path, default=DEFAULT_FLASHBIND_DATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RUN)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    parent = json.loads(args.parent_manifest.read_text())
    parent_by_id = {sample["sample_id"]: sample for sample in parent["samples"]}
    ligand_lmdb = args.flashbind_root / "ligand_sdf.lmdb"
    if not ligand_lmdb.exists():
        raise FileNotFoundError(ligand_lmdb)

    results: list[dict[str, Any]] = []
    ligand_environment = lmdb.open(str(ligand_lmdb), readonly=True, lock=False)
    try:
        for sample in manifest["samples"]:
            print(f"Analyzing {sample['sample_id']} against {sample['pdb_id']}...")
            try:
                result = analyze_sample(
                    sample,
                    parent_by_id[sample["sample_id"]],
                    args.raw_root,
                    args.flashbind_root,
                    ligand_environment,
                )
            except Exception as error:  # preserve every per-pair failure
                result = {
                    "sample_id": sample["sample_id"],
                    "target_id": sample["target_id"],
                    "pdb_id": sample["pdb_id"],
                    "ccd_id": sample["ccd_id"],
                    "status": "failed",
                    "failure": f"{type(error).__name__}: {error}",
                }
            results.append(result)
    finally:
        ligand_environment.close()

    summary = summarize(results)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPOSITORY_ROOT, text=True
    ).strip()
    report = {
        "experiment_id": "exp010_flashbind_crystal_pose",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_git_revision": revision,
        "manifest_id": manifest["manifest_id"],
        "prediction_protocol": "authors_released_FABind_plus_pose_archive",
        "released_pose_generation_seed": "unknown_not_provided",
        "reference_coordinates_used_during_prediction": False,
        "model_source_revision": "f161268176237ab6ce5757031a8c1b93937e0d37",
        "released_data_revision": "50b1511e080236d80f9b8ff1e4d0cd38bb2480b9",
        "released_ligand_lmdb_sha256": sha256(ligand_lmdb / "data.mdb"),
        "method": {
            "protein_mapping": "global sequence alignment; identical observed residues with CA atoms",
            "global_fit": "all mapped kinase CA atoms",
            "pocket_fit": "mapped kinase CA atoms whose reference residue has a heavy atom within 5 A of ligand",
            "ligand_metric": "minimum heavy-atom RMSD over chirality-compatible graph mappings after fixed protein transform",
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
