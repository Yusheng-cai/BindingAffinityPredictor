#!/usr/bin/env python3
"""Score Boltz-2 structures against withheld Runs N' Poses coordinates.

The model receives sequence, SMILES, and an MSA only. This script is run after
inference and uses the experimental coordinates solely as evaluation labels.
"""

from __future__ import annotations

import argparse
import csv
import json
import pickle
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
from rdkit import Chem

from affinity_benchmark.metrics.distogram import precision_recall_f1
from affinity_benchmark.metrics.pose import kabsch, minimum_mapped_rmsd
from analyze_boltz2_crystal_pose import (
    best_atoms,
    choose_predicted_pair,
    graph_mappings,
    molecule_coordinates,
    one_letter,
)
from analyze_nesso_rnp_distograms import (
    _cluster_bootstrap_interval,
    _physical_pocket_labels,
    _protein_ligand_distances,
    _spearman,
    familiarity_metadata,
    read_reference_protein,
    remove_nesso_tokenization_hydrogens,
)


COMMON_METRICS = (
    "interface_expected_distance_mae_angstrom",
    "token_contact_precision_at_6a",
    "token_contact_recall_at_6a",
    "token_contact_f1_at_6a",
    "physical_pocket_precision_at_6a",
    "physical_pocket_recall_at_6a",
    "physical_pocket_f1_at_6a",
)


def chain_coordinates(
    chain: gemmi.Chain, sequence: str
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray | None]]:
    """Return sequence-indexed C-beta/C-alpha, C-alpha, and heavy coordinates."""

    token = np.full((len(sequence), 3), np.nan, dtype=float)
    ca = np.full((len(sequence), 3), np.nan, dtype=float)
    heavy: list[np.ndarray | None] = [None] * len(sequence)
    for residue in chain:
        index = residue.seqid.num - 1
        if index < 0 or index >= len(sequence):
            continue
        # Gemmi uses lowercase one-letter codes for several modified residues;
        # compare their parent identity case-insensitively.
        observed = one_letter(residue).upper()
        if observed != "X" and observed != sequence[index]:
            raise ValueError(
                f"sequence mismatch at {chain.name}:{index + 1}: "
                f"expected={sequence[index]} observed={observed}"
            )
        atoms = best_atoms(residue, heavy_only=True)
        if "CA" in atoms:
            ca[index] = (atoms["CA"].pos.x, atoms["CA"].pos.y, atoms["CA"].pos.z)
        token_name = "CA" if observed == "G" else "CB"
        if token_name in atoms:
            atom = atoms[token_name]
            token[index] = (atom.pos.x, atom.pos.y, atom.pos.z)
        if atoms:
            heavy[index] = np.asarray(
                [(atom.pos.x, atom.pos.y, atom.pos.z) for atom in atoms.values()],
                dtype=float,
            )
    return token, ca, heavy


def reference_chain(model: gemmi.Model, chain_id: str) -> gemmi.Chain:
    matches = [chain for chain in model if chain.name == chain_id]
    if len(matches) != 1:
        raise ValueError(f"expected one reference chain {chain_id}, found {len(matches)}")
    return matches[0]


def reference_ligand_molecule(path: Path) -> tuple[Chem.Mol, np.ndarray]:
    molecule = Chem.MolFromMolFile(str(path), removeHs=True, sanitize=True)
    if molecule is None:
        raise ValueError(f"RDKit could not parse {path}")
    molecule = Chem.RemoveHs(molecule)
    conformer = molecule.GetConformer()
    coordinates = np.asarray(
        [
            tuple(conformer.GetAtomPosition(index))
            for index in range(molecule.GetNumAtoms())
        ],
        dtype=float,
    )
    return molecule, coordinates


def score_sample(
    sample: dict[str, Any],
    prediction_root: Path,
    ground_truth_root: Path,
    *,
    interface_cutoff: float = 15.0,
    contact_cutoff: float = 6.0,
) -> dict[str, Any]:
    sample_id = sample["sample_id"]
    system_id = sample["runs_n_poses"]["system_id"]
    sequence = sample["protein"]["chains"][0]["sequence"]
    reference_chain_id = sample["protein"]["chains"][0]["id"]

    prediction_dir = prediction_root / "predictions" / sample_id
    prediction_cif = prediction_dir / f"{sample_id}_model_0.cif"
    confidence_path = prediction_dir / f"confidence_{sample_id}_model_0.json"
    molecule_path = prediction_root / "processed" / "mols" / f"{sample_id}.pkl"
    for path in (prediction_cif, molecule_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    predicted_structure = gemmi.read_structure(str(prediction_cif))
    predicted_chain, predicted_ligand, predicted_selection = choose_predicted_pair(
        predicted_structure[0], sequence
    )
    predicted_token, predicted_ca, predicted_heavy = chain_coordinates(
        predicted_chain, sequence
    )

    system_root = ground_truth_root / system_id
    reference_token, reference_heavy = read_reference_protein(
        system_root / "receptor.cif", reference_chain_id, sequence
    )
    reference_structure = gemmi.read_structure(str(system_root / "receptor.cif"))
    _, reference_ca, _ = chain_coordinates(
        reference_chain(reference_structure[0], reference_chain_id), sequence
    )

    with molecule_path.open("rb") as handle:
        predicted_molecule = remove_nesso_tokenization_hydrogens(
            pickle.load(handle)["LIG1"]  # noqa: S301
        )
    predicted_ligand_coordinates = molecule_coordinates(
        predicted_molecule, predicted_ligand
    )
    ligand_sdf = (
        system_root
        / "ligand_files"
        / f"{sample['ligand']['instance_chain']}.sdf"
    )
    reference_molecule, reference_ligand_coordinates = reference_ligand_molecule(
        ligand_sdf
    )
    mappings, mapping_mode = graph_mappings(predicted_molecule, reference_molecule)

    predicted_distances = _protein_ligand_distances(
        predicted_token, predicted_ligand_coordinates
    )
    candidates: list[dict[str, Any]] = []
    for mapping in mappings:
        mapped_reference_ligand = reference_ligand_coordinates[
            np.asarray(mapping, dtype=int)
        ]
        reference_distances = _protein_ligand_distances(
            reference_token, mapped_reference_ligand
        )
        valid = np.isfinite(predicted_distances) & np.isfinite(reference_distances)
        interface = valid & (reference_distances <= interface_cutoff)
        if not np.any(interface):
            continue
        mae = float(
            np.mean(np.abs(predicted_distances[interface] - reference_distances[interface]))
        )
        true_contacts = reference_distances[valid] <= contact_cutoff
        predicted_contacts = predicted_distances[valid] <= contact_cutoff
        contact_precision, contact_recall, contact_f1 = precision_recall_f1(
            true_contacts, predicted_contacts
        )
        candidates.append(
            {
                "mapping": tuple(mapping),
                "reference_distances": reference_distances,
                "reference_ligand": mapped_reference_ligand,
                "mae": mae,
                "contact_precision": contact_precision,
                "contact_recall": contact_recall,
                "contact_f1": contact_f1,
            }
        )
    if not candidates:
        raise ValueError("no resolved reference interface pairs for any ligand mapping")
    mae_choice = min(candidates, key=lambda item: (item["mae"], item["mapping"]))
    contact_choice = min(
        candidates, key=lambda item: (-item["contact_f1"], item["mapping"])
    )
    interface_mae = mae_choice["mae"]
    mapped_reference_ligand = mae_choice["reference_ligand"]

    true_pocket, reference_resolved = _physical_pocket_labels(
        reference_heavy, mapped_reference_ligand, contact_cutoff
    )
    predicted_pocket, predicted_resolved = _physical_pocket_labels(
        predicted_heavy, predicted_ligand_coordinates, contact_cutoff
    )
    comparable_residues = reference_resolved & predicted_resolved
    pocket_precision, pocket_recall, pocket_f1 = precision_recall_f1(
        true_pocket[comparable_residues], predicted_pocket[comparable_residues]
    )

    ca_mask = np.all(np.isfinite(predicted_ca), axis=1) & np.all(
        np.isfinite(reference_ca), axis=1
    )
    if np.sum(ca_mask) < 3:
        raise ValueError("fewer than three mapped C-alpha atoms")
    rotation, translation, protein_ca_rmsd = kabsch(
        predicted_ca[ca_mask], reference_ca[ca_mask]
    )
    ligand_rmsd, ligand_rmsd_mapping = minimum_mapped_rmsd(
        predicted_ligand_coordinates,
        reference_ligand_coordinates,
        mappings,
        rotation=rotation,
        translation=translation,
    )

    familiarity_score, familiarity_bin = familiarity_metadata(sample, "primary")
    confidence = (
        json.loads(confidence_path.read_text()) if confidence_path.is_file() else {}
    )
    return {
        "sample_id": sample_id,
        "system_id": system_id,
        "runs_n_poses_cluster": sample["runs_n_poses"]["cluster"],
        "familiarity_source": "primary_june2023",
        "familiarity_score_0_to_100": familiarity_score,
        "familiarity_bin": familiarity_bin,
        "prediction_file": str(prediction_cif),
        "reference_system_directory": str(system_root),
        "protein_tokens": len(sequence),
        "resolved_reference_token_atoms": int(
            np.sum(np.all(np.isfinite(reference_token), axis=1))
        ),
        "resolved_predicted_token_atoms": int(
            np.sum(np.all(np.isfinite(predicted_token), axis=1))
        ),
        "ligand_heavy_atoms": predicted_molecule.GetNumAtoms(),
        "predicted_selection": predicted_selection,
        "ligand_mapping": {
            "mode": mapping_mode,
            "symmetry_mappings_considered": len(mappings),
            "selection_rule": "metric-specific optimum over exact graph symmetries; ties lexicographic",
            "interface_mae_mapping": list(mae_choice["mapping"]),
            "contact_f1_mapping": list(contact_choice["mapping"]),
            "protein_aligned_rmsd_mapping": list(ligand_rmsd_mapping),
        },
        "metrics": {
            "interface_expected_distance_mae_angstrom": interface_mae,
            "token_contact_precision_at_6a": contact_choice["contact_precision"],
            "token_contact_recall_at_6a": contact_choice["contact_recall"],
            "token_contact_f1_at_6a": contact_choice["contact_f1"],
            "physical_pocket_precision_at_6a": pocket_precision,
            "physical_pocket_recall_at_6a": pocket_recall,
            "physical_pocket_f1_at_6a": pocket_f1,
            "protein_global_ca_rmsd_angstrom": protein_ca_rmsd,
            "ligand_heavy_atom_rmsd_after_global_protein_alignment_angstrom": ligand_rmsd,
        },
        "boltz_confidence": confidence,
    }


def aggregate_results(
    records: list[dict[str, Any]],
    failures: list[dict[str, Any]],
    *,
    iterations: int = 2000,
    seed: int = 2026082802,
) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "systems_requested": len(records) + len(failures),
        "systems_scored": len(records),
        "systems_failed": len(failures),
        "overall": {},
        "by_familiarity_bin": {},
        "spearman_vs_familiarity": {},
        "bootstrap": {
            "unit": "Runs N' Poses structural cluster",
            "iterations": iterations,
            "seed": seed,
        },
    }
    if not records:
        return summary
    rng = np.random.default_rng(seed)
    familiarity = np.asarray(
        [record["familiarity_score_0_to_100"] for record in records], dtype=float
    )
    clusters = np.asarray(
        [record["runs_n_poses_cluster"] for record in records], dtype=object
    )
    metrics = COMMON_METRICS + (
        "protein_global_ca_rmsd_angstrom",
        "ligand_heavy_atom_rmsd_after_global_protein_alignment_angstrom",
    )
    for metric in metrics:
        values = np.asarray([record["metrics"][metric] for record in records], dtype=float)
        finite = np.isfinite(values) & np.isfinite(familiarity)
        summary["overall"][metric] = {
            "n": int(np.sum(finite)),
            "mean": float(np.mean(values[finite])),
            "median": float(np.median(values[finite])),
        }
        paired = np.column_stack((familiarity[finite], values[finite]))

        def statistic(sampled: np.ndarray) -> float:
            return _spearman(sampled[:, 0], sampled[:, 1])

        interval = _cluster_bootstrap_interval(
            paired,
            clusters[finite],
            statistic,
            iterations=iterations,
            rng=rng,
        )
        rho = _spearman(familiarity[finite], values[finite])
        summary["spearman_vs_familiarity"][metric] = {
            "n": int(np.sum(finite)),
            "rho": rho if np.isfinite(rho) else None,
            **interval,
        }
    for bin_name in sorted(
        {record["familiarity_bin"] for record in records},
        key=lambda name: float(name.split("-")[0]),
    ):
        chosen = [record for record in records if record["familiarity_bin"] == bin_name]
        summary["by_familiarity_bin"][bin_name] = {
            "n": len(chosen),
            "metrics": {
                metric: {
                    "median": float(
                        np.median([record["metrics"][metric] for record in chosen])
                    )
                }
                for metric in metrics
            },
        }
    return summary


def write_csv(records: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "sample_id",
        "system_id",
        "runs_n_poses_cluster",
        "familiarity_score_0_to_100",
        "familiarity_bin",
        "protein_tokens",
        "ligand_heavy_atoms",
        *COMMON_METRICS,
        "protein_global_ca_rmsd_angstrom",
        "ligand_heavy_atom_rmsd_after_global_protein_alignment_angstrom",
        "confidence_score",
        "ligand_iptm",
        "complex_iplddt",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            row = {key: record.get(key) for key in fields}
            row.update(record["metrics"])
            for key in ("confidence_score", "ligand_iptm", "complex_iplddt"):
                row[key] = record["boltz_confidence"].get(key)
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ground-truth-root", type=Path, required=True)
    parser.add_argument("--prediction-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample-id", action="append", default=[])
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    requested = set(args.sample_id)
    samples = [
        sample
        for sample in manifest["samples"]
        if not requested or sample["sample_id"] in requested
    ]
    missing = requested - {sample["sample_id"] for sample in samples}
    if missing:
        raise ValueError(f"sample IDs absent from manifest: {sorted(missing)}")

    records: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for sample in samples:
        try:
            records.append(
                score_sample(sample, args.prediction_root, args.ground_truth_root)
            )
        except Exception as error:
            failures.append(
                {
                    "sample_id": sample["sample_id"],
                    "stage": "scoring",
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "schema_version": 1,
        "experiment_id": "exp015_boltz2_nesso1_rnp_postcutoff50",
        "reference_coordinates_used_during_inference": False,
        "metric_definitions": {
            "token_distance": "protein C-beta (C-alpha for glycine) to ligand heavy atom",
            "interface_mae": "mean absolute token-distance error for resolved reference pairs within 15 A",
            "contacts": "binary token distance <=6 A",
            "physical_pocket": "protein residue with any heavy atom <=6 A from ligand",
            "protein_rmsd": "C-alpha RMSD after sequence-indexed global Kabsch fit",
            "ligand_rmsd": "symmetry-corrected heavy-atom RMSD after the fixed global protein transform; ligand atoms are not fitted",
        },
        "records": records,
        "failures": failures,
        "aggregate": aggregate_results(records, failures),
    }
    (args.output_dir / "per_system_metrics.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    write_csv(records, args.output_dir / "per_system_metrics.csv")
    print(f"Scored {len(records)} systems with {len(failures)} failures")


if __name__ == "__main__":
    main()
