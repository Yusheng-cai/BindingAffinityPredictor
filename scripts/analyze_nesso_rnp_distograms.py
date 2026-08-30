#!/usr/bin/env python3
"""Score Nesso-1 distograms against selected Runs N' Poses complexes.

Run this script in the pinned Nesso environment because it uses gemmi, RDKit,
safetensors, and torch. Experimental coordinates are used only here, after
inference; they are never passed to Nesso.
"""

from __future__ import annotations

import argparse
import csv
import json
import warnings
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
from rdkit import Chem
from safetensors.torch import load_file
from scipy.stats import rankdata, spearmanr

from affinity_benchmark.metrics.distogram import (
    average_precision,
    binary_auroc,
    contact_probabilities,
    distogram_negative_log_likelihood,
    expected_distances,
    precision_recall_f1,
)


AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

FAMILIARITY_BINS = (
    "0-20",
    "20-30",
    "30-40",
    "40-50",
    "50-60",
    "60-70",
    "70-80",
    "80-100",
)

STRUCTURAL_METRICS = (
    "interface_distogram_nll",
    "interface_expected_distance_mae_angstrom",
    "token_contact_average_precision_at_6a",
    "token_contact_auroc_at_6a",
    "token_contact_precision_at_6a",
    "token_contact_recall_at_6a",
    "token_contact_f1_at_6a",
    "physical_pocket_precision_at_6a",
    "physical_pocket_recall_at_6a",
    "physical_pocket_f1_at_6a",
)


def familiarity_metadata(
    sample: dict[str, Any], source: str = "primary"
) -> tuple[float, str]:
    """Resolve model-appropriate Runs N' Poses familiarity metadata."""

    rnp = sample["runs_n_poses"]
    if source == "primary":
        return float(rnp["similarity_score_0_to_100"]), rnp["similarity_bin"]
    if source != "nesso1_sep2021":
        raise ValueError(f"unknown familiarity source: {source}")
    score = float(rnp["nesso1_sep2021_familiarity"]["similarity_score_0_to_100"])
    for bin_name in FAMILIARITY_BINS:
        lower, upper = (float(value) for value in bin_name.split("-"))
        if lower <= score < upper or (score == 100.0 and upper == 100.0):
            return score, bin_name
    raise ValueError(f"familiarity score outside [0, 100]: {score}")


def _spearman(first: np.ndarray, second: np.ndarray) -> float:
    if len(first) < 2:
        return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spearmanr(first, second).statistic)


def _partial_spearman(
    first: np.ndarray, second: np.ndarray, controls: list[np.ndarray]
) -> float:
    """Spearman correlation after linear residualization of ranked controls."""

    if len(first) < 3:
        return float("nan")
    ranked_first = rankdata(first)
    ranked_second = rankdata(second)
    design = np.column_stack(
        [np.ones(len(first), dtype=float)] + [rankdata(control) for control in controls]
    )
    residual_first = ranked_first - design @ np.linalg.lstsq(
        design, ranked_first, rcond=None
    )[0]
    residual_second = ranked_second - design @ np.linalg.lstsq(
        design, ranked_second, rcond=None
    )[0]
    if np.std(residual_first) == 0 or np.std(residual_second) == 0:
        return float("nan")
    return float(np.corrcoef(residual_first, residual_second)[0, 1])


def _float_or_default(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_reference_protein(
    cif_path: Path, chain_id: str, sequence: str
) -> tuple[np.ndarray, list[np.ndarray | None]]:
    """Read C-beta/C-alpha token coordinates and per-residue heavy atoms."""

    block = gemmi.cif.read(str(cif_path)).sole_block()
    names = (
        "label_asym_id",
        "label_seq_id",
        "label_comp_id",
        "label_atom_id",
        "type_symbol",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "label_alt_id",
        "pdbx_PDB_model_num",
    )
    columns = {name: block.find_values(f"_atom_site.{name}") for name in names}
    if not columns["label_asym_id"]:
        raise ValueError(f"no atom_site records in {cif_path}")

    # Keep the highest-occupancy alternate for each named atom.
    atoms: dict[int, dict[str, tuple[float, str, np.ndarray, str]]] = {}
    for row in range(len(columns["label_asym_id"])):
        if columns["label_asym_id"][row] != chain_id:
            continue
        model = columns["pdbx_PDB_model_num"][row]
        if model not in (".", "?", "1"):
            continue
        raw_seq_id = columns["label_seq_id"][row]
        if raw_seq_id in (".", "?"):
            continue
        seq_index = int(raw_seq_id) - 1
        if seq_index < 0 or seq_index >= len(sequence):
            raise ValueError(f"sequence index {seq_index + 1} is outside {cif_path}")
        atom_name = columns["label_atom_id"][row]
        element = columns["type_symbol"][row].upper()
        if element == "H":
            continue
        coordinate = np.array(
            [
                float(columns["Cartn_x"][row]),
                float(columns["Cartn_y"][row]),
                float(columns["Cartn_z"][row]),
            ],
            dtype=float,
        )
        occupancy = _float_or_default(columns["occupancy"][row], 1.0)
        alt_id = columns["label_alt_id"][row]
        residue_name = columns["label_comp_id"][row]
        previous = atoms.setdefault(seq_index, {}).get(atom_name)
        rank = (occupancy, alt_id in (".", "?", "A"))
        if previous is None or rank > (previous[0], previous[1] in (".", "?", "A")):
            atoms[seq_index][atom_name] = (
                occupancy,
                alt_id,
                coordinate,
                residue_name,
            )

    token_coordinates = np.full((len(sequence), 3), np.nan, dtype=float)
    heavy_coordinates: list[np.ndarray | None] = [None] * len(sequence)
    for seq_index, residue_atoms in atoms.items():
        residue_names = {record[3] for record in residue_atoms.values()}
        if len(residue_names) != 1:
            raise ValueError(f"mixed residue identities at position {seq_index + 1}")
        residue_name = next(iter(residue_names))
        observed_aa = AA3_TO_AA1.get(residue_name)
        if observed_aa is not None and observed_aa != sequence[seq_index]:
            raise ValueError(
                f"sequence mismatch at {chain_id}:{seq_index + 1}: "
                f"manifest={sequence[seq_index]} CIF={observed_aa}"
            )
        distogram_atom = "CA" if residue_name == "GLY" else "CB"
        if distogram_atom in residue_atoms:
            token_coordinates[seq_index] = residue_atoms[distogram_atom][2]
        heavy_coordinates[seq_index] = np.stack(
            [record[2] for record in residue_atoms.values()]
        )

    return token_coordinates, heavy_coordinates


def remove_nesso_tokenization_hydrogens(molecule: Chem.Mol) -> Chem.Mol:
    """Remove all atoms Nesso skips while retaining heavy-atom token order."""

    remove_hydrogens = Chem.RemoveHsParameters()
    remove_hydrogens.removeDefiningBondStereo = True
    heavy_molecule = Chem.RemoveHs(molecule, remove_hydrogens, sanitize=False)
    Chem.SanitizeMol(heavy_molecule)
    return heavy_molecule


def read_ligand_and_mappings(
    sdf_path: Path, input_smiles: str
) -> tuple[np.ndarray, list[tuple[int, ...]], str]:
    """Return reference coordinates and input-token to SDF atom mappings."""

    input_mol = Chem.MolFromSmiles(input_smiles)
    reference_mol = Chem.MolFromMolFile(str(sdf_path), removeHs=True, sanitize=True)
    if input_mol is None or reference_mol is None:
        raise ValueError(f"RDKit could not parse input SMILES or {sdf_path}")
    # Nesso's nonpolymer tokenizer skips every hydrogen atom, including an
    # explicit [H] used only to define double-bond stereochemistry in a SMILES.
    # RDKit retains those stereo-defining hydrogens by default, so request their
    # removal explicitly while preserving the order of the remaining atoms.
    input_mol = remove_nesso_tokenization_hydrogens(input_mol)
    reference_mol = Chem.RemoveHs(reference_mol)
    if input_mol.GetNumAtoms() != reference_mol.GetNumAtoms():
        raise ValueError(
            f"ligand heavy-atom mismatch: input={input_mol.GetNumAtoms()} "
            f"reference={reference_mol.GetNumAtoms()}"
        )

    mapping_mode = "exact_graph_with_chirality"
    mappings = list(
        reference_mol.GetSubstructMatches(
            input_mol, uniquify=False, useChirality=True, maxMatches=10000
        )
    )
    if not mappings:
        mapping_mode = "exact_graph_without_chirality"
        mappings = list(
            reference_mol.GetSubstructMatches(
                input_mol, uniquify=False, useChirality=False, maxMatches=10000
            )
        )
    mappings = sorted(set(tuple(int(index) for index in mapping) for mapping in mappings))
    if not mappings or any(len(mapping) != input_mol.GetNumAtoms() for mapping in mappings):
        raise ValueError(f"no exact ligand graph isomorphism for {sdf_path}")

    conformer = reference_mol.GetConformer()
    coordinates = np.array(
        [list(conformer.GetAtomPosition(index)) for index in range(reference_mol.GetNumAtoms())],
        dtype=float,
    )
    return coordinates, mappings, mapping_mode


def _protein_ligand_distances(
    protein_coordinates: np.ndarray, ligand_coordinates: np.ndarray
) -> np.ndarray:
    return np.linalg.norm(
        protein_coordinates[:, None, :] - ligand_coordinates[None, :, :], axis=-1
    )


def _physical_pocket_labels(
    residue_heavy_coordinates: list[np.ndarray | None],
    ligand_coordinates: np.ndarray,
    cutoff: float,
) -> tuple[np.ndarray, np.ndarray]:
    labels = np.zeros(len(residue_heavy_coordinates), dtype=bool)
    resolved = np.zeros(len(residue_heavy_coordinates), dtype=bool)
    for index, coordinates in enumerate(residue_heavy_coordinates):
        if coordinates is None or len(coordinates) == 0:
            continue
        resolved[index] = True
        minimum = np.min(
            np.linalg.norm(coordinates[:, None, :] - ligand_coordinates[None, :, :], axis=-1)
        )
        labels[index] = minimum <= cutoff
    return labels, resolved


def score_sample(
    sample: dict[str, Any],
    tensor_path: Path,
    ground_truth_root: Path,
    interface_cutoff: float = 15.0,
    contact_cutoff: float = 6.0,
    familiarity_source: str = "primary",
    symmetry_policy: str = "legacy_nll",
) -> dict[str, Any]:
    """Score one Nesso output, choosing only among exact ligand symmetries."""

    system_id = sample["runs_n_poses"]["system_id"]
    chain = sample["protein"]["chains"][0]
    sequence = chain["sequence"]
    system_root = ground_truth_root / system_id
    protein_coordinates, residue_heavy_coordinates = read_reference_protein(
        system_root / "receptor.cif", chain["id"], sequence
    )
    ligand_coordinates_sdf, mappings, mapping_mode = read_ligand_and_mappings(
        system_root / "ligand_files" / f"{sample['ligand']['instance_chain']}.sdf",
        sample["ligand"]["input_smiles"],
    )

    tensors = load_file(str(tensor_path), device="cpu")
    logits = tensors["pdistogram"].float().numpy()
    refine_mask = tensors["refine_mask"].bool().numpy()
    pocket_mask = tensors["pocket_mask"].bool().numpy()
    token_pad_mask = tensors["token_pad_mask"].bool().numpy()
    protein_tokens = len(sequence)
    ligand_tokens = len(mappings[0])
    expected_tokens = protein_tokens + ligand_tokens
    if logits.shape != (expected_tokens, expected_tokens, 64):
        raise ValueError(
            f"unexpected pdistogram shape {logits.shape}; expected "
            f"({expected_tokens}, {expected_tokens}, 64)"
        )
    for name, array in (
        ("refine_mask", refine_mask),
        ("pocket_mask", pocket_mask),
        ("token_pad_mask", token_pad_mask),
    ):
        if array.shape != (expected_tokens,):
            raise ValueError(f"unexpected {name} shape {array.shape}")
    if not np.all(token_pad_mask):
        raise ValueError("unpadded single-record inference should have all tokens valid")

    pl_logits = logits[:protein_tokens, protein_tokens:, :]
    predicted_distances = expected_distances(pl_logits)
    predicted_contacts = contact_probabilities(pl_logits, cutoff=contact_cutoff)

    mapping_scores: list[dict[str, Any]] = []
    for mapping in mappings:
        ligand_in_input_order = ligand_coordinates_sdf[np.asarray(mapping, dtype=int)]
        true_distances = _protein_ligand_distances(
            protein_coordinates, ligand_in_input_order
        )
        valid_pairs = np.isfinite(true_distances)
        interface_pairs = valid_pairs & (true_distances <= interface_cutoff)
        nll = distogram_negative_log_likelihood(
            pl_logits, true_distances, mask=interface_pairs
        )
        expected_distance_mae = float(
            np.mean(
                np.abs(
                    predicted_distances[interface_pairs]
                    - true_distances[interface_pairs]
                )
            )
        )
        contact_labels = true_distances[valid_pairs] <= contact_cutoff
        contact_scores = predicted_contacts[valid_pairs]
        contact_predictions = predicted_distances[valid_pairs] <= contact_cutoff
        contact_precision, contact_recall, contact_f1 = precision_recall_f1(
            contact_labels, contact_predictions
        )
        mapping_scores.append(
            {
                "mapping": tuple(mapping),
                "true_distances": true_distances,
                "ligand_coordinates": ligand_in_input_order,
                "nll": nll,
                "mae": expected_distance_mae,
                "contact_average_precision": average_precision(
                    contact_labels, contact_scores
                ),
                "contact_auroc": binary_auroc(contact_labels, contact_scores),
                "contact_precision": contact_precision,
                "contact_recall": contact_recall,
                "contact_f1": contact_f1,
            }
        )
    nll_choice = min(mapping_scores, key=lambda item: (item["nll"], item["mapping"]))
    if symmetry_policy == "legacy_nll":
        # Preserve the accepted exp012-exp014 protocol and its published result.
        mae_choice = ap_choice = auroc_choice = contact_choice = nll_choice
    elif symmetry_policy == "metric_specific":
        mae_choice = min(
            mapping_scores, key=lambda item: (item["mae"], item["mapping"])
        )
        ap_choice = min(
            mapping_scores,
            key=lambda item: (-item["contact_average_precision"], item["mapping"]),
        )
        auroc_choice = min(
            mapping_scores,
            key=lambda item: (-item["contact_auroc"], item["mapping"]),
        )
        contact_choice = min(
            mapping_scores, key=lambda item: (-item["contact_f1"], item["mapping"])
        )
    else:
        raise ValueError(f"unknown ligand symmetry policy: {symmetry_policy}")
    true_distances = mae_choice["true_distances"]
    ligand_in_input_order = mae_choice["ligand_coordinates"]
    valid_pairs = np.isfinite(true_distances)
    interface_pairs = valid_pairs & (true_distances <= interface_cutoff)
    contact_labels = true_distances[valid_pairs] <= contact_cutoff

    true_pocket, resolved_residues = _physical_pocket_labels(
        residue_heavy_coordinates, ligand_in_input_order, contact_cutoff
    )
    # Match Nesso's own pocket construction (minimum expected token distance),
    # but apply the benchmark's 6 A cutoff. The saved native pocket_mask uses
    # the affinity-head cutoff (15 A) and therefore is not a 6 A prediction.
    predicted_pocket = np.min(predicted_distances, axis=1) <= contact_cutoff
    pocket_precision, pocket_recall, pocket_f1 = precision_recall_f1(
        true_pocket[resolved_residues], predicted_pocket[resolved_residues]
    )

    affinity_path = tensor_path.with_name("affinity.json")
    native_outputs = json.loads(affinity_path.read_text()) if affinity_path.exists() else {}
    rnp = sample["runs_n_poses"]
    familiarity_score, familiarity_bin = familiarity_metadata(
        sample, familiarity_source
    )
    return {
        "sample_id": sample["sample_id"],
        "system_id": system_id,
        "runs_n_poses_cluster": rnp["cluster"],
        "familiarity_source": familiarity_source,
        "familiarity_score_0_to_100": familiarity_score,
        "familiarity_bin": familiarity_bin,
        "prediction_file": str(tensor_path),
        "reference_system_directory": str(system_root),
        "protein_tokens": protein_tokens,
        "resolved_protein_distogram_atoms": int(
            np.sum(np.all(np.isfinite(protein_coordinates), axis=1))
        ),
        "ligand_tokens": ligand_tokens,
        "refined_protein_tokens": int(np.sum(refine_mask[:protein_tokens])),
        "predicted_token_pocket_residues_at_6a": int(np.sum(predicted_pocket)),
        "native_affinity_pocket_mask_protein_residues_at_15a": int(
            np.sum(pocket_mask[:protein_tokens])
        ),
        "reference_pocket_residues_6a_any_heavy_atom": int(np.sum(true_pocket)),
        "reference_interface_pairs_within_15a": int(np.sum(interface_pairs)),
        "reference_token_contacts_within_6a": int(np.sum(contact_labels)),
        "ligand_mapping": {
            "mode": mapping_mode,
            "symmetry_mappings_considered": len(mappings),
            "symmetry_policy": symmetry_policy,
            "selection_rule": (
                "metric-specific optimum over exact graph symmetries; ties lexicographic"
                if symmetry_policy == "metric_specific"
                else "minimum interface distogram NLL for all metrics; ties lexicographic"
            ),
            "distogram_nll_mapping": list(nll_choice["mapping"]),
            "interface_mae_mapping": list(mae_choice["mapping"]),
            "contact_ap_mapping": list(ap_choice["mapping"]),
            "contact_auroc_mapping": list(auroc_choice["mapping"]),
            "contact_f1_mapping": list(contact_choice["mapping"]),
        },
        "metrics": {
            "interface_distogram_nll": nll_choice["nll"],
            "interface_expected_distance_mae_angstrom": mae_choice["mae"],
            "token_contact_average_precision_at_6a": ap_choice[
                "contact_average_precision"
            ],
            "token_contact_auroc_at_6a": auroc_choice["contact_auroc"],
            "token_contact_precision_at_6a": contact_choice["contact_precision"],
            "token_contact_recall_at_6a": contact_choice["contact_recall"],
            "token_contact_f1_at_6a": contact_choice["contact_f1"],
            "physical_pocket_precision_at_6a": pocket_precision,
            "physical_pocket_recall_at_6a": pocket_recall,
            "physical_pocket_f1_at_6a": pocket_f1,
        },
        "native_nesso_outputs": native_outputs,
    }


def _write_csv(records: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "sample_id",
        "system_id",
        "runs_n_poses_cluster",
        "familiarity_score_0_to_100",
        "familiarity_bin",
        "protein_tokens",
        "resolved_protein_distogram_atoms",
        "ligand_tokens",
        "refined_protein_tokens",
        "predicted_token_pocket_residues_at_6a",
        "native_affinity_pocket_mask_protein_residues_at_15a",
        "reference_pocket_residues_6a_any_heavy_atom",
        "reference_interface_pairs_within_15a",
        "reference_token_contacts_within_6a",
        "interface_distogram_nll",
        "interface_expected_distance_mae_angstrom",
        "token_contact_average_precision_at_6a",
        "token_contact_auroc_at_6a",
        "token_contact_precision_at_6a",
        "token_contact_recall_at_6a",
        "token_contact_f1_at_6a",
        "physical_pocket_precision_at_6a",
        "physical_pocket_recall_at_6a",
        "physical_pocket_f1_at_6a",
        "entropy_crop_pl",
        "affinity_pred_value",
        "affinity_pred_value1",
        "affinity_pred_value2",
        "affinity_probability_binary",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            row = {key: record.get(key) for key in fields}
            row.update(record["metrics"])
            for key in fields:
                if key in record["native_nesso_outputs"]:
                    row[key] = record["native_nesso_outputs"][key]
            writer.writerow(row)


def _bootstrap_interval(
    values: np.ndarray,
    statistic,
    *,
    iterations: int,
    rng: np.random.Generator,
) -> dict[str, float | int | None]:
    estimates = []
    for _ in range(iterations):
        indices = rng.integers(0, len(values), size=len(values))
        estimate = float(statistic(values[indices]))
        if np.isfinite(estimate):
            estimates.append(estimate)
    array = np.asarray(estimates, dtype=float)
    if len(array) == 0:
        return {
            "iterations_requested": iterations,
            "iterations_finite": 0,
            "lower_95": None,
            "upper_95": None,
        }
    return {
        "iterations_requested": iterations,
        "iterations_finite": int(len(array)),
        "lower_95": float(np.quantile(array, 0.025)),
        "upper_95": float(np.quantile(array, 0.975)),
    }


def _cluster_bootstrap_interval(
    values: np.ndarray,
    clusters: np.ndarray,
    statistic,
    *,
    iterations: int,
    rng: np.random.Generator,
) -> dict[str, float | int | None]:
    """Resample complete Runs N' Poses clusters with replacement."""

    unique_clusters = np.unique(clusters)
    estimates = []
    for _ in range(iterations):
        sampled_clusters = rng.choice(
            unique_clusters, size=len(unique_clusters), replace=True
        )
        indices = np.concatenate(
            [np.flatnonzero(clusters == cluster) for cluster in sampled_clusters]
        )
        estimate = float(statistic(values[indices]))
        if np.isfinite(estimate):
            estimates.append(estimate)
    array = np.asarray(estimates, dtype=float)
    if len(array) == 0:
        return {
            "iterations_requested": iterations,
            "iterations_finite": 0,
            "lower_95": None,
            "upper_95": None,
        }
    return {
        "iterations_requested": iterations,
        "iterations_finite": int(len(array)),
        "lower_95": float(np.quantile(array, 0.025)),
        "upper_95": float(np.quantile(array, 0.975)),
    }


def aggregate_results(
    records: list[dict[str, Any]],
    failures: list[dict[str, Any]],
    *,
    iterations: int = 2000,
    seed: int = 20260820,
) -> dict[str, Any]:
    """Summarize metrics by familiarity and bootstrap systems/clusters."""

    rng = np.random.default_rng(seed)
    summary: dict[str, Any] = {
        "systems_requested": len(records) + len(failures),
        "systems_scored": len(records),
        "systems_failed": len(failures),
        "bootstrap": {
            "unit": "Runs N' Poses cluster (all records in a sampled cluster retained)",
            "iterations": iterations,
            "seed": seed,
            "interval": "percentile 95%",
        },
        "overall": {},
        "by_familiarity_bin": {},
        "spearman_vs_familiarity": {},
        "entropy_error_spearman": {},
        "exploratory_entropy_error_partial_spearman": {},
    }
    if not records:
        return summary

    familiarity = np.asarray(
        [record["familiarity_score_0_to_100"] for record in records], dtype=float
    )
    clusters = np.asarray(
        [record["runs_n_poses_cluster"] for record in records], dtype=object
    )
    for metric in STRUCTURAL_METRICS:
        values = np.asarray([record["metrics"][metric] for record in records], dtype=float)
        finite = np.isfinite(values) & np.isfinite(familiarity)
        values_finite = values[finite]
        familiarity_finite = familiarity[finite]
        summary["overall"][metric] = {
            "n": int(len(values_finite)),
            "mean": float(np.mean(values_finite)),
            "median": float(np.median(values_finite)),
        }
        observed_rho = _spearman(familiarity_finite, values_finite)

        paired = np.column_stack((familiarity_finite, values_finite))

        def bootstrap_spearman(resampled: np.ndarray) -> float:
            return _spearman(resampled[:, 0], resampled[:, 1])

        interval = _cluster_bootstrap_interval(
            paired,
            clusters[finite],
            bootstrap_spearman,
            iterations=iterations,
            rng=rng,
        )
        summary["spearman_vs_familiarity"][metric] = {
            "n": int(len(values_finite)),
            "rho": observed_rho if np.isfinite(observed_rho) else None,
            **interval,
        }

    for bin_name in FAMILIARITY_BINS:
        selected = [record for record in records if record["familiarity_bin"] == bin_name]
        selected_clusters = np.asarray(
            [record["runs_n_poses_cluster"] for record in selected], dtype=object
        )
        bin_result: dict[str, Any] = {"n": len(selected), "metrics": {}}
        for metric in STRUCTURAL_METRICS:
            values = np.asarray(
                [record["metrics"][metric] for record in selected], dtype=float
            )
            finite = np.isfinite(values)
            values = values[finite]
            metric_clusters = selected_clusters[finite]
            if len(values) == 0:
                continue
            interval = _cluster_bootstrap_interval(
                values,
                metric_clusters,
                np.median,
                iterations=iterations,
                rng=rng,
            )
            bin_result["metrics"][metric] = {
                "n": int(len(values)),
                "median": float(np.median(values)),
                **interval,
            }
        summary["by_familiarity_bin"][bin_name] = bin_result

    entropy = np.asarray(
        [record["native_nesso_outputs"].get("entropy_crop_pl", np.nan) for record in records],
        dtype=float,
    )
    protein_tokens = np.asarray([record["protein_tokens"] for record in records], dtype=float)
    ligand_tokens = np.asarray([record["ligand_tokens"] for record in records], dtype=float)
    for error_metric in (
        "interface_distogram_nll",
        "interface_expected_distance_mae_angstrom",
    ):
        errors = np.asarray(
            [record["metrics"][error_metric] for record in records], dtype=float
        )
        finite = np.isfinite(entropy) & np.isfinite(errors)
        paired = np.column_stack((entropy[finite], errors[finite]))

        def bootstrap_entropy_error(resampled: np.ndarray) -> float:
            return _spearman(resampled[:, 0], resampled[:, 1])

        interval = _cluster_bootstrap_interval(
            paired,
            clusters[finite],
            bootstrap_entropy_error,
            iterations=iterations,
            rng=rng,
        )
        summary["entropy_error_spearman"][error_metric] = {
            "n": int(np.sum(finite)),
            "rho": (
                value
                if np.isfinite(value := _spearman(entropy[finite], errors[finite]))
                else None
            ),
            **interval,
        }
        control_sets = {
            "controlling_familiarity": [familiarity],
            "controlling_familiarity_protein_and_ligand_size": [
                familiarity,
                protein_tokens,
                ligand_tokens,
            ],
        }
        partial_results = {}
        for control_name, controls in control_sets.items():
            finite_partial = np.isfinite(entropy) & np.isfinite(errors)
            for control in controls:
                finite_partial &= np.isfinite(control)
            columns = [entropy[finite_partial], errors[finite_partial]] + [
                control[finite_partial] for control in controls
            ]
            paired_partial = np.column_stack(columns)

            def bootstrap_partial(resampled: np.ndarray) -> float:
                return _partial_spearman(
                    resampled[:, 0],
                    resampled[:, 1],
                    [resampled[:, index] for index in range(2, resampled.shape[1])],
                )

            partial_interval = _cluster_bootstrap_interval(
                paired_partial,
                clusters[finite_partial],
                bootstrap_partial,
                iterations=iterations,
                rng=rng,
            )
            partial_value = _partial_spearman(
                columns[0], columns[1], columns[2:]
            )
            partial_results[control_name] = {
                "n": int(np.sum(finite_partial)),
                "rho": partial_value if np.isfinite(partial_value) else None,
                **partial_interval,
            }
        summary["exploratory_entropy_error_partial_spearman"][error_metric] = (
            partial_results
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ground-truth-root", type=Path, required=True)
    parser.add_argument("--predictions-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--experiment-id",
        default="exp012_nesso1_rnp_distogram_generalization",
    )
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument(
        "--familiarity-source",
        choices=("primary", "nesso1_sep2021"),
        default="primary",
    )
    parser.add_argument(
        "--symmetry-policy",
        choices=("legacy_nll", "metric_specific"),
        default="legacy_nll",
        help=(
            "legacy_nll preserves exp012-exp014; metric_specific applies the "
            "same per-metric symmetry correction used in the paired exp015 analysis"
        ),
    )
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    requested = set(args.sample_id)
    samples = {
        sample["sample_id"]: sample
        for sample in manifest["samples"]
        if not requested or sample["sample_id"] in requested
    }
    if requested - set(samples):
        raise ValueError(f"sample IDs absent from manifest: {sorted(requested - set(samples))}")

    prediction_files = {
        path.parent.name: path
        for path in args.predictions_root.rglob("predictions.safetensors")
    }
    records = []
    failures = []
    for sample_id in sorted(samples):
        if sample_id not in prediction_files:
            failures.append(
                {
                    "sample_id": sample_id,
                    "stage": "prediction_discovery",
                    "error": "predictions.safetensors not found",
                }
            )
            continue
        try:
            records.append(
                score_sample(
                    samples[sample_id],
                    prediction_files[sample_id],
                    args.ground_truth_root,
                    familiarity_source=args.familiarity_source,
                    symmetry_policy=args.symmetry_policy,
                )
            )
        except Exception as error:  # keep an auditable per-system failure record
            failures.append(
                {
                    "sample_id": sample_id,
                    "stage": "scoring",
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "schema_version": 1,
        "experiment_id": args.experiment_id,
        "manifest": str(args.manifest),
        "metric_definitions": {
            "interface_distogram_nll": (
                "mean negative log probability of the observed Nesso distance bin "
                "for resolved protein-token/ligand-atom pairs at true distance <=15 A"
            ),
            "token_contacts": (
                "C-beta (C-alpha for glycine) to ligand-heavy-atom distance <=6 A"
            ),
            "physical_pocket": (
                "protein residue with any resolved protein heavy atom <=6 A from ligand"
            ),
            "ligand_symmetry": (
                "minimum NLL over exact graph-isomorphic ligand atom mappings"
            ),
        },
        "records": records,
        "failures": failures,
        "aggregate": aggregate_results(records, failures),
    }
    (args.output_dir / "per_system_metrics.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    _write_csv(records, args.output_dir / "per_system_metrics.csv")
    print(
        f"Scored {len(records)} systems with {len(failures)} failures into "
        f"{args.output_dir}"
    )


if __name__ == "__main__":
    main()
