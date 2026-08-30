#!/usr/bin/env python3
"""Select a blinded, familiarity-balanced Runs N' Poses pilot for Nesso-1."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from affinity_benchmark.data.manifest import validate_manifest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ANNOTATIONS = REPOSITORY_ROOT / "data/raw/runs_n_poses_metadata/annotations.csv"
DEFAULT_INPUTS = REPOSITORY_ROOT / "data/raw/runs_n_poses_metadata/inputs.json"
DEFAULT_OUTPUT = REPOSITORY_ROOT / "data/manifests/rnp_nesso1_pilot100.json"

SOURCE_DOI = "10.5281/zenodo.14794785"
SOURCE_RECORD = "18366081"
SOURCE_REPOSITORY_REVISION = "197fafc60a5cef3a9e7f8a4b0dac7c965eed3839"
EXPECTED_SHA256 = {
    "annotations.csv": "259aa0a8f5ea6008d4036886f1ad5ae255689c223eb839fa1a2e7cb54adae609",
    "inputs.json": "85ded25555c1efba0aeb7dc90df2fbf24d9550f390ed3e8ac9b32cf11a0ede4d",
}

# These are the intervals used by the Runs N' Poses plotting code.  The first
# four receive one extra sample because low-familiarity examples are both rarer
# and central to the generalization question.
SIMILARITY_EDGES = (0.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0)
BIN_QUOTAS = (13, 13, 13, 13, 12, 12, 12, 12)
SELECTION_SEED = 20260820
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
PROPERTY_COLUMNS = (
    "protein_total_length",
    "ligand_num_heavy_atoms",
    "ligand_num_rot_bonds",
    "ligand_tpsa",
    "ligand_num_pocket_residues",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def similarity_bin_index(value: float) -> int:
    """Return the Runs N' Poses similarity-bin index for a score in [0, 100]."""

    if not math.isfinite(value) or value < SIMILARITY_EDGES[0] or value > SIMILARITY_EDGES[-1]:
        raise ValueError(f"similarity must be finite and in [0, 100], got {value!r}")
    if value == SIMILARITY_EDGES[-1]:
        return len(SIMILARITY_EDGES) - 2
    return int(np.searchsorted(SIMILARITY_EDGES, value, side="right") - 1)


def similarity_bin_label(index: int) -> str:
    return f"{int(SIMILARITY_EDGES[index])}-{int(SIMILARITY_EDGES[index + 1])}"


def sequence_kmers(sequence: str, k: int = 3) -> frozenset[str]:
    return frozenset(sequence[i : i + k] for i in range(max(0, len(sequence) - k + 1)))


def jaccard(left: frozenset[str], right: frozenset[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def _canonical_smiles(smiles: str) -> tuple[str, Any]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES {smiles!r}")
    return Chem.MolToSmiles(molecule, isomericSmiles=True), MORGAN_GENERATOR.GetFingerprint(molecule)


def prepare_candidates(annotations: pd.DataFrame, inputs: dict[str, Any]) -> pd.DataFrame:
    """Apply the preregistered Nesso-compatible metadata filters."""

    candidates = annotations[
        annotations["ligand_is_proper"]
        & annotations["sucos_shape_pocket_qcov"].notna()
        & annotations["system_id"].isin(inputs)
        & (annotations["num_ligand_chains"] == 1)
        & (annotations["num_proper_ligand_chains"] == 1)
        & (annotations["num_protein_chains"] == 1)
        & annotations["ligand_num_heavy_atoms"].between(6, 50)
    ].copy()

    candidates["input_num_sequences"] = candidates["system_id"].map(
        lambda system_id: len(inputs[system_id]["sequences"])
    )
    candidates["input_num_ligands"] = candidates["system_id"].map(
        lambda system_id: len(inputs[system_id]["smiles"])
    )
    candidates["protein_sequence"] = candidates["system_id"].map(
        lambda system_id: next(iter(inputs[system_id]["sequences"].values()))
        if len(inputs[system_id]["sequences"]) == 1
        else ""
    )
    candidates["protein_chain_id"] = candidates["system_id"].map(
        lambda system_id: next(iter(inputs[system_id]["sequences"]))
        if len(inputs[system_id]["sequences"]) == 1
        else ""
    )
    candidates["protein_total_length"] = candidates["protein_sequence"].str.len()
    candidates["input_ligand_smiles"] = candidates["system_id"].map(
        lambda system_id: inputs[system_id]["smiles"][0]
        if len(inputs[system_id]["smiles"]) == 1
        else ""
    )
    candidates = candidates[
        (candidates["input_num_sequences"] == 1)
        & (candidates["input_num_ligands"] == 1)
        & candidates["protein_total_length"].between(40, 900)
        & candidates["protein_sequence"].map(
            lambda sequence: bool(sequence) and not (set(sequence) - STANDARD_AMINO_ACIDS)
        )
    ].copy()

    canonical_and_fp = candidates["input_ligand_smiles"].map(_canonical_smiles)
    candidates["canonical_smiles"] = canonical_and_fp.map(lambda item: item[0])
    candidates["fingerprint"] = canonical_and_fp.map(lambda item: item[1])
    candidates["protein_kmers"] = candidates["protein_sequence"].map(sequence_kmers)
    candidates["similarity_bin_index"] = candidates["sucos_shape_pocket_qcov"].map(
        similarity_bin_index
    )
    candidates["similarity_bin"] = candidates["similarity_bin_index"].map(
        similarity_bin_label
    )

    properties = candidates.loc[:, PROPERTY_COLUMNS].astype(float)
    scale = properties.std(axis=0, ddof=0).replace(0, 1.0)
    normalized = (properties - properties.mean(axis=0)) / scale
    candidates["property_vector"] = list(normalized.to_numpy())
    return candidates.sort_values("group_key").reset_index(drop=True)


def _combined_distance(left: pd.Series, right: pd.Series) -> float:
    ligand_distance = 1.0 - float(
        DataStructs.TanimotoSimilarity(left["fingerprint"], right["fingerprint"])
    )
    protein_distance = 1.0 - jaccard(left["protein_kmers"], right["protein_kmers"])
    property_distance = min(
        float(np.linalg.norm(left["property_vector"] - right["property_vector"])) / 4.0,
        1.0,
    )
    return 0.50 * ligand_distance + 0.35 * protein_distance + 0.15 * property_distance


def _tie_break(group_key: str, seed: int) -> str:
    return hashlib.sha256(f"{seed}:{group_key}".encode()).hexdigest()


def select_diverse(
    candidates: pd.DataFrame,
    seed: int = SELECTION_SEED,
    bin_quotas: tuple[int, ...] = BIN_QUOTAS,
) -> pd.DataFrame:
    """Round-robin max-min selection with unique cluster, PDB and ligand constraints."""

    if len(bin_quotas) != len(SIMILARITY_EDGES) - 1:
        raise ValueError("bin_quotas must have one entry per familiarity interval")

    selected_indices: list[int] = []
    selected_clusters: set[str] = set()
    selected_pdbs: set[str] = set()
    selected_smiles: set[str] = set()
    counts = [0] * len(bin_quotas)

    while counts != list(bin_quotas):
        made_progress = False
        for bin_index, quota in enumerate(bin_quotas):
            if counts[bin_index] >= quota:
                continue
            eligible = candidates[
                (candidates["similarity_bin_index"] == bin_index)
                & ~candidates.index.isin(selected_indices)
                & ~candidates["cluster"].isin(selected_clusters)
                & ~candidates["entry_pdb_id"].str.upper().isin(selected_pdbs)
                & ~candidates["canonical_smiles"].isin(selected_smiles)
            ]
            if eligible.empty:
                raise RuntimeError(
                    f"cannot fill bin {similarity_bin_label(bin_index)} without relaxing diversity constraints"
                )

            if not selected_indices:
                bin_properties = np.vstack(eligible["property_vector"])
                center = np.median(bin_properties, axis=0)
                scores = {
                    index: -float(np.linalg.norm(row["property_vector"] - center))
                    for index, row in eligible.iterrows()
                }
            else:
                selected_rows = [candidates.loc[index] for index in selected_indices]
                scores = {
                    index: min(_combined_distance(row, chosen) for chosen in selected_rows)
                    for index, row in eligible.iterrows()
                }

            chosen_index = max(
                eligible.index,
                key=lambda index: (scores[index], _tie_break(eligible.loc[index, "group_key"], seed)),
            )
            chosen = candidates.loc[chosen_index]
            selected_indices.append(int(chosen_index))
            selected_clusters.add(str(chosen["cluster"]))
            selected_pdbs.add(str(chosen["entry_pdb_id"]).upper())
            selected_smiles.add(str(chosen["canonical_smiles"]))
            counts[bin_index] += 1
            made_progress = True
        if not made_progress:
            raise RuntimeError("selection made no progress")

    selected = candidates.loc[selected_indices].copy()
    selected["selection_order"] = range(1, len(selected) + 1)
    return selected.sort_values(
        ["similarity_bin_index", "sucos_shape_pocket_qcov", "group_key"]
    ).reset_index(drop=True)


def build_manifest(
    selected: pd.DataFrame,
    source_files: dict[str, str],
    *,
    manifest_id: str = "rnp_nesso1_pilot100_sep2021",
    purpose: str = "Blinded 100-system pilot for Nesso-1 distogram and pocket generalization across Runs N' Poses familiarity bins",
    seed: int = SELECTION_SEED,
    bin_quotas: tuple[int, ...] = BIN_QUOTAS,
    excluded_manifest: str | None = None,
) -> dict[str, Any]:
    samples: list[dict[str, Any]] = []
    for row in selected.itertuples(index=False):
        samples.append(
            {
                "sample_id": row.group_key,
                "target_id": row.system_id,
                "protein": {
                    "name": "Runs N' Poses deposited protein construct",
                    "chains": [
                        {
                            "id": row.protein_chain_id,
                            "sequence": row.protein_sequence,
                            "sequence_length": len(row.protein_sequence),
                        }
                    ],
                },
                "ligand": {
                    "name": row.ligand_ccd_code,
                    "ccd_id": row.ligand_ccd_code,
                    "instance_chain": row.ligand_instance_chain,
                    "input_smiles": row.canonical_smiles,
                    "source_smiles": row.input_ligand_smiles,
                    "heavy_atoms": int(row.ligand_num_heavy_atoms),
                    "molecular_weight_da": float(row.ligand_molecular_weight),
                    "rotatable_bonds": int(row.ligand_num_rot_bonds),
                    "tpsa_angstrom2": float(row.ligand_tpsa),
                },
                "measurement": None,
                "structure_reference": {
                    "pdb_id": str(row.entry_pdb_id).upper(),
                    "use_as_model_input": False,
                    "use_for_scoring_after_prediction": True,
                    "ground_truth_coordinates_downloaded": False,
                },
                "runs_n_poses": {
                    "system_id": row.system_id,
                    "group_key": row.group_key,
                    "cluster": row.cluster,
                    "release_date": row.release_date,
                    "closest_pre_cutoff_system": row.target_system,
                    "closest_system_release_date": row.target_release_date,
                    "similarity_metric": "sucos_shape_pocket_qcov",
                    "similarity_score_0_to_100": float(row.sucos_shape_pocket_qcov),
                    "similarity_bin": row.similarity_bin,
                    "pocket_residues": int(row.ligand_num_pocket_residues),
                    "unique_interactions": int(row.ligand_num_unique_interactions),
                    "selection_order": int(row.selection_order),
                },
                "included": True,
            }
        )

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "manifest_id": manifest_id,
        "created_on": "2026-08-20",
        "purpose": purpose,
        "source_dataset": {
            "name": "Runs N' Poses",
            "doi": SOURCE_DOI,
            "zenodo_record_id": SOURCE_RECORD,
            "repository_url": "https://github.com/plinder-org/runs-n-poses",
            "repository_revision": SOURCE_REPOSITORY_REVISION,
            "source_files": [
                {"name": name, "sha256": digest} for name, digest in sorted(source_files.items())
            ],
        },
        "selection": {
            "status": "frozen_before_nesso_inference",
            "structural_training_cutoff": "2021-09-30",
            "familiarity_metric": "sucos_shape_pocket_qcov",
            "similarity_edges": list(SIMILARITY_EDGES),
            "bin_quotas": {
                similarity_bin_label(index): quota for index, quota in enumerate(bin_quotas)
            },
            "seed": seed,
            "eligibility": {
                "proper_ligand": True,
                "exactly_one_ligand_chain": True,
                "exactly_one_proper_ligand_chain": True,
                "exactly_one_protein_chain": True,
                "protein_length_residues": [40, 900],
                "protein_alphabet": "20_standard_amino_acids",
                "ligand_heavy_atoms": [6, 50],
                "non_null_sep2021_familiarity": True,
                "present_in_inputs_json": True,
            },
            "diversity": {
                "hard_constraints": [
                    "unique_runs_n_poses_cluster",
                    "unique_pdb_entry",
                    "unique_canonical_ligand_smiles",
                ],
                "algorithm": "round_robin_across_bins_then_greedy_max_min",
                "distance_weights": {
                    "ligand_morgan_radius2_2048_tanimoto_distance": 0.50,
                    "protein_3mer_jaccard_distance": 0.35,
                    "standardized_property_distance": 0.15,
                },
                "property_columns": list(PROPERTY_COLUMNS),
            },
            "blinding": "Nesso predictions and experimental coordinates were not used for selection",
        },
        "samples": samples,
    }
    if excluded_manifest is not None:
        manifest["selection"]["excluded_discovery_manifest"] = excluded_manifest
        manifest["selection"]["discovery_exclusion_rule"] = "sample_id_only"
    validate_manifest(manifest)
    validate_selection_manifest(manifest)
    return manifest


def validate_selection_manifest(manifest: dict[str, Any]) -> None:
    samples = manifest["samples"]
    expected_count = sum(manifest["selection"]["bin_quotas"].values())
    if len(samples) != expected_count:
        raise ValueError(
            f"expected {expected_count} selected samples, found {len(samples)}"
        )
    for field, getter in (
        ("cluster", lambda sample: sample["runs_n_poses"]["cluster"]),
        ("PDB", lambda sample: sample["structure_reference"]["pdb_id"]),
        ("ligand", lambda sample: sample["ligand"]["input_smiles"]),
    ):
        values = [getter(sample) for sample in samples]
        if len(set(values)) != len(values):
            raise ValueError(f"selected {field} values are not unique")
    observed = {
        similarity_bin_label(index): 0
        for index in range(len(SIMILARITY_EDGES) - 1)
    }
    for sample in samples:
        observed[sample["runs_n_poses"]["similarity_bin"]] += 1
    expected = manifest["selection"]["bin_quotas"]
    if observed != expected:
        raise ValueError(f"bin counts differ: observed={observed}, expected={expected}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--exclude-manifest", type=Path)
    parser.add_argument(
        "--bin-quotas",
        type=int,
        nargs=len(BIN_QUOTAS),
        default=BIN_QUOTAS,
    )
    parser.add_argument("--seed", type=int, default=SELECTION_SEED)
    parser.add_argument("--manifest-id", default="rnp_nesso1_pilot100_sep2021")
    parser.add_argument(
        "--purpose",
        default="Blinded 100-system pilot for Nesso-1 distogram and pocket generalization across Runs N' Poses familiarity bins",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_files = {args.annotations.name: sha256(args.annotations), args.inputs.name: sha256(args.inputs)}
    for name, expected in EXPECTED_SHA256.items():
        if source_files.get(name) != expected:
            raise ValueError(f"unexpected checksum for {name}: {source_files.get(name)}")

    annotations = pd.read_csv(args.annotations)
    with args.inputs.open(encoding="utf-8") as handle:
        inputs = json.load(handle)
    candidates = prepare_candidates(annotations, inputs)
    excluded_manifest_label = None
    if args.exclude_manifest is not None:
        excluded_manifest = json.loads(args.exclude_manifest.read_text())
        excluded_sample_ids = {
            sample["sample_id"] for sample in excluded_manifest["samples"]
        }
        candidates = candidates[~candidates["group_key"].isin(excluded_sample_ids)].copy()
        excluded_manifest_label = str(args.exclude_manifest)
    bin_quotas = tuple(args.bin_quotas)
    selected = select_diverse(
        candidates,
        seed=args.seed,
        bin_quotas=bin_quotas,
    )
    manifest = build_manifest(
        selected,
        source_files,
        manifest_id=args.manifest_id,
        purpose=args.purpose,
        seed=args.seed,
        bin_quotas=bin_quotas,
        excluded_manifest=excluded_manifest_label,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"eligible candidates: {len(candidates)}")
    print(f"selected samples: {len(selected)}")
    print(selected["similarity_bin"].value_counts().sort_index().to_string())
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
