#!/usr/bin/env python3
"""Freeze a diverse post-Boltz-cutoff Runs N' Poses cohort for two models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from affinity_benchmark.data.manifest import validate_manifest
from select_rnp_nesso_pilot import (
    PROPERTY_COLUMNS,
    _combined_distance,
    _tie_break,
    prepare_candidates,
    sha256,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ANNOTATIONS = REPOSITORY_ROOT / "data/raw/runs_n_poses_metadata/annotations.csv"
DEFAULT_INPUTS = REPOSITORY_ROOT / "data/raw/runs_n_poses_metadata/inputs.json"
DEFAULT_SIMILARITY = (
    REPOSITORY_ROOT / "data/processed/runs_n_poses/boltz2_2023_similarity.csv"
)
DEFAULT_OUTPUT = (
    REPOSITORY_ROOT / "data/manifests/rnp_boltz2_nesso1_postcutoff50.json"
)
SOURCE_PARQUET = (
    REPOSITORY_ROOT
    / "data/raw/runs_n_poses_metadata/all_similarity_scores.parquet"
)

SOURCE_DOI = "10.5281/zenodo.14794785"
SOURCE_RECORD = "18366081"
SOURCE_REPOSITORY_REVISION = "197fafc60a5cef3a9e7f8a4b0dac7c965eed3839"
CUTOFF = "2023-06-01"
SELECTION_SEED = 2026082801
SIMILARITY_EDGES = (0.0, 20.0, 40.0, 60.0, 80.0, 100.000001)
BIN_QUOTAS = (10, 10, 10, 10, 10)


def bin_index(value: float) -> int:
    if not np.isfinite(value) or value < 0 or value > 100:
        raise ValueError(f"familiarity must be finite and in [0, 100], got {value}")
    return min(
        int(np.searchsorted(SIMILARITY_EDGES, value, side="right") - 1),
        len(SIMILARITY_EDGES) - 2,
    )


def bin_label(index: int) -> str:
    right = 100 if index == len(SIMILARITY_EDGES) - 2 else int(SIMILARITY_EDGES[index + 1])
    return f"{int(SIMILARITY_EDGES[index])}-{right}"


def prepare_eligible(
    annotations: pd.DataFrame,
    inputs: dict[str, Any],
    cutoff_similarity: pd.DataFrame,
) -> pd.DataFrame:
    candidates = prepare_candidates(annotations, inputs)
    candidates = candidates[
        pd.to_datetime(candidates["release_date"]) > pd.Timestamp(CUTOFF)
    ].copy()

    old_similarity = candidates.loc[
        :,
        ["group_key", "target_system", "target_release_date", "sucos_shape_pocket_qcov"],
    ].rename(
        columns={
            "target_system": "nesso1_closest_pre_cutoff_system",
            "target_release_date": "nesso1_closest_system_release_date",
            "sucos_shape_pocket_qcov": "nesso1_sep2021_similarity_score_0_to_100",
        }
    )
    candidates = candidates.drop(
        columns=[
            "target_system",
            "target_release_date",
            "sucos_shape_pocket_qcov",
            "similarity_bin_index",
            "similarity_bin",
        ]
    )
    candidates = candidates.merge(old_similarity, on="group_key", validate="one_to_one")
    candidates = candidates.merge(
        cutoff_similarity.loc[
            :,
            [
                "group_key",
                "closest_pre_cutoff_system",
                "closest_system_release_date",
                "similarity_score_0_to_100",
            ],
        ],
        on="group_key",
        how="inner",
        validate="one_to_one",
    )
    candidates["similarity_bin_index"] = candidates[
        "similarity_score_0_to_100"
    ].map(bin_index)
    candidates["similarity_bin"] = candidates["similarity_bin_index"].map(bin_label)

    properties = candidates.loc[:, PROPERTY_COLUMNS].astype(float)
    scale = properties.std(axis=0, ddof=0).replace(0, 1.0)
    normalized = (properties - properties.mean(axis=0)) / scale
    candidates["property_vector"] = list(normalized.to_numpy())
    return candidates.sort_values("group_key").reset_index(drop=True)


def select_diverse(candidates: pd.DataFrame) -> pd.DataFrame:
    selected_indices: list[int] = []
    selected_clusters: set[str] = set()
    selected_pdbs: set[str] = set()
    selected_smiles: set[str] = set()
    counts = [0] * len(BIN_QUOTAS)

    while counts != list(BIN_QUOTAS):
        made_progress = False
        for current_bin, quota in enumerate(BIN_QUOTAS):
            if counts[current_bin] >= quota:
                continue
            eligible = candidates[
                (candidates["similarity_bin_index"] == current_bin)
                & ~candidates.index.isin(selected_indices)
                & ~candidates["cluster"].isin(selected_clusters)
                & ~candidates["entry_pdb_id"].str.upper().isin(selected_pdbs)
                & ~candidates["canonical_smiles"].isin(selected_smiles)
            ]
            if eligible.empty:
                raise RuntimeError(
                    f"cannot fill bin {bin_label(current_bin)} with frozen diversity rules"
                )

            if not selected_indices:
                vectors = np.vstack(eligible["property_vector"])
                center = np.median(vectors, axis=0)
                scores = {
                    index: -float(np.linalg.norm(row["property_vector"] - center))
                    for index, row in eligible.iterrows()
                }
            else:
                chosen = [candidates.loc[index] for index in selected_indices]
                scores = {
                    index: min(_combined_distance(row, prior) for prior in chosen)
                    for index, row in eligible.iterrows()
                }
            selected_index = max(
                eligible.index,
                key=lambda index: (
                    scores[index],
                    _tie_break(eligible.loc[index, "group_key"], SELECTION_SEED),
                ),
            )
            row = candidates.loc[selected_index]
            selected_indices.append(int(selected_index))
            selected_clusters.add(str(row["cluster"]))
            selected_pdbs.add(str(row["entry_pdb_id"]).upper())
            selected_smiles.add(str(row["canonical_smiles"]))
            counts[current_bin] += 1
            made_progress = True
        if not made_progress:
            raise RuntimeError("selection made no progress")

    selected = candidates.loc[selected_indices].copy()
    selected["selection_order"] = range(1, len(selected) + 1)
    return selected.sort_values(
        ["similarity_bin_index", "similarity_score_0_to_100", "group_key"]
    ).reset_index(drop=True)


def build_manifest(selected: pd.DataFrame, source_files: dict[str, str]) -> dict[str, Any]:
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
                    "ground_truth_coordinates_downloaded": True,
                },
                "runs_n_poses": {
                    "system_id": row.system_id,
                    "group_key": row.group_key,
                    "cluster": row.cluster,
                    "release_date": row.release_date,
                    "structural_training_cutoff": CUTOFF,
                    "closest_pre_cutoff_system": row.closest_pre_cutoff_system,
                    "closest_system_release_date": str(row.closest_system_release_date),
                    "similarity_metric": "sucos_shape_pocket_qcov",
                    "similarity_score_0_to_100": float(row.similarity_score_0_to_100),
                    "similarity_bin": row.similarity_bin,
                    "nesso1_sep2021_familiarity": {
                        "structural_training_cutoff": "2021-09-30",
                        "closest_pre_cutoff_system": row.nesso1_closest_pre_cutoff_system,
                        "closest_system_release_date": str(row.nesso1_closest_system_release_date),
                        "similarity_score_0_to_100": float(
                            row.nesso1_sep2021_similarity_score_0_to_100
                        ),
                    },
                    "pocket_residues": int(row.ligand_num_pocket_residues),
                    "unique_interactions": int(row.ligand_num_unique_interactions),
                    "selection_order": int(row.selection_order),
                },
                "included": True,
            }
        )

    manifest = {
        "schema_version": 1,
        "manifest_id": "rnp_boltz2_nesso1_postcutoff50",
        "created_on": "2026-08-28",
        "purpose": (
            "Frozen paired 50-system structural benchmark for Boltz-2 and Nesso-1 "
            "across June-2023 structural familiarity"
        ),
        "source_dataset": {
            "name": "Runs N' Poses",
            "doi": SOURCE_DOI,
            "zenodo_record_id": SOURCE_RECORD,
            "repository_url": "https://github.com/plinder-org/runs-n-poses",
            "repository_revision": SOURCE_REPOSITORY_REVISION,
            "source_files": [
                {"name": name, "sha256": digest}
                for name, digest in sorted(source_files.items())
            ],
        },
        "selection": {
            "status": "frozen_before_nesso_or_boltz_inference",
            "structural_training_cutoff": CUTOFF,
            "test_release_rule": f"release_date > {CUTOFF}",
            "familiarity_metric": "sucos_shape_pocket_qcov",
            "similarity_edges": [0, 20, 40, 60, 80, 100],
            "bin_quotas": {bin_label(i): quota for i, quota in enumerate(BIN_QUOTAS)},
            "seed": SELECTION_SEED,
            "eligible_candidates": 205,
            "eligibility": {
                "proper_ligand": True,
                "exactly_one_ligand_chain": True,
                "exactly_one_proper_ligand_chain": True,
                "exactly_one_protein_chain": True,
                "protein_length_residues": [40, 900],
                "protein_alphabet": "20_standard_amino_acids",
                "ligand_heavy_atoms": [6, 50],
                "present_in_inputs_json": True,
                "non_null_june2023_familiarity": True,
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
            "blinding": (
                "Nesso predictions, Boltz predictions, and experimental coordinates "
                "were not used for selection"
            ),
        },
        "samples": samples,
    }
    validate_manifest(manifest)
    validate_selection(manifest)
    return manifest


def validate_selection(manifest: dict[str, Any]) -> None:
    samples = manifest["samples"]
    if len(samples) != 50:
        raise ValueError(f"expected 50 samples, found {len(samples)}")
    for label, values in {
        "cluster": [sample["runs_n_poses"]["cluster"] for sample in samples],
        "PDB": [sample["structure_reference"]["pdb_id"] for sample in samples],
        "ligand": [sample["ligand"]["input_smiles"] for sample in samples],
    }.items():
        if len(values) != len(set(values)):
            raise ValueError(f"selected {label} values are not unique")
    observed = {bin_label(i): 0 for i in range(len(BIN_QUOTAS))}
    for sample in samples:
        observed[sample["runs_n_poses"]["similarity_bin"]] += 1
        if sample["runs_n_poses"]["release_date"] <= CUTOFF:
            raise ValueError(f"pre-cutoff test sample: {sample['sample_id']}")
    if observed != manifest["selection"]["bin_quotas"]:
        raise ValueError(f"bin mismatch: {observed}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--similarity", type=Path, default=DEFAULT_SIMILARITY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    annotations = pd.read_csv(args.annotations)
    inputs = json.loads(args.inputs.read_text())
    cutoff_similarity = pd.read_csv(args.similarity)
    candidates = prepare_eligible(annotations, inputs, cutoff_similarity)
    if len(candidates) != 205:
        raise ValueError(f"expected 205 eligible candidates, found {len(candidates)}")
    selected = select_diverse(candidates)
    source_files = {
        args.annotations.name: sha256(args.annotations),
        args.inputs.name: sha256(args.inputs),
        args.similarity.name: sha256(args.similarity),
        SOURCE_PARQUET.name: sha256(SOURCE_PARQUET),
    }
    manifest = build_manifest(selected, source_files)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Eligible candidates: {len(candidates)}")
    print(f"Selected systems: {len(selected)}")
    print(selected["similarity_bin"].value_counts().sort_index().to_string())
    print(f"Unique sequences: {selected['protein_sequence'].nunique()}")
    print(f"Protein length range: {selected['protein_total_length'].min()}--{selected['protein_total_length'].max()}")
    print(f"Ligand heavy-atom range: {selected['ligand_num_heavy_atoms'].min()}--{selected['ligand_num_heavy_atoms'].max()}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
