#!/usr/bin/env python3
"""Analyze a small RFdiffusion backbone ensemble with explicit metric semantics.

The script compares intrinsic binder trace shape separately from target-relative
pose. It accepts unequal binder lengths and writes compact, reviewable CSV/JSON
outputs plus a browser-native JavaScript asset for the local tutorial.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pickle
import re
from pathlib import Path
from typing import Any

import numpy as np


BACKBONE_ATOMS = {"N", "CA", "C", "O"}
HOTSPOTS = (59, 83, 91)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_pdb(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            atoms.append(
                {
                    "atom": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "residue": int(line[22:26]),
                    "coord": np.asarray(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                        dtype=float,
                    ),
                }
            )
    if not atoms:
        raise ValueError(f"No atoms found in {path}")
    return atoms


def coords(
    atoms: list[dict[str, Any]], chain: str, atom_names: set[str] | None = None
) -> np.ndarray:
    selected = [
        atom["coord"]
        for atom in atoms
        if atom["chain"] == chain
        and (atom_names is None or atom["atom"] in atom_names)
    ]
    if not selected:
        raise ValueError(f"No requested atoms found for chain {chain}")
    return np.asarray(selected, dtype=float)


def ca_residue_map(atoms: list[dict[str, Any]], chain: str) -> dict[int, np.ndarray]:
    return {
        atom["residue"]: atom["coord"]
        for atom in atoms
        if atom["chain"] == chain and atom["atom"] == "CA"
    }


def fit_transform(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return row-vector rotation and translation mapping mobile to reference."""
    if mobile.shape != reference.shape:
        raise ValueError(f"Alignment shapes differ: {mobile.shape} and {reference.shape}")
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right_t))
    rotation = left @ correction @ right_t
    translation = reference_center - mobile_center @ rotation
    return rotation, translation


def apply_transform(xyz: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return xyz @ rotation + translation


def rmsd(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((first - second) ** 2, axis=1))))


def radius_of_gyration(xyz: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((xyz - xyz.mean(axis=0)) ** 2, axis=1))))


def resample_trace(xyz: np.ndarray, count: int = 64) -> np.ndarray:
    segment_lengths = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    if cumulative[-1] <= 0:
        raise ValueError("Cannot resample a zero-length trace")
    sample_positions = np.linspace(0.0, cumulative[-1], count)
    return np.column_stack(
        [np.interp(sample_positions, cumulative, xyz[:, axis]) for axis in range(3)]
    )


def intrinsic_trace_rmsd(first: np.ndarray, second: np.ndarray) -> float:
    first_sampled = resample_trace(first)
    second_sampled = resample_trace(second)
    rotation, translation = fit_transform(second_sampled, first_sampled)
    return rmsd(first_sampled, apply_transform(second_sampled, rotation, translation))


def symmetric_chamfer(first: np.ndarray, second: np.ndarray) -> float:
    distances = np.linalg.norm(first[:, None, :] - second[None, :, :], axis=-1)
    return float((distances.min(axis=0).mean() + distances.min(axis=1).mean()) / 2.0)


def jaccard(first: set[int], second: set[int]) -> float:
    union = first | second
    return float(len(first & second) / len(union)) if union else 1.0


def classical_mds(distance_matrix: np.ndarray) -> list[list[float]]:
    count = distance_matrix.shape[0]
    centering = np.eye(count) - np.ones((count, count)) / count
    gram = -0.5 * centering @ (distance_matrix**2) @ centering
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = np.maximum(values[order[:2]], 0)
    embedded = vectors[:, order[:2]] * np.sqrt(values)
    # Eigenvectors have arbitrary signs. Fix each axis for stable exports.
    for axis in range(embedded.shape[1]):
        pivot = int(np.argmax(np.abs(embedded[:, axis])))
        if embedded[pivot, axis] < 0:
            embedded[:, axis] *= -1
    return np.round(embedded, 4).tolist()


def off_diagonal_summary(matrix: np.ndarray) -> dict[str, float]:
    values = matrix[np.triu_indices_from(matrix, k=1)]
    return {
        "minimum": round(float(values.min()), 4),
        "median": round(float(np.median(values)), 4),
        "mean": round(float(values.mean()), 4),
        "maximum": round(float(values.max()), 4),
    }


def extreme_pair(matrix: np.ndarray, seeds: list[int], mode: str) -> dict[str, Any]:
    triangle = np.triu_indices_from(matrix, k=1)
    values = matrix[triangle]
    flat_index = int(np.argmin(values) if mode == "minimum" else np.argmax(values))
    first, second = int(triangle[0][flat_index]), int(triangle[1][flat_index])
    return {"seeds": [seeds[first], seeds[second]], "value": round(float(matrix[first, second]), 4)}


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed42-root", type=Path, required=True)
    parser.add_argument("--ensemble-root", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--web-output", type=Path, required=True)
    args = parser.parse_args()

    paths: list[tuple[int, Path, Path]] = [
        (
            42,
            args.seed42_root / "design_ppi_42.pdb",
            args.seed42_root / "design_ppi_42.trb",
        )
    ]
    paths.extend(
        (seed, args.ensemble_root / f"design_ppi_{seed}.pdb", args.ensemble_root / f"design_ppi_{seed}.trb")
        for seed in range(43, 52)
    )
    for _, pdb_path, trb_path in paths:
        if not pdb_path.is_file() or not trb_path.is_file():
            raise FileNotFoundError(f"Missing PDB/TRB pair: {pdb_path}, {trb_path}")

    seeds = [item[0] for item in paths]
    parsed = [parse_pdb(item[1]) for item in paths]
    reference_target = coords(parsed[0], "A", {"CA"})
    if len(reference_target) != 150:
        raise ValueError(f"Seed 42 target has {len(reference_target)} rather than 150 residues")

    designs: list[dict[str, Any]] = []
    aligned_binders: list[np.ndarray] = []
    contact_sets: list[set[int]] = []
    for (seed, pdb_path, trb_path), atoms in zip(paths, parsed):
        target_map = ca_residue_map(atoms, "A")
        binder_map = ca_residue_map(atoms, "B")
        target_ca = np.asarray(list(target_map.values()))
        binder_ca = np.asarray(list(binder_map.values()))
        if len(target_ca) != 150 or not 70 <= len(binder_ca) <= 100:
            raise ValueError(f"Seed {seed} has target/binder lengths {len(target_ca)}/{len(binder_ca)}")

        rotation, translation = fit_transform(target_ca, reference_target)
        aligned_target = apply_transform(target_ca, rotation, translation)
        aligned_binder = apply_transform(binder_ca, rotation, translation)
        aligned_binders.append(aligned_binder)
        target_alignment_rmsd = rmsd(reference_target, aligned_target)

        ca_distances = np.linalg.norm(binder_ca[:, None, :] - target_ca[None, :, :], axis=-1)
        target_residue_numbers = list(target_map.keys())
        contact_set = {
            target_residue_numbers[index]
            for index in np.where(ca_distances.min(axis=0) <= 10.0)[0]
        }
        contact_sets.append(contact_set)

        hotspot_distances = {
            f"A{residue}": float(np.linalg.norm(binder_ca - target_map[residue], axis=1).min())
            for residue in HOTSPOTS
        }
        target_backbone = coords(atoms, "A", BACKBONE_ATOMS)
        binder_backbone = coords(atoms, "B", BACKBONE_ATOMS)
        backbone_distances = np.linalg.norm(
            binder_backbone[:, None, :] - target_backbone[None, :, :], axis=-1
        )

        # This is a trusted, locally generated file from the pinned checkout.
        # Never unpickle untrusted external model artifacts.
        with trb_path.open("rb") as handle:
            metadata = pickle.load(handle)
        confidence = np.asarray(metadata["plddt"], dtype=float)
        if confidence.shape[0] != 50 or confidence.shape[1] != len(target_ca) + len(binder_ca):
            raise ValueError(f"Seed {seed} has unexpected confidence shape {confidence.shape}")

        designs.append(
            {
                "seed": seed,
                "length": len(binder_ca),
                "binderConfidence": round(float(confidence[-1, : len(binder_ca)].mean()), 4),
                "radiusOfGyrationAngstrom": round(radius_of_gyration(binder_ca), 3),
                "contactResiduesWithin10A": len(contact_set),
                "contactResidues": sorted(contact_set),
                "hotspotMinCaDistanceAngstrom": {
                    label: round(value, 3) for label, value in hotspot_distances.items()
                },
                "backbonePairsBelow2A": int(np.count_nonzero(backbone_distances < 2.0)),
                "minimumInterchainBackboneDistanceAngstrom": round(float(backbone_distances.min()), 3),
                "targetAlignmentRmsdAngstrom": round(target_alignment_rmsd, 4),
                "modelReportedSeconds": round(float(metadata["time"]), 2),
                "gpu": metadata["device"],
                "pdbSha256": sha256(pdb_path),
                "trbSha256": sha256(trb_path),
                "alignedBinderCa": np.round(aligned_binder, 3).tolist(),
            }
        )

    count = len(designs)
    center_distance = np.zeros((count, count), dtype=float)
    pose_chamfer = np.zeros((count, count), dtype=float)
    contact_jaccard = np.eye(count, dtype=float)
    shape_rmsd = np.zeros((count, count), dtype=float)
    pair_rows: list[dict[str, Any]] = []
    for first in range(count):
        for second in range(first + 1, count):
            center_value = float(
                np.linalg.norm(aligned_binders[first].mean(axis=0) - aligned_binders[second].mean(axis=0))
            )
            chamfer_value = symmetric_chamfer(aligned_binders[first], aligned_binders[second])
            jaccard_value = jaccard(contact_sets[first], contact_sets[second])
            shape_value = intrinsic_trace_rmsd(aligned_binders[first], aligned_binders[second])
            for matrix, value in (
                (center_distance, center_value),
                (pose_chamfer, chamfer_value),
                (contact_jaccard, jaccard_value),
                (shape_rmsd, shape_value),
            ):
                matrix[first, second] = matrix[second, first] = value
            pair_rows.append(
                {
                    "seed_i": seeds[first],
                    "seed_j": seeds[second],
                    "center_of_mass_distance_A": round(center_value, 4),
                    "pose_chamfer_distance_A": round(chamfer_value, 4),
                    "contact_jaccard_similarity": round(jaccard_value, 4),
                    "intrinsic_trace_rmsd_A": round(shape_value, 4),
                }
            )

    matrices = {
        "centerOfMassDistanceAngstrom": np.round(center_distance, 4).tolist(),
        "poseChamferDistanceAngstrom": np.round(pose_chamfer, 4).tolist(),
        "contactJaccardSimilarity": np.round(contact_jaccard, 4).tolist(),
        "intrinsicTraceRmsdAngstrom": np.round(shape_rmsd, 4).tolist(),
    }
    summary = {
        "length": {
            "minimum": min(design["length"] for design in designs),
            "maximum": max(design["length"] for design in designs),
            "mean": round(float(np.mean([design["length"] for design in designs])), 2),
        },
        "binderConfidence": {
            "minimum": min(design["binderConfidence"] for design in designs),
            "maximum": max(design["binderConfidence"] for design in designs),
            "mean": round(float(np.mean([design["binderConfidence"] for design in designs])), 4),
        },
        "shapeRmsdAngstrom": off_diagonal_summary(shape_rmsd),
        "poseChamferAngstrom": off_diagonal_summary(pose_chamfer),
        "centerOfMassDistanceAngstrom": off_diagonal_summary(center_distance),
        "contactJaccardSimilarity": off_diagonal_summary(contact_jaccard),
        "closestShapePair": extreme_pair(shape_rmsd, seeds, "minimum"),
        "mostDifferentShapePair": extreme_pair(shape_rmsd, seeds, "maximum"),
        "closestPosePair": extreme_pair(pose_chamfer, seeds, "minimum"),
        "mostDifferentPosePair": extreme_pair(pose_chamfer, seeds, "maximum"),
        "designsWithBackbonePairsBelow2A": sum(design["backbonePairsBelow2A"] > 0 for design in designs),
    }

    export = {
        "schemaVersion": 1,
        "experiment": "exp004_rfdiffusion_insr_structural_diversity",
        "status": "ten RFdiffusion backbones complete; no sequence or affinity validation for seeds 43-51",
        "seeds": seeds,
        "target": {"chain": "A", "ca": np.round(reference_target, 3).tolist(), "hotspots": list(HOTSPOTS)},
        "designs": designs,
        "matrices": matrices,
        "embeddings": {
            "intrinsicShape": classical_mds(shape_rmsd),
            "targetRelativePose": classical_mds(pose_chamfer),
        },
        "summary": summary,
        "definitions": {
            "intrinsicTraceRmsd": "64-point arc-length-resampled C-alpha trace RMSD after independent Kabsch superposition; descriptive and not TM-score",
            "poseChamfer": "symmetric nearest-neighbor binder C-alpha distance after aligning the 150-residue target",
            "contactJaccard": "Jaccard similarity between target-residue sets with a binder C-alpha within 10 angstrom",
        },
    }

    per_design_rows = []
    for design in designs:
        per_design_rows.append(
            {
                "seed": design["seed"],
                "length_residues": design["length"],
                "binder_confidence": design["binderConfidence"],
                "ca_radius_of_gyration_A": design["radiusOfGyrationAngstrom"],
                "target_contact_residues_10A": design["contactResiduesWithin10A"],
                "A59_min_ca_distance_A": design["hotspotMinCaDistanceAngstrom"]["A59"],
                "A83_min_ca_distance_A": design["hotspotMinCaDistanceAngstrom"]["A83"],
                "A91_min_ca_distance_A": design["hotspotMinCaDistanceAngstrom"]["A91"],
                "backbone_pairs_below_2A": design["backbonePairsBelow2A"],
                "minimum_interchain_backbone_distance_A": design["minimumInterchainBackboneDistanceAngstrom"],
                "target_alignment_rmsd_A": design["targetAlignmentRmsdAngstrom"],
                "model_reported_seconds": design["modelReportedSeconds"],
            }
        )
    write_csv(args.report_dir / "per_design.csv", list(per_design_rows[0]), per_design_rows)
    write_csv(args.report_dir / "pairwise.csv", list(pair_rows[0]), pair_rows)
    (args.report_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    args.web_output.parent.mkdir(parents=True, exist_ok=True)
    args.web_output.write_text(
        "// Generated by scripts/analyze_rfdiffusion_diversity.py\n"
        f"window.RFDIFFUSION_DIVERSITY={json.dumps(export, separators=(',', ':'), allow_nan=False)};\n"
    )
    print(json.dumps(summary, indent=2))
    print(f"Wrote {args.report_dir / 'per_design.csv'} ({len(designs)} rows)")
    print(f"Wrote {args.report_dir / 'pairwise.csv'} ({len(pair_rows)} rows)")
    print(f"Wrote {args.web_output} ({args.web_output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
