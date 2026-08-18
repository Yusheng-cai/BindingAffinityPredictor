#!/usr/bin/env python3
"""Compare matched RFdiffusion backbones with and without hotspot conditioning.

The target is aligned before target-relative pose comparisons. Intrinsic binder
shape is evaluated separately. These are geometry metrics, not affinity scores.
"""

from __future__ import annotations

import argparse
import csv
import json
import pickle
from pathlib import Path
from typing import Any

import numpy as np

from analyze_rfdiffusion_diversity import (
    BACKBONE_ATOMS,
    HOTSPOTS,
    apply_transform,
    ca_residue_map,
    coords,
    fit_transform,
    intrinsic_trace_rmsd,
    jaccard,
    parse_pdb,
    rmsd,
    sha256,
    symmetric_chamfer,
)


def contact_set_and_distances(
    target_map: dict[int, np.ndarray], binder_ca: np.ndarray, threshold: float = 10.0
) -> tuple[set[int], np.ndarray]:
    """Return contacted target residues and the full binder-target C-alpha matrix."""
    target_residues = list(target_map)
    target_ca = np.asarray(list(target_map.values()), dtype=float)
    distances = np.linalg.norm(binder_ca[:, None, :] - target_ca[None, :, :], axis=-1)
    contacts = {
        target_residues[index]
        for index in np.where(distances.min(axis=0) <= threshold)[0]
    }
    return contacts, distances


def summarize_design(
    pdb_path: Path,
    trb_path: Path,
    reference_target: np.ndarray,
    contact_threshold: float = 10.0,
) -> dict[str, Any]:
    atoms = parse_pdb(pdb_path)
    target_map = ca_residue_map(atoms, "A")
    binder_map = ca_residue_map(atoms, "B")
    target_ca = np.asarray(list(target_map.values()), dtype=float)
    binder_ca = np.asarray(list(binder_map.values()), dtype=float)
    if len(target_ca) != 150 or not 70 <= len(binder_ca) <= 100:
        raise ValueError(
            f"Unexpected target/binder lengths in {pdb_path}: "
            f"{len(target_ca)}/{len(binder_ca)}"
        )
    missing_hotspots = sorted(set(HOTSPOTS) - set(target_map))
    if missing_hotspots:
        raise ValueError(f"Missing nominal hotspots {missing_hotspots} in {pdb_path}")

    rotation, translation = fit_transform(target_ca, reference_target)
    aligned_target = apply_transform(target_ca, rotation, translation)
    aligned_binder = apply_transform(binder_ca, rotation, translation)
    contacts, ca_distances = contact_set_and_distances(
        target_map, binder_ca, threshold=contact_threshold
    )
    hotspot_distances = {
        f"A{residue}": float(np.linalg.norm(binder_ca - target_map[residue], axis=1).min())
        for residue in HOTSPOTS
    }
    hotspot_centroid = np.asarray([target_map[residue] for residue in HOTSPOTS]).mean(axis=0)
    nearest_flat = int(np.argmin(ca_distances))
    _, nearest_target_index = np.unravel_index(nearest_flat, ca_distances.shape)
    target_residues = list(target_map)

    target_backbone = coords(atoms, "A", BACKBONE_ATOMS)
    binder_backbone = coords(atoms, "B", BACKBONE_ATOMS)
    backbone_distances = np.linalg.norm(
        binder_backbone[:, None, :] - target_backbone[None, :, :], axis=-1
    )

    # Trusted local artifact from the pinned RFdiffusion checkout only.
    with trb_path.open("rb") as handle:
        metadata = pickle.load(handle)
    confidence = np.asarray(metadata["plddt"], dtype=float)
    if confidence.shape != (50, len(target_ca) + len(binder_ca)):
        raise ValueError(f"Unexpected confidence shape {confidence.shape} in {trb_path}")

    return {
        "length": len(binder_ca),
        "contacts": contacts,
        "hotspot_distances": hotspot_distances,
        "minimum_hotspot_distance": min(hotspot_distances.values()),
        "binder_center_to_hotspot_centroid": float(
            np.linalg.norm(binder_ca.mean(axis=0) - hotspot_centroid)
        ),
        "nearest_target_residue": target_residues[nearest_target_index],
        "minimum_ca_distance": float(ca_distances.min()),
        "minimum_backbone_distance": float(backbone_distances.min()),
        "backbone_pairs_below_2A": int(np.count_nonzero(backbone_distances < 2.0)),
        "final_binder_confidence": float(confidence[-1, : len(binder_ca)].mean()),
        "target_alignment_rmsd": rmsd(reference_target, aligned_target),
        "aligned_binder": aligned_binder,
        "pdb_sha256": sha256(pdb_path),
        "trb_sha256": sha256(trb_path),
    }


def locate_guided(seed: int, seed42_root: Path, ensemble_root: Path) -> tuple[Path, Path]:
    root = seed42_root if seed == 42 else ensemble_root
    return root / f"design_ppi_{seed}.pdb", root / f"design_ppi_{seed}.trb"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def rounded_summary(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "minimum": round(float(array.min()), 4),
        "median": round(float(np.median(array)), 4),
        "mean": round(float(array.mean()), 4),
        "maximum": round(float(array.max()), 4),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--guided-seed42-root", type=Path, required=True)
    parser.add_argument("--guided-ensemble-root", type=Path, required=True)
    parser.add_argument("--unguided-root", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()

    seeds = list(range(42, 52))
    guided_paths = [locate_guided(seed, args.guided_seed42_root, args.guided_ensemble_root) for seed in seeds]
    unguided_paths = [
        (
            args.unguided_root / f"design_ppi_{seed}.pdb",
            args.unguided_root / f"design_ppi_{seed}.trb",
        )
        for seed in seeds
    ]
    for pair in guided_paths + unguided_paths:
        if not all(path.is_file() for path in pair):
            raise FileNotFoundError(f"Missing PDB/TRB pair: {pair}")

    reference_atoms = parse_pdb(guided_paths[0][0])
    reference_target = coords(reference_atoms, "A", {"CA"})
    guided = [summarize_design(*paths, reference_target) for paths in guided_paths]
    unguided = [summarize_design(*paths, reference_target) for paths in unguided_paths]

    rows: list[dict[str, Any]] = []
    for seed, with_hint, without_hint in zip(seeds, guided, unguided):
        rows.append(
            {
                "seed": seed,
                "guided_length": with_hint["length"],
                "unguided_length": without_hint["length"],
                "guided_min_hotspot_CA_distance_A": round(with_hint["minimum_hotspot_distance"], 4),
                "unguided_min_hotspot_CA_distance_A": round(without_hint["minimum_hotspot_distance"], 4),
                "guided_center_to_hotspot_centroid_A": round(with_hint["binder_center_to_hotspot_centroid"], 4),
                "unguided_center_to_hotspot_centroid_A": round(without_hint["binder_center_to_hotspot_centroid"], 4),
                "guided_contact_count_10A": len(with_hint["contacts"]),
                "unguided_contact_count_10A": len(without_hint["contacts"]),
                "matched_contact_jaccard": round(jaccard(with_hint["contacts"], without_hint["contacts"]), 4),
                "matched_pose_chamfer_A": round(symmetric_chamfer(with_hint["aligned_binder"], without_hint["aligned_binder"]), 4),
                "matched_center_displacement_A": round(float(np.linalg.norm(with_hint["aligned_binder"].mean(axis=0) - without_hint["aligned_binder"].mean(axis=0))), 4),
                "matched_intrinsic_trace_rmsd_A": round(intrinsic_trace_rmsd(with_hint["aligned_binder"], without_hint["aligned_binder"]), 4),
                "guided_nearest_target_residue": with_hint["nearest_target_residue"],
                "unguided_nearest_target_residue": without_hint["nearest_target_residue"],
                "guided_min_backbone_distance_A": round(with_hint["minimum_backbone_distance"], 4),
                "unguided_min_backbone_distance_A": round(without_hint["minimum_backbone_distance"], 4),
                "guided_backbone_pairs_below_2A": with_hint["backbone_pairs_below_2A"],
                "unguided_backbone_pairs_below_2A": without_hint["backbone_pairs_below_2A"],
                "guided_final_confidence": round(with_hint["final_binder_confidence"], 4),
                "unguided_final_confidence": round(without_hint["final_binder_confidence"], 4),
                "guided_pdb_sha256": with_hint["pdb_sha256"],
                "unguided_pdb_sha256": without_hint["pdb_sha256"],
                "guided_trb_sha256": with_hint["trb_sha256"],
                "unguided_trb_sha256": without_hint["trb_sha256"],
            }
        )

    summary = {
        "paired_seeds": seeds,
        "pairs_with_guided_minimum_hotspot_distance_smaller": sum(
            with_hint["minimum_hotspot_distance"] < without_hint["minimum_hotspot_distance"]
            for with_hint, without_hint in zip(guided, unguided)
        ),
        "pairs_total": len(seeds),
        "guided_min_hotspot_CA_distance_A": rounded_summary([d["minimum_hotspot_distance"] for d in guided]),
        "unguided_min_hotspot_CA_distance_A": rounded_summary([d["minimum_hotspot_distance"] for d in unguided]),
        "guided_center_to_hotspot_centroid_A": rounded_summary([d["binder_center_to_hotspot_centroid"] for d in guided]),
        "unguided_center_to_hotspot_centroid_A": rounded_summary([d["binder_center_to_hotspot_centroid"] for d in unguided]),
        "matched_contact_jaccard": rounded_summary([row["matched_contact_jaccard"] for row in rows]),
        "matched_pose_chamfer_A": rounded_summary([row["matched_pose_chamfer_A"] for row in rows]),
        "matched_center_displacement_A": rounded_summary([row["matched_center_displacement_A"] for row in rows]),
        "matched_intrinsic_trace_rmsd_A": rounded_summary([row["matched_intrinsic_trace_rmsd_A"] for row in rows]),
        "guided_designs_with_backbone_pairs_below_2A": sum(d["backbone_pairs_below_2A"] > 0 for d in guided),
        "unguided_designs_with_backbone_pairs_below_2A": sum(d["backbone_pairs_below_2A"] > 0 for d in unguided),
        "interpretation": "Geometry-only hotspot-conditioning ablation; no sequence, affinity, or experimental validation.",
    }

    args.report_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.report_dir / "paired_designs.csv", rows)
    (args.report_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {args.report_dir / 'paired_designs.csv'} ({len(rows)} paired rows)")


if __name__ == "__main__":
    main()
