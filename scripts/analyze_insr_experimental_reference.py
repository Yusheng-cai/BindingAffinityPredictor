#!/usr/bin/env python3
"""Compare RFdiffusion binders with the experimental 4OGA site-1 complex.

The experimental receptor L1 C-alpha trace (chain E, residues 6-155) is
least-squares aligned to the 150-residue RFdiffusion target. Insulin and alphaCT
are carried through that rigid transform. All reported overlaps are coarse
C-alpha geometry descriptors, not energetic or affinity estimates.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from analyze_rfdiffusion_diversity import (
    apply_transform,
    ca_residue_map,
    coords,
    fit_transform,
    jaccard,
    parse_pdb,
    rmsd,
    sha256,
)
from analyze_rfdiffusion_hotspot_ablation import contact_set_and_distances, locate_guided


EXPERIMENTAL_L1_RESIDUES = tuple(range(6, 156))
CONTACT_THRESHOLD_ANGSTROM = 10.0
COVERAGE_THRESHOLD_ANGSTROM = 8.0


def selected_ca(
    atoms: list[dict[str, Any]], chain: str, residues: set[int] | None = None
) -> np.ndarray:
    """Return ordered C-alpha coordinates, optionally restricted by residue number."""
    selected = [
        atom["coord"]
        for atom in atoms
        if atom["chain"] == chain
        and atom["atom"] == "CA"
        and (residues is None or atom["residue"] in residues)
    ]
    if not selected:
        raise ValueError(f"No requested C-alpha atoms for chain {chain}")
    return np.asarray(selected, dtype=float)


def component_contact_set(
    target_map: dict[int, np.ndarray], component_ca: np.ndarray, threshold: float
) -> set[int]:
    """Target residues with a component C-alpha no farther than ``threshold``."""
    contacts, _ = contact_set_and_distances(target_map, component_ca, threshold)
    return contacts


def component_coverage(
    component_ca: np.ndarray, binder_ca: np.ndarray, threshold: float
) -> float:
    """Fraction of experimental component C-alpha sites approached by the binder."""
    distances = np.linalg.norm(component_ca[:, None, :] - binder_ca[None, :, :], axis=-1)
    return float(np.mean(distances.min(axis=1) <= threshold))


def minimum_cross_distance(first: np.ndarray, second: np.ndarray) -> float:
    distances = np.linalg.norm(first[:, None, :] - second[None, :, :], axis=-1)
    return float(distances.min())


def load_4oga_aligned(
    pdb_path: Path, reference_target_ca: np.ndarray
) -> dict[str, Any]:
    """Load 4OGA and rigidly place its biological components in target coordinates."""
    atoms = parse_pdb(pdb_path)
    chain_e_map = ca_residue_map(atoms, "E")
    missing = sorted(set(EXPERIMENTAL_L1_RESIDUES) - set(chain_e_map))
    if missing:
        raise ValueError(f"4OGA is missing L1 C-alpha residues: {missing}")
    experimental_l1 = np.asarray(
        [chain_e_map[residue] for residue in EXPERIMENTAL_L1_RESIDUES], dtype=float
    )
    if reference_target_ca.shape != experimental_l1.shape:
        raise ValueError(
            f"Target/reference shapes differ: {reference_target_ca.shape} and "
            f"{experimental_l1.shape}"
        )
    rotation, translation = fit_transform(experimental_l1, reference_target_ca)
    aligned_l1 = apply_transform(experimental_l1, rotation, translation)
    insulin_a = selected_ca(atoms, "A")
    insulin_b = selected_ca(atoms, "B")
    insulin = np.vstack((insulin_a, insulin_b))
    alpha_ct = selected_ca(atoms, "F")
    return {
        "l1": aligned_l1,
        "insulin_a": apply_transform(insulin_a, rotation, translation),
        "insulin_b": apply_transform(insulin_b, rotation, translation),
        "insulin": apply_transform(insulin, rotation, translation),
        "alpha_ct": apply_transform(alpha_ct, rotation, translation),
        "alignment_rmsd": rmsd(reference_target_ca, aligned_l1),
        "experimental_target_map": {
            residue - 5: chain_e_map[residue] for residue in EXPERIMENTAL_L1_RESIDUES
        },
        "unaligned_insulin": insulin,
        "unaligned_alpha_ct": alpha_ct,
        "pdb_sha256": sha256(pdb_path),
    }


def summarize(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "minimum": round(float(array.min()), 4),
        "median": round(float(np.median(array)), 4),
        "mean": round(float(array.mean()), 4),
        "maximum": round(float(array.max()), 4),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experimental-pdb", type=Path, required=True)
    parser.add_argument("--guided-seed42-root", type=Path, required=True)
    parser.add_argument("--guided-ensemble-root", type=Path, required=True)
    parser.add_argument("--unguided-root", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()

    seeds = list(range(42, 52))
    guided_paths = [
        locate_guided(seed, args.guided_seed42_root, args.guided_ensemble_root)[0]
        for seed in seeds
    ]
    unguided_paths = [args.unguided_root / f"design_ppi_{seed}.pdb" for seed in seeds]
    for path in [args.experimental_pdb, *guided_paths, *unguided_paths]:
        if not path.is_file():
            raise FileNotFoundError(path)

    reference_atoms = parse_pdb(guided_paths[0])
    reference_target = coords(reference_atoms, "A", {"CA"})
    experimental = load_4oga_aligned(args.experimental_pdb, reference_target)
    insulin_contacts = component_contact_set(
        experimental["experimental_target_map"],
        experimental["unaligned_insulin"],
        CONTACT_THRESHOLD_ANGSTROM,
    )
    alpha_ct_contacts = component_contact_set(
        experimental["experimental_target_map"],
        experimental["unaligned_alpha_ct"],
        CONTACT_THRESHOLD_ANGSTROM,
    )
    combined_contacts = insulin_contacts | alpha_ct_contacts

    rows: list[dict[str, Any]] = []
    for arm, paths in (("guided", guided_paths), ("unguided", unguided_paths)):
        for seed, path in zip(seeds, paths):
            atoms = parse_pdb(path)
            target_map = ca_residue_map(atoms, "A")
            target_ca = np.asarray(list(target_map.values()), dtype=float)
            binder_ca = coords(atoms, "B", {"CA"})
            rotation, translation = fit_transform(target_ca, reference_target)
            aligned_binder = apply_transform(binder_ca, rotation, translation)
            binder_contacts, _ = contact_set_and_distances(
                target_map, binder_ca, CONTACT_THRESHOLD_ANGSTROM
            )
            rows.append(
                {
                    "arm": arm,
                    "seed": seed,
                    "binder_length": len(binder_ca),
                    "minimum_CA_distance_to_experimental_insulin_A": round(
                        minimum_cross_distance(aligned_binder, experimental["insulin"]), 4
                    ),
                    "minimum_CA_distance_to_experimental_alphaCT_A": round(
                        minimum_cross_distance(aligned_binder, experimental["alpha_ct"]), 4
                    ),
                    "experimental_insulin_CA_coverage_within_8A": round(
                        component_coverage(
                            experimental["insulin"], aligned_binder, COVERAGE_THRESHOLD_ANGSTROM
                        ),
                        4,
                    ),
                    "experimental_alphaCT_CA_coverage_within_8A": round(
                        component_coverage(
                            experimental["alpha_ct"], aligned_binder, COVERAGE_THRESHOLD_ANGSTROM
                        ),
                        4,
                    ),
                    "binder_contact_jaccard_with_experimental_insulin": round(
                        jaccard(binder_contacts, insulin_contacts), 4
                    ),
                    "binder_contact_jaccard_with_experimental_alphaCT": round(
                        jaccard(binder_contacts, alpha_ct_contacts), 4
                    ),
                    "binder_contact_jaccard_with_experimental_combined": round(
                        jaccard(binder_contacts, combined_contacts), 4
                    ),
                    "pdb_sha256": sha256(path),
                }
            )

    summary: dict[str, Any] = {
        "experimental_reference": {
            "pdb_id": "4OGA",
            "method": "X-ray diffraction",
            "resolution_angstrom": 3.5,
            "l1_chain_and_residues": "E6-E155",
            "insulin_chains": ["A", "B"],
            "alpha_ct_chain_and_modeled_residues": "F705-F719",
            "l1_alignment_rmsd_angstrom": round(experimental["alignment_rmsd"], 4),
            "pdb_sha256": experimental["pdb_sha256"],
            "insulin_modeled_CA_count": len(experimental["insulin"]),
            "alpha_ct_modeled_CA_count": len(experimental["alpha_ct"]),
            "experimental_insulin_contact_residues_local_numbering": sorted(insulin_contacts),
            "experimental_alphaCT_contact_residues_local_numbering": sorted(alpha_ct_contacts),
            "contact_definition": "L1 C-alpha within 10 A of a component C-alpha",
        },
        "metric_interpretation": {
            "coverage": "Fraction of experimental component C-alpha sites within 8 A of a generated binder C-alpha.",
            "jaccard": "Overlap between generated-binder and experimental-component L1 contact-residue sets.",
            "affinity": "None of these geometry descriptors estimates binding affinity or functional activity.",
        },
        "arms": {},
    }
    for arm in ("guided", "unguided"):
        arm_rows = [row for row in rows if row["arm"] == arm]
        summary["arms"][arm] = {
            key: summarize([float(row[key]) for row in arm_rows])
            for key in (
                "minimum_CA_distance_to_experimental_insulin_A",
                "minimum_CA_distance_to_experimental_alphaCT_A",
                "experimental_insulin_CA_coverage_within_8A",
                "experimental_alphaCT_CA_coverage_within_8A",
                "binder_contact_jaccard_with_experimental_insulin",
                "binder_contact_jaccard_with_experimental_alphaCT",
                "binder_contact_jaccard_with_experimental_combined",
            )
        }

    args.report_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.report_dir / "experimental_reference_4oga.csv", rows)
    (args.report_dir / "experimental_reference_4oga.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
