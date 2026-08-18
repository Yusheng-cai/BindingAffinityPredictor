#!/usr/bin/env python3
"""Build a paired pose-RMSD comparison from two complete compact CSV files."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def load_complete(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    incomplete = [row["sample_id"] for row in rows if row.get("status", "complete") != "complete"]
    if incomplete:
        raise ValueError(f"{path} contains incomplete samples: {incomplete}")
    result = {row["sample_id"]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"{path} contains duplicate sample IDs")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flashbind", type=Path, required=True)
    parser.add_argument("--boltz2", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    flashbind = load_complete(args.flashbind)
    boltz2 = load_complete(args.boltz2)
    if set(flashbind) != set(boltz2):
        raise ValueError("FlashBind and Boltz-2 sample IDs are not identical")

    rows = []
    for sample_id, flashbind_row in flashbind.items():
        boltz_row = boltz2[sample_id]
        flashbind_rmsd = float(flashbind_row["ligand_rmsd_pocket_A"])
        boltz2_rmsd = float(boltz_row["ligand_rmsd_pocket_A"])
        rows.append(
            {
                "sample_id": sample_id,
                "target_id": flashbind_row["target_id"],
                "pdb_id": flashbind_row["pdb_id"],
                "flashbind_rmsd_A": flashbind_rmsd,
                "boltz2_rmsd_A": boltz2_rmsd,
                "flashbind_minus_boltz2_A": flashbind_rmsd - boltz2_rmsd,
                "lower_rmsd": "flashbind"
                if flashbind_rmsd < boltz2_rmsd
                else "boltz2",
            }
        )

    flashbind_values = np.array([row["flashbind_rmsd_A"] for row in rows])
    boltz2_values = np.array([row["boltz2_rmsd_A"] for row in rows])
    differences = flashbind_values - boltz2_values
    targets = sorted({row["target_id"] for row in rows})
    summary = {
        "metric": "symmetry-corrected ligand heavy-atom RMSD after reference-pocket protein CA alignment",
        "unit": "angstrom",
        "paired_samples": len(rows),
        "descriptive_success_threshold_angstrom": 2.0,
        "flashbind_released_FABind_plus": {
            "mean": float(np.mean(flashbind_values)),
            "median": float(np.median(flashbind_values)),
            "below_2A": int(np.sum(flashbind_values < 2.0)),
        },
        "boltz2_msa1024_seed42": {
            "mean": float(np.mean(boltz2_values)),
            "median": float(np.median(boltz2_values)),
            "below_2A": int(np.sum(boltz2_values < 2.0)),
        },
        "paired_difference_flashbind_minus_boltz2_angstrom": {
            "mean": float(np.mean(differences)),
            "median": float(np.median(differences)),
        },
        "lower_rmsd_counts": {
            "flashbind": int(np.sum(differences < 0)),
            "boltz2": int(np.sum(differences > 0)),
            "ties": int(np.sum(differences == 0)),
        },
        "by_target": {},
        "interpretation": (
            "Descriptive retrospective paired comparison; the two pose sources use "
            "different receptor inputs and neither is a prospective time split."
        ),
    }
    for target in targets:
        target_rows = [row for row in rows if row["target_id"] == target]
        flash = np.array([row["flashbind_rmsd_A"] for row in target_rows])
        boltz = np.array([row["boltz2_rmsd_A"] for row in target_rows])
        summary["by_target"][target] = {
            "n": len(target_rows),
            "flashbind_median_rmsd_A": float(np.median(flash)),
            "boltz2_median_rmsd_A": float(np.median(boltz)),
        }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "paired_pose_rmsd.csv"
    json_path = args.output_dir / "paired_pose_comparison.json"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
