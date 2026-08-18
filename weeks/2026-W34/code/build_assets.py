#!/usr/bin/env python3
"""Build compact, tracked assets for the 2026-W34 report."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
from pathlib import Path


TARGET_ORDER = ("cdk2", "tyk2", "jnk1", "p38")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def zscores(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / len(values)
    scale = math.sqrt(variance)
    if scale == 0:
        raise ValueError("cannot standardize a constant series")
    return [(value - mean) / scale for value in values]


def write_standardized_scatter(
    *,
    model: str,
    manifest: dict,
    prediction_path: Path,
    figures: Path,
) -> tuple[list[Path], int]:
    with prediction_path.open(newline="", encoding="utf-8") as handle:
        predictions = {row["sample_id"]: row for row in csv.DictReader(handle)}

    rows: list[dict[str, object]] = []
    for target in TARGET_ORDER:
        samples = [sample for sample in manifest["samples"] if sample["target_id"] == target]
        observed = [float(sample["measurement"]["log10_value_uM"]) for sample in samples]
        predicted = [float(predictions[sample["sample_id"]]["affinity_pred_value"]) for sample in samples]
        for sample, observed_z, predicted_z in zip(
            samples, zscores(observed), zscores(predicted), strict=True
        ):
            rows.append(
                {
                    "sample_id": sample["sample_id"],
                    "target": target,
                    "observed_z": f"{observed_z:.8f}",
                    "predicted_z": f"{predicted_z:.8f}",
                }
            )

    output_paths = []
    scatter_path = figures / f"{model}_standardized_scatter.csv"
    with scatter_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    output_paths.append(scatter_path)

    for target in TARGET_ORDER:
        target_path = figures / f"{model}_{target}.csv"
        target_rows = [row for row in rows if row["target"] == target]
        with target_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(target_rows)
        output_paths.append(target_path)
    return output_paths, len(rows)


def write_pose_comparison(*, source: Path, figures: Path) -> tuple[list[Path], int]:
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    output_paths = []
    for target in TARGET_ORDER:
        target_rows = [row for row in rows if row["target_id"] == target]
        output = figures / f"pose_rmsd_{target}.csv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "sample_id",
                    "target_id",
                    "pdb_id",
                    "boltz2_rmsd_A",
                    "flashbind_rmsd_A",
                ),
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(target_rows)
        output_paths.append(output)
    return output_paths, len(rows)


def main() -> int:
    root = Path(__file__).resolve().parents[3]
    report = root / "weeks/2026-W34/report"
    figures = report / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    manifest_path = root / "data/manifests/fepplus4_87.json"
    nesso_prediction_path = (
        root
        / "runs/exp007_fepplus4_boltz2_nesso1/nesso1/seed42/predictions.csv"
    )
    boltz_prediction_path = (
        root
        / "runs/exp007_fepplus4_boltz2_nesso1/boltz2_msa1024/seed42/predictions.csv"
    )
    flashbind_prediction_path = (
        root
        / "runs/exp009_flashbind_fepplus4_released_poses/flashbind_released_poses/seed42/full87/predictions.csv"
    )
    nesso_metrics_path = (
        root / "reports/exp007_fepplus4_boltz2_nesso1/nesso1_metrics.json"
    )
    comparison_metrics_path = (
        root
        / "reports/exp007_fepplus4_boltz2_nesso1/boltz2_msa1024_vs_nesso1_metrics.json"
    )
    flashbind_metrics_path = (
        root / "reports/exp009_flashbind_fepplus4_released_poses/flashbind_metrics.json"
    )
    pose_comparison_path = (
        root / "reports/exp010_flashbind_crystal_pose/paired_pose_rmsd.csv"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    nesso_scatter_paths, nesso_rows = write_standardized_scatter(
        model="nesso1",
        manifest=manifest,
        prediction_path=nesso_prediction_path,
        figures=figures,
    )
    boltz_scatter_paths, boltz_rows = write_standardized_scatter(
        model="boltz2_msa1024",
        manifest=manifest,
        prediction_path=boltz_prediction_path,
        figures=figures,
    )
    flashbind_scatter_paths, flashbind_rows = write_standardized_scatter(
        model="flashbind",
        manifest=manifest,
        prediction_path=flashbind_prediction_path,
        figures=figures,
    )
    pose_paths, pose_rows = write_pose_comparison(
        source=pose_comparison_path,
        figures=figures,
    )

    runtime_path = figures / "nesso1_runtime.csv"
    runtime_rows = []
    for target in TARGET_ORDER:
        run_path = (
            root
            / f"runs/exp007_fepplus4_boltz2_nesso1/nesso1/seed42/{target}/run.json"
        )
        run = json.loads(run_path.read_text(encoding="utf-8"))
        n = sum(sample["target_id"] == target for sample in manifest["samples"])
        runtime_rows.append(
            {
                "target": target,
                "compounds": n,
                "wall_seconds": f"{run['wall_time_seconds']:.3f}",
                "seconds_per_compound": f"{run['wall_time_seconds'] / n:.3f}",
                "peak_gpu_delta_mib": run["peak_gpu_memory_above_baseline_mib"],
                "max_host_rss_kib": run["max_host_rss_kib"],
            }
        )
    with runtime_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=runtime_rows[0].keys())
        writer.writeheader()
        writer.writerows(runtime_rows)

    bibliography = root / "literature/references.bib"
    bibliography_snapshot = report / "references.bib"
    shutil.copy2(bibliography, bibliography_snapshot)

    source_paths = [
        manifest_path,
        nesso_prediction_path,
        boltz_prediction_path,
        flashbind_prediction_path,
        nesso_metrics_path,
        comparison_metrics_path,
        flashbind_metrics_path,
        pose_comparison_path,
        bibliography,
    ]
    output_paths = [
        *nesso_scatter_paths,
        *boltz_scatter_paths,
        *flashbind_scatter_paths,
        *pose_paths,
        runtime_path,
        bibliography_snapshot,
    ]
    records = {
        "sources": [
            {"path": str(path.relative_to(root)), "sha256": sha256(path)}
            for path in source_paths
        ],
        "outputs": [
            {"path": str(path.relative_to(root)), "sha256": sha256(path)}
            for path in output_paths
        ],
    }
    (figures / "assets.json").write_text(
        json.dumps(records, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Wrote {nesso_rows} Nesso, {boltz_rows} Boltz, and "
        f"{flashbind_rows} FlashBind scatter rows "
        f"plus {pose_rows} paired pose rows and {len(runtime_rows)} runtime rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
