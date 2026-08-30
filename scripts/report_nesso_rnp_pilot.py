#!/usr/bin/env python3
"""Promote compact Nesso Runs N' Poses results and make review figures."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 13,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 13,
    }
)


BIN_ORDER = ("0-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100")


def _bin_series(aggregate: dict, metric: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    median = []
    lower = []
    upper = []
    for bin_name in BIN_ORDER:
        result = aggregate["by_familiarity_bin"][bin_name]["metrics"][metric]
        median.append(result["median"])
        lower.append(result["lower_95"])
        upper.append(result["upper_95"])
    values = np.asarray(median, dtype=float)
    return values, values - np.asarray(lower), np.asarray(upper) - values


def plot_familiarity(aggregate: dict, output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11.2, 4.3), constrained_layout=True)
    x = np.arange(len(BIN_ORDER))
    panels = (
        ("interface_distogram_nll", "Interface distogram NLL", "lower is better", "#4267B2"),
        ("physical_pocket_f1_at_6a", "6 Å physical-pocket F1", "higher is better", "#2B8C68"),
    )
    for axis, (metric, title, direction, color) in zip(axes, panels):
        values, low_error, high_error = _bin_series(aggregate, metric)
        axis.bar(x, values, color=color, alpha=0.88, width=0.74)
        axis.errorbar(
            x,
            values,
            yerr=np.vstack((low_error, high_error)),
            fmt="none",
            ecolor="#202020",
            elinewidth=1.1,
            capsize=3,
        )
        axis.set_title(title, loc="left", fontweight="bold")
        axis.text(0.99, 0.98, direction, transform=axis.transAxes, ha="right", va="top", fontsize=9)
        axis.set_xticks(x, BIN_ORDER, rotation=35, ha="right")
        axis.set_xlabel("Runs N' Poses familiarity score bin")
        axis.grid(axis="y", alpha=0.22)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Median across systems (95% bootstrap CI)")
    figure.savefig(output, dpi=220)
    plt.close(figure)


def plot_entropy(records: list[dict], aggregate: dict, output: Path) -> None:
    entropy = np.asarray([r["native_nesso_outputs"]["entropy_crop_pl"] for r in records])
    error = np.asarray([r["metrics"]["interface_expected_distance_mae_angstrom"] for r in records])
    familiarity = np.asarray([r["familiarity_score_0_to_100"] for r in records])
    rho = aggregate["entropy_error_spearman"]["interface_expected_distance_mae_angstrom"]

    figure, axis = plt.subplots(figsize=(6.3, 4.7), constrained_layout=True)
    points = axis.scatter(
        entropy,
        error,
        c=familiarity,
        cmap="viridis",
        s=38,
        alpha=0.82,
        edgecolors="white",
        linewidths=0.35,
    )
    axis.set_xlabel("Nesso native protein–ligand crop entropy")
    axis.set_ylabel("Interface expected-distance MAE (Å)")
    axis.set_title("Nesso entropy detects structural error", loc="left", fontweight="bold")
    axis.text(
        0.02,
        0.96,
        f"Spearman $\\rho$ = {rho['rho']:.3f}\n"
        f"95% CI [{rho['lower_95']:.3f}, {rho['upper_95']:.3f}]",
        transform=axis.transAxes,
        va="top",
        bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.9},
    )
    axis.grid(alpha=0.22)
    axis.spines[["top", "right"]].set_visible(False)
    colorbar = figure.colorbar(points, ax=axis)
    colorbar.set_label("Structural familiarity score")
    figure.savefig(output, dpi=220)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-json", type=Path, required=True)
    parser.add_argument("--analysis-csv", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    analysis = json.loads(args.analysis_json.read_text())
    run = json.loads(args.run_json.read_text())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.analysis_csv, args.output_dir / "per_system_metrics.csv")

    compact = {
        "schema_version": 1,
        "experiment_id": analysis["experiment_id"],
        "run": {
            "status": run["status"],
            "return_code": run["return_code"],
            "command": run["command"],
            "project_git": run["project_git"],
            "started_utc": run["started_utc"],
            "ended_utc": run["ended_utc"],
            "wall_time_seconds": run["wall_time_seconds"],
            "max_host_rss_kib": run["max_host_rss_kib"],
            "peak_total_gpu_memory_mib": run["peak_total_gpu_memory_mib"],
            "peak_gpu_memory_above_baseline_mib": run["peak_gpu_memory_above_baseline_mib"],
            "gpu_before": run["gpu_before"],
        },
        "metric_definitions": analysis["metric_definitions"],
        "failures": analysis["failures"],
        "aggregate": analysis["aggregate"],
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(compact, indent=2, allow_nan=False) + "\n"
    )
    plot_familiarity(analysis["aggregate"], args.output_dir / "familiarity_metrics.png")
    plot_entropy(
        analysis["records"], analysis["aggregate"], args.output_dir / "entropy_error.png"
    )
    print(f"Wrote compact results and figures to {args.output_dir}")


if __name__ == "__main__":
    main()
