#!/usr/bin/env python3
"""Build compact figures for the Week 35 paired Nesso-1/Boltz-2 report."""

from __future__ import annotations

import json
import csv
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
SUMMARY_PATH = (
    REPOSITORY_ROOT
    / "reports/exp015_boltz2_nesso1_rnp_postcutoff50/results/summary.json"
)
REVIEWED_FAMILIARITY_FIGURE = (
    REPOSITORY_ROOT
    / "reports/exp015_boltz2_nesso1_rnp_postcutoff50/results/familiarity_bar_charts.pdf"
)
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "weeks/2026-W35/report/figures"


def write_familiarity_figure() -> Path:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIRECTORY / "nesso_boltz_familiarity.pdf"
    shutil.copyfile(REVIEWED_FAMILIARITY_FIGURE, output)
    print(f"Wrote {output}")
    return output


def write_runtime_comparison() -> tuple[Path, Path]:
    """Plot observed wall times from the paired 50-system structure run."""

    summary = json.loads(SUMMARY_PATH.read_text())
    runtime = summary["runtime"]
    rows = [
        {
            "model": "Nesso-1",
            "attempted": summary["completion"]["nesso1"]["attempted"],
            "completed": summary["completion"]["nesso1"]["completed"],
            "wall_seconds": runtime["nesso1_wall_seconds"],
            "seconds_per_attempted_system": runtime[
                "nesso1_seconds_per_attempted_system"
            ],
        },
        {
            "model": "Boltz-2 MSA-1024",
            "attempted": summary["completion"]["boltz2"]["attempted"],
            "completed": summary["completion"]["boltz2"]["completed"],
            "wall_seconds": runtime["boltz2_wall_seconds"],
            "seconds_per_attempted_system": runtime[
                "boltz2_seconds_per_attempted_system"
            ],
        },
    ]

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIRECTORY / "nesso_boltz_runtime.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    style = {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 9,
        "ytick.labelsize": 8,
        "xtick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.size": 3,
        "ytick.major.width": 0.8,
        "font.family": "sans-serif",
    }
    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=(6.8, 4.2), constrained_layout=True)
        bars = ax.bar(
            [row["model"] for row in rows],
            [row["seconds_per_attempted_system"] for row in rows],
            color=["#2B8C68", "#4267B2"],
            alpha=0.9,
            width=0.62,
        )
        ax.bar_label(bars, fmt="%.1f s", padding=3, fontsize=9)
        speedup = runtime["boltz2_to_nesso1_ratio_per_attempted_system"]
        ax.set_title(
            f"Same RTX 3080: Nesso-1 was {speedup:.1f}× faster",
            loc="left",
            weight="bold",
        )
        ax.set_ylabel("Seconds per attempted system")
        ax.grid(axis="y", alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylim(
            0,
            max(row["seconds_per_attempted_system"] for row in rows) * 1.22,
        )

        figure_path = OUTPUT_DIRECTORY / "nesso_boltz_runtime.pdf"
        fig.savefig(
            figure_path,
            bbox_inches="tight",
            metadata={
                "Creator": "BindingAffinityPredictor Week 35 build_assets.py",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        plt.close(fig)
    print(f"Wrote {csv_path}")
    print(f"Wrote {figure_path}")
    return csv_path, figure_path


def main() -> None:
    write_familiarity_figure()
    write_runtime_comparison()


if __name__ == "__main__":
    main()
