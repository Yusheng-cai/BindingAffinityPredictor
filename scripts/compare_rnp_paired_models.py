#!/usr/bin/env python3
"""Create the paired Nesso-1/Boltz-2 Runs N' Poses comparison."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


METRICS = {
    "interface_expected_distance_mae_angstrom": {
        "label": "Interface distance MAE (Å)",
        "higher_is_better": False,
    },
    "token_contact_f1_at_6a": {
        "label": "6 Å token-pair contact F1",
        "higher_is_better": True,
    },
    "physical_pocket_f1_at_6a": {
        "label": "6 Å residue-pocket F1",
        "higher_is_better": True,
    },
}

FAMILIARITY_BINS = ((0, 20), (20, 40), (40, 60), (60, 80), (80, 100))


def percentile_interval(values: np.ndarray) -> dict[str, float]:
    return {
        "lower_95": float(np.quantile(values, 0.025)),
        "upper_95": float(np.quantile(values, 0.975)),
    }


def paired_comparison(
    nesso: np.ndarray,
    boltz: np.ndarray,
    *,
    higher_is_better: bool,
    rng: np.random.Generator,
    iterations: int,
) -> dict[str, Any]:
    if higher_is_better:
        improvement = boltz - nesso
    else:
        improvement = nesso - boltz
    ties = np.isclose(improvement, 0.0, atol=1e-12, rtol=1e-9)
    boot_medians = np.empty(iterations, dtype=float)
    boot_means = np.empty(iterations, dtype=float)
    for index in range(iterations):
        chosen = rng.integers(0, len(improvement), size=len(improvement))
        boot_medians[index] = np.median(improvement[chosen])
        boot_means[index] = np.mean(improvement[chosen])
    return {
        "n": len(improvement),
        "nesso1_mean": float(np.mean(nesso)),
        "nesso1_median": float(np.median(nesso)),
        "boltz2_mean": float(np.mean(boltz)),
        "boltz2_median": float(np.median(boltz)),
        "improvement_definition": (
            "Boltz-2 minus Nesso-1" if higher_is_better else "Nesso-1 minus Boltz-2"
        ),
        "positive_improvement_favors": "Boltz-2",
        "paired_median_improvement": float(np.median(improvement)),
        "paired_median_improvement_interval": percentile_interval(boot_medians),
        "paired_mean_improvement": float(np.mean(improvement)),
        "paired_mean_improvement_interval": percentile_interval(boot_means),
        "boltz2_better": int(np.sum(improvement > 0)),
        "nesso1_better": int(np.sum(improvement < 0)),
        "ties": int(np.sum(ties)),
    }


def binned_familiarity_summary(
    records: list[dict[str, Any]],
    metric: str,
    *,
    rng: np.random.Generator,
    iterations: int,
) -> list[dict[str, Any]]:
    """Summarize a metric in fixed familiarity bins with bootstrap median CIs."""

    summaries = []
    for lower, upper in FAMILIARITY_BINS:
        values = np.asarray(
            [
                record["metrics"][metric]
                for record in records
                if lower <= record["familiarity_score_0_to_100"]
                and (
                    record["familiarity_score_0_to_100"] < upper
                    or (
                        upper == FAMILIARITY_BINS[-1][1]
                        and record["familiarity_score_0_to_100"] <= upper
                    )
                )
            ],
            dtype=float,
        )
        if len(values) == 0:
            summaries.append(
                {
                    "bin": f"{lower}-{upper}",
                    "bin_midpoint": (lower + upper) / 2,
                    "n": 0,
                    "median": None,
                    "lower_95": None,
                    "upper_95": None,
                }
            )
            continue
        boot_medians = np.empty(iterations, dtype=float)
        for index in range(iterations):
            chosen = rng.integers(0, len(values), size=len(values))
            boot_medians[index] = np.median(values[chosen])
        summaries.append(
            {
                "bin": f"{lower}-{upper}",
                "bin_midpoint": (lower + upper) / 2,
                "n": len(values),
                "median": float(np.median(values)),
                **percentile_interval(boot_medians),
            }
        )
    return summaries


def plot_paired(
    rows: list[dict[str, Any]], output_pdf: Path, output_png: Path
) -> None:
    style = {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.titlesize": 12,
        "font.family": "sans-serif",
    }
    with plt.rc_context(style):
        figure, axes = plt.subplots(
            1, 3, figsize=(10.8, 3.45), constrained_layout=True
        )
        familiarity = np.asarray(
            [row["boltz2_june2023_familiarity"] for row in rows]
        )
        scatter = None
        for axis, (metric, definition) in zip(axes, METRICS.items()):
            nesso = np.asarray([row[f"nesso1_{metric}"] for row in rows])
            boltz = np.asarray([row[f"boltz2_{metric}"] for row in rows])
            scatter = axis.scatter(
                nesso,
                boltz,
                c=familiarity,
                cmap="viridis",
                vmin=0,
                vmax=100,
                s=30,
                edgecolor="white",
                linewidth=0.35,
                alpha=0.9,
            )
            low = min(float(np.min(nesso)), float(np.min(boltz)))
            high = max(float(np.max(nesso)), float(np.max(boltz)))
            padding = max(0.04 * (high - low), 0.02)
            axis.plot(
                [low - padding, high + padding],
                [low - padding, high + padding],
                "--",
                color="#777777",
                linewidth=1,
            )
            axis.set_xlim(low - padding, high + padding)
            axis.set_ylim(low - padding, high + padding)
            axis.set_xlabel("Nesso-1")
            axis.set_ylabel("Boltz-2")
            axis.grid(alpha=0.2)
            favored = (
                "below diagonal favors Boltz-2"
                if not definition["higher_is_better"]
                else "above diagonal favors Boltz-2"
            )
            axis.set_title(f"{definition['label']}\n{favored}")
        colorbar = figure.colorbar(scatter, ax=axes, shrink=0.82, pad=0.02)
        colorbar.set_label("June 2023 familiarity")
        figure.suptitle(
            "Paired structure recovery on 49 shared Runs N’ Poses systems"
        )
        output_pdf.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_pdf, bbox_inches="tight")
        figure.savefig(output_png, dpi=180, bbox_inches="tight")
        plt.close(figure)


def plot_familiarity_trends(
    model_records: dict[str, list[dict[str, Any]]],
    familiarity_results: dict[str, Any],
    binned_results: dict[str, dict[str, list[dict[str, Any]]]],
    output_pdf: Path,
    output_png: Path,
) -> None:
    """Plot each model against familiarity defined at its own cutoff date."""

    model_definitions = (
        ("nesso1", "Nesso-1", "30 Sep 2021", "#2B8C68"),
        ("boltz2", "Boltz-2", "1 Jun 2023", "#4267B2"),
    )
    style = {
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.size": 3,
        "ytick.major.width": 0.8,
        "legend.fontsize": 8.5,
        "figure.titlesize": 12,
        "font.family": "sans-serif",
        "text.usetex": False,
    }
    with plt.rc_context(style):
        figure, axes = plt.subplots(
            len(METRICS),
            len(model_definitions),
            figsize=(9.6, 9.0),
            sharex=True,
            constrained_layout=True,
        )
        for row_index, (metric, definition) in enumerate(METRICS.items()):
            all_metric_values = np.asarray(
                [
                    record["metrics"][metric]
                    for model_key, *_ in model_definitions
                    for record in model_records[model_key]
                ],
                dtype=float,
            )
            if definition["higher_is_better"]:
                shared_limits = (-0.04, 1.04)
            else:
                span = float(np.max(all_metric_values) - np.min(all_metric_values))
                shared_limits = (
                    max(0.0, float(np.min(all_metric_values)) - 0.04 * span),
                    float(np.max(all_metric_values)) + 0.06 * span,
                )
            for column_index, (model_key, model_label, cutoff, color) in enumerate(
                model_definitions
            ):
                axis = axes[row_index, column_index]
                records = model_records[model_key]
                x = np.asarray(
                    [record["familiarity_score_0_to_100"] for record in records],
                    dtype=float,
                )
                y = np.asarray(
                    [record["metrics"][metric] for record in records], dtype=float
                )
                axis.scatter(
                    x,
                    y,
                    s=22,
                    color=color,
                    alpha=0.38,
                    linewidth=0,
                    label="individual system",
                )
                bins = binned_results[model_key][metric]
                centers = np.asarray([entry["bin_midpoint"] for entry in bins])
                medians = np.asarray([entry["median"] for entry in bins], dtype=float)
                lower = np.asarray([entry["lower_95"] for entry in bins], dtype=float)
                upper = np.asarray([entry["upper_95"] for entry in bins], dtype=float)
                axis.errorbar(
                    centers,
                    medians,
                    yerr=np.vstack((medians - lower, upper - medians)),
                    fmt="o-",
                    color=color,
                    markeredgecolor="white",
                    markeredgewidth=0.7,
                    markersize=5.5,
                    linewidth=1.8,
                    capsize=2.5,
                    label="bin median with 95% CI",
                )
                for entry in bins:
                    axis.text(
                        entry["bin_midpoint"],
                        0.025,
                        f"n={entry['n']}",
                        transform=axis.get_xaxis_transform(),
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        color="#555555",
                    )
                rho = familiarity_results[model_key][metric]["rho"]
                axis.set_title(
                    f"{model_label}: structures before {cutoff}\n"
                    f"Spearman $\\rho$={rho:+.3f}",
                    loc="left",
                )
                axis.set_xlim(-3, 103)
                axis.set_ylim(*shared_limits)
                axis.grid(alpha=0.18)
                axis.spines[["top", "right"]].set_visible(False)
                if column_index == 0:
                    axis.set_ylabel(definition["label"])
                if row_index == len(METRICS) - 1:
                    axis.set_xlabel("Runs N’ Poses familiarity score")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        figure.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.05),
            ncol=2,
            frameon=False,
        )
        figure.suptitle(
            "Both models improve as test complexes become more structurally familiar",
            fontweight="bold",
        )
        output_pdf.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_pdf, bbox_inches="tight")
        figure.savefig(output_png, dpi=180, bbox_inches="tight")
        plt.close(figure)


def _familiarity_plot_style() -> dict[str, Any]:
    return {
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.size": 3,
        "ytick.major.width": 0.8,
        "legend.fontsize": 8.5,
        "figure.titlesize": 12,
        "font.family": "sans-serif",
        "text.usetex": False,
    }


def _shared_metric_limits(
    model_records: dict[str, list[dict[str, Any]]],
    metric: str,
    *,
    higher_is_better: bool,
) -> tuple[float, float]:
    if higher_is_better:
        return (-0.04, 1.04)
    values = np.asarray(
        [
            record["metrics"][metric]
            for records in model_records.values()
            for record in records
        ],
        dtype=float,
    )
    span = float(np.max(values) - np.min(values))
    return (
        max(0.0, float(np.min(values)) - 0.04 * span),
        float(np.max(values)) + 0.06 * span,
    )


def plot_familiarity_bars(
    model_records: dict[str, list[dict[str, Any]]],
    binned_results: dict[str, dict[str, list[dict[str, Any]]]],
    output_pdf: Path,
    output_png: Path,
) -> None:
    """Plot bin medians and bootstrap intervals as presentation-style bars."""

    model_definitions = (
        ("nesso1", "Nesso-1", "PDB structural-training cutoff: 30 Sep 2021", "#2B8C68"),
        ("boltz2", "Boltz-2", "PDB structural-training cutoff: 1 Jun 2023", "#4267B2"),
    )
    bin_labels = [f"{lower}–{upper}" for lower, upper in FAMILIARITY_BINS]
    positions = np.arange(len(FAMILIARITY_BINS))
    with plt.rc_context(_familiarity_plot_style()):
        figure, axes = plt.subplots(
            len(METRICS),
            len(model_definitions),
            figsize=(9.6, 8.4),
            sharex=True,
            constrained_layout=True,
        )
        for row_index, (metric, definition) in enumerate(METRICS.items()):
            limits = _shared_metric_limits(
                model_records,
                metric,
                higher_is_better=definition["higher_is_better"],
            )
            for column_index, (model_key, model_label, cutoff, color) in enumerate(
                model_definitions
            ):
                axis = axes[row_index, column_index]
                bins = binned_results[model_key][metric]
                medians = np.asarray([entry["median"] for entry in bins], dtype=float)
                lower = np.asarray([entry["lower_95"] for entry in bins], dtype=float)
                upper = np.asarray([entry["upper_95"] for entry in bins], dtype=float)
                bars = axis.bar(
                    positions,
                    medians,
                    color=color,
                    alpha=0.86,
                    width=0.72,
                    yerr=np.vstack((medians - lower, upper - medians)),
                    error_kw={"ecolor": "#333333", "elinewidth": 1, "capsize": 2.5},
                )
                for bar, entry in zip(bars, bins):
                    axis.text(
                        bar.get_x() + bar.get_width() / 2,
                        0.025,
                        f"n={entry['n']}",
                        transform=axis.get_xaxis_transform(),
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        color="#555555",
                    )
                axis.set_ylim(*limits)
                axis.grid(axis="y", alpha=0.18)
                axis.spines[["top", "right"]].set_visible(False)
                if row_index == 0:
                    axis.set_title(f"{model_label}\n{cutoff}", loc="left")
                if column_index == 0:
                    axis.set_ylabel(definition["label"])
                if row_index == len(METRICS) - 1:
                    axis.set_xticks(positions, bin_labels)
                    axis.set_xlabel("Familiarity interval")
        figure.suptitle(
            "Median structural accuracy by model-specific familiarity",
            fontweight="bold",
        )
        output_pdf.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_pdf, bbox_inches="tight")
        figure.savefig(output_png, dpi=180, bbox_inches="tight")
        plt.close(figure)


def plot_familiarity_violins(
    model_records: dict[str, list[dict[str, Any]]],
    output_pdf: Path,
    output_png: Path,
    *,
    jitter_seed: int = 2026082901,
) -> None:
    """Show each bin's full distribution with violins, boxes, and raw points."""

    model_definitions = (
        ("nesso1", "Nesso-1", "structures before 30 Sep 2021", "#2B8C68"),
        ("boltz2", "Boltz-2", "structures before 1 Jun 2023", "#4267B2"),
    )
    bin_labels = [f"{lower}–{upper}" for lower, upper in FAMILIARITY_BINS]
    positions = np.arange(1, len(FAMILIARITY_BINS) + 1)
    rng = np.random.default_rng(jitter_seed)
    with plt.rc_context(_familiarity_plot_style()):
        figure, axes = plt.subplots(
            len(METRICS),
            len(model_definitions),
            figsize=(9.6, 8.4),
            sharex=True,
            constrained_layout=True,
        )
        for row_index, (metric, definition) in enumerate(METRICS.items()):
            limits = _shared_metric_limits(
                model_records,
                metric,
                higher_is_better=definition["higher_is_better"],
            )
            for column_index, (model_key, model_label, cutoff, color) in enumerate(
                model_definitions
            ):
                axis = axes[row_index, column_index]
                records = model_records[model_key]
                groups = []
                for lower, upper in FAMILIARITY_BINS:
                    groups.append(
                        np.asarray(
                            [
                                record["metrics"][metric]
                                for record in records
                                if lower <= record["familiarity_score_0_to_100"]
                                and (
                                    record["familiarity_score_0_to_100"] < upper
                                    or (upper == 100 and record["familiarity_score_0_to_100"] <= 100)
                                )
                            ],
                            dtype=float,
                        )
                    )
                violins = axis.violinplot(
                    groups,
                    positions=positions,
                    widths=0.78,
                    showmeans=False,
                    showmedians=False,
                    showextrema=False,
                )
                for body in violins["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor(color)
                    body.set_alpha(0.20)
                box = axis.boxplot(
                    groups,
                    positions=positions,
                    widths=0.22,
                    showfliers=False,
                    patch_artist=True,
                    medianprops={"color": "white", "linewidth": 1.4},
                    boxprops={"facecolor": color, "edgecolor": color, "linewidth": 1},
                    whiskerprops={"color": color, "linewidth": 1},
                    capprops={"color": color, "linewidth": 1},
                )
                del box
                for position, values in zip(positions, groups):
                    jitter = rng.uniform(-0.13, 0.13, size=len(values))
                    axis.scatter(
                        position + jitter,
                        values,
                        s=14,
                        color=color,
                        alpha=0.58,
                        edgecolor="white",
                        linewidth=0.25,
                        zorder=3,
                    )
                    axis.text(
                        position,
                        0.025,
                        f"n={len(values)}",
                        transform=axis.get_xaxis_transform(),
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        color="#555555",
                    )
                axis.set_ylim(*limits)
                axis.grid(axis="y", alpha=0.18)
                axis.spines[["top", "right"]].set_visible(False)
                if row_index == 0:
                    axis.set_title(f"{model_label}\n{cutoff}", loc="left")
                if column_index == 0:
                    axis.set_ylabel(definition["label"])
                if row_index == len(METRICS) - 1:
                    axis.set_xticks(positions, bin_labels)
                    axis.set_xlabel("Familiarity interval")
        figure.suptitle(
            "Distribution of structural accuracy within familiarity intervals",
            fontweight="bold",
        )
        output_pdf.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_pdf, bbox_inches="tight")
        figure.savefig(output_png, dpi=180, bbox_inches="tight")
        plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--nesso-metrics", type=Path, required=True)
    parser.add_argument("--boltz-metrics", type=Path, required=True)
    parser.add_argument("--nesso-run", type=Path, required=True)
    parser.add_argument("--boltz-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=2026082803)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    samples = {sample["sample_id"]: sample for sample in manifest["samples"]}
    nesso_result = json.loads(args.nesso_metrics.read_text())
    boltz_result = json.loads(args.boltz_metrics.read_text())
    nesso = {record["sample_id"]: record for record in nesso_result["records"]}
    boltz = {record["sample_id"]: record for record in boltz_result["records"]}
    paired_ids = sorted(set(nesso) & set(boltz))
    rows = []
    for sample_id in paired_ids:
        sample = samples[sample_id]
        row: dict[str, Any] = {
            "sample_id": sample_id,
            "runs_n_poses_cluster": sample["runs_n_poses"]["cluster"],
            "boltz2_june2023_familiarity": sample["runs_n_poses"]["similarity_score_0_to_100"],
            "nesso1_sep2021_familiarity": sample["runs_n_poses"]["nesso1_sep2021_familiarity"]["similarity_score_0_to_100"],
        }
        for metric in METRICS:
            row[f"nesso1_{metric}"] = nesso[sample_id]["metrics"][metric]
            row[f"boltz2_{metric}"] = boltz[sample_id]["metrics"][metric]
        rows.append(row)

    rng = np.random.default_rng(args.bootstrap_seed)
    comparisons = {}
    for metric, definition in METRICS.items():
        comparisons[metric] = paired_comparison(
            np.asarray([row[f"nesso1_{metric}"] for row in rows], dtype=float),
            np.asarray([row[f"boltz2_{metric}"] for row in rows], dtype=float),
            higher_is_better=definition["higher_is_better"],
            rng=rng,
            iterations=args.bootstrap_iterations,
        )

    model_records = {
        "nesso1": nesso_result["records"],
        "boltz2": boltz_result["records"],
    }
    familiarity_results = {
        "nesso1": {
            metric: nesso_result["aggregate"]["spearman_vs_familiarity"][metric]
            for metric in METRICS
        },
        "boltz2": {
            metric: boltz_result["aggregate"]["spearman_vs_familiarity"][metric]
            for metric in METRICS
        },
    }
    binned_results: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for model_key, records in model_records.items():
        binned_results[model_key] = {}
        for metric in METRICS:
            binned_results[model_key][metric] = binned_familiarity_summary(
                records,
                metric,
                rng=rng,
                iterations=args.bootstrap_iterations,
            )

    nesso_run = json.loads(args.nesso_run.read_text())
    boltz_run = json.loads(args.boltz_run.read_text())
    ligand_rmsd = np.asarray(
        [
            record["metrics"][
                "ligand_heavy_atom_rmsd_after_global_protein_alignment_angstrom"
            ]
            for record in boltz.values()
        ],
        dtype=float,
    )
    protein_rmsd = np.asarray(
        [record["metrics"]["protein_global_ca_rmsd_angstrom"] for record in boltz.values()],
        dtype=float,
    )
    summary = {
        "schema_version": 1,
        "experiment_id": "exp015_boltz2_nesso1_rnp_postcutoff50",
        "paired_systems": len(rows),
        "completion": {
            "nesso1": {"attempted": 50, "completed": len(nesso), "failed": 50 - len(nesso)},
            "boltz2": {
                "attempted": 50,
                "completed": len(boltz),
                "failed": 50 - len(boltz),
                "failure_sample_ids": [failure["sample_id"] for failure in boltz_result["failures"]],
            },
        },
        "runtime": {
            "scope": "local input processing plus inference; excludes remote ColabFold MSA generation",
            "nesso1_wall_seconds": nesso_run["wall_time_seconds"],
            "boltz2_wall_seconds": boltz_run["wall_time_seconds"],
            "nesso1_seconds_per_attempted_system": nesso_run["wall_time_seconds"] / 50,
            "boltz2_seconds_per_attempted_system": boltz_run["wall_time_seconds"] / 50,
            "boltz2_to_nesso1_ratio_per_attempted_system": boltz_run["wall_time_seconds"] / nesso_run["wall_time_seconds"],
        },
        "paired_metrics": comparisons,
        "model_specific_familiarity": {
            "nesso1_cutoff": "2021-09-30",
            "nesso1": {
                metric: nesso_result["aggregate"]["spearman_vs_familiarity"][metric]
                for metric in METRICS
            },
            "boltz2_cutoff": "2023-06-01",
            "boltz2": {
                metric: boltz_result["aggregate"]["spearman_vs_familiarity"][metric]
                for metric in METRICS
            },
            "binned_medians": binned_results,
        },
        "boltz2_coordinate_metrics": {
            "n": len(ligand_rmsd),
            "protein_global_ca_rmsd_angstrom_median": float(np.median(protein_rmsd)),
            "ligand_rmsd_angstrom_median": float(np.median(ligand_rmsd)),
            "ligand_rmsd_below_2_angstrom": int(np.sum(ligand_rmsd < 2.0)),
            "ligand_rmsd_below_2_angstrom_fraction": float(np.mean(ligand_rmsd < 2.0)),
            "ligand_rmsd_below_5_angstrom": int(np.sum(ligand_rmsd < 5.0)),
            "ligand_rmsd_below_5_angstrom_fraction": float(np.mean(ligand_rmsd < 5.0)),
        },
        "bootstrap": {
            "iterations": args.bootstrap_iterations,
            "seed": args.bootstrap_seed,
            "unit": "paired system (all 49 systems have unique Runs N' Poses clusters)",
        },
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    with (args.output_dir / "paired_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    familiarity_rows = []
    for model_key, metrics in binned_results.items():
        cutoff = "2021-09-30" if model_key == "nesso1" else "2023-06-01"
        for metric, bins in metrics.items():
            for entry in bins:
                familiarity_rows.append(
                    {
                        "model": model_key,
                        "cutoff": cutoff,
                        "metric": metric,
                        **entry,
                    }
                )
    with (args.output_dir / "familiarity_binned_metrics.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(familiarity_rows[0]))
        writer.writeheader()
        writer.writerows(familiarity_rows)
    plot_paired(
        rows,
        args.output_dir / "paired_structure_metrics.pdf",
        args.output_dir / "paired_structure_metrics.png",
    )
    plot_familiarity_trends(
        model_records,
        familiarity_results,
        binned_results,
        args.output_dir / "familiarity_vs_accuracy.pdf",
        args.output_dir / "familiarity_vs_accuracy.png",
    )
    plot_familiarity_bars(
        model_records,
        binned_results,
        args.output_dir / "familiarity_bar_charts.pdf",
        args.output_dir / "familiarity_bar_charts.png",
    )
    plot_familiarity_violins(
        model_records,
        args.output_dir / "familiarity_violin_plots.pdf",
        args.output_dir / "familiarity_violin_plots.png",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
