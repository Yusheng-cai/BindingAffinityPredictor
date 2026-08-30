#!/usr/bin/env python3
"""Compare the Nesso Runs N' Poses discovery and confirmation cohorts."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from analyze_nesso_rnp_distograms import STRUCTURAL_METRICS, aggregate_results


LOWER_IS_BETTER = {
    "interface_distogram_nll",
    "interface_expected_distance_mae_angstrom",
}


def boundary_comparison(
    records: list[dict[str, Any]], metric: str, rng: np.random.Generator
) -> dict[str, Any]:
    relevant = [
        record
        for record in records
        if record["familiarity_bin"] in {"40-50", "50-60"}
    ]
    below = np.asarray(
        [r["metrics"][metric] for r in relevant if r["familiarity_bin"] == "40-50"],
        dtype=float,
    )
    above = np.asarray(
        [r["metrics"][metric] for r in relevant if r["familiarity_bin"] == "50-60"],
        dtype=float,
    )
    clusters = np.asarray(
        [record["runs_n_poses_cluster"] for record in relevant], dtype=object
    )
    unique_clusters = np.unique(clusters)
    differences = []
    for _ in range(2000):
        sampled_clusters = rng.choice(
            unique_clusters, size=len(unique_clusters), replace=True
        )
        sampled = [
            record
            for cluster in sampled_clusters
            for record in relevant
            if record["runs_n_poses_cluster"] == cluster
        ]
        below_sample = [
            r["metrics"][metric] for r in sampled if r["familiarity_bin"] == "40-50"
        ]
        above_sample = [
            r["metrics"][metric] for r in sampled if r["familiarity_bin"] == "50-60"
        ]
        if below_sample and above_sample:
            differences.append(np.median(above_sample) - np.median(below_sample))
    differences = np.asarray(differences, dtype=float)
    direction = -1.0 if metric in LOWER_IS_BETTER else 1.0
    return {
        "metric": metric,
        "difference_definition": "median(50-60) - median(40-50)",
        "direction_for_improvement": "negative" if metric in LOWER_IS_BETTER else "positive",
        "n_40_50": int(len(below)),
        "n_50_60": int(len(above)),
        "median_40_50": float(np.median(below)),
        "median_50_60": float(np.median(above)),
        "observed_difference": float(np.median(above) - np.median(below)),
        "bootstrap_lower_95": float(np.quantile(differences, 0.025)),
        "bootstrap_upper_95": float(np.quantile(differences, 0.975)),
        "bootstrap_iterations_finite": int(len(differences)),
        "bootstrap_probability_of_improvement": float(np.mean(direction * differences > 0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery-analysis", type=Path, required=True)
    parser.add_argument("--confirmation-analysis", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    discovery = json.loads(args.discovery_analysis.read_text())
    confirmation = json.loads(args.confirmation_analysis.read_text())
    cohorts = {
        "discovery_100": discovery["records"],
        "confirmation_50": confirmation["records"],
    }
    discovery_ids = {record["sample_id"] for record in cohorts["discovery_100"]}
    confirmation_ids = {record["sample_id"] for record in cohorts["confirmation_50"]}
    overlap = sorted(discovery_ids & confirmation_ids)
    if overlap:
        raise ValueError(f"cohorts overlap: {overlap}")

    combined = cohorts["discovery_100"] + cohorts["confirmation_50"]
    overlap_counts = {
        "sample_id": len(overlap),
        "runs_n_poses_cluster": len(
            {record["runs_n_poses_cluster"] for record in cohorts["discovery_100"]}
            & {
                record["runs_n_poses_cluster"]
                for record in cohorts["confirmation_50"]
            }
        ),
    }
    rng = np.random.default_rng(2026082051)
    boundary_metrics = (
        "interface_distogram_nll",
        "interface_expected_distance_mae_angstrom",
        "token_contact_average_precision_at_6a",
        "physical_pocket_f1_at_6a",
    )
    summary = {
        "schema_version": 1,
        "experiment_id": "exp013_nesso1_rnp_confirmation50",
        "cohort_overlap_count": len(overlap),
        "cross_cohort_overlap_counts": overlap_counts,
        "cohort_sizes": {name: len(records) for name, records in cohorts.items()},
        "aggregates": {
            "discovery_100": discovery["aggregate"],
            "confirmation_50": confirmation["aggregate"],
            "combined_150": aggregate_results(
                combined, [], iterations=2000, seed=2026082051
            ),
        },
        "boundary_40_50_vs_50_60": {
            name: {
                metric: boundary_comparison(records, metric, rng)
                for metric in boundary_metrics
            }
            for name, records in {**cohorts, "combined_150": combined}.items()
        },
        "interpretation_policy": {
            "familiarity_trend": "assess Spearman direction and interval separately in discovery and confirmation cohorts",
            "score_50_boundary": "exploratory; call replicated only when metric direction agrees in both cohorts and uncertainty is reported",
        },
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "comparison.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )

    fields = [
        "cohort",
        "sample_id",
        "familiarity_score_0_to_100",
        "familiarity_bin",
        *STRUCTURAL_METRICS,
        "entropy_crop_pl",
    ]
    with (args.output_dir / "combined_per_system_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for cohort, records in cohorts.items():
            for record in records:
                row = {
                    "cohort": cohort,
                    "sample_id": record["sample_id"],
                    "familiarity_score_0_to_100": record[
                        "familiarity_score_0_to_100"
                    ],
                    "familiarity_bin": record["familiarity_bin"],
                    "entropy_crop_pl": record["native_nesso_outputs"].get(
                        "entropy_crop_pl"
                    ),
                }
                row.update(record["metrics"])
                writer.writerow(row)
    print(f"Wrote discovery/confirmation comparison to {args.output_dir}")


if __name__ == "__main__":
    main()
