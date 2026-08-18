#!/usr/bin/env python3
"""Analyze Boltz-2 and Nesso-1 predictions on the frozen FEP+4 benchmark."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from affinity_benchmark.data.manifest import load_manifest
from affinity_benchmark.metrics.affinity_regression import (
    assay_metrics,
    bootstrap_weighted_metric,
    compound_weighted_average,
    log10_micromolar_to_kcal_per_mol,
    paired_bootstrap_weighted_metric_difference,
)


PAPER_REFERENCE = {
    "boltz2": {"weighted_pearson_r": 0.66, "weighted_kendall_tau": 0.48},
    "nesso1": {"weighted_pearson_r": 0.80},
    "flashbind": {"weighted_pearson_r": 0.53, "weighted_kendall_tau": 0.38},
}


def read_predictions(path: Path) -> dict[str, dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["sample_id"]: row for row in csv.DictReader(handle)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260817)
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    samples = {sample["sample_id"]: sample for sample in manifest["samples"]}
    results = {
        "manifest_id": manifest["manifest_id"],
        "aggregation": "per-target metrics; compound-count-weighted target average",
        "temperature_k": 298.15,
        "bootstrap": {"iterations": args.bootstrap_iterations, "seed": args.bootstrap_seed},
        "models": {},
        "paper_reference": PAPER_REFERENCE,
    }

    observations_by_model = {}
    for prediction_path in args.predictions:
        predictions = read_predictions(prediction_path)
        completed = [row for row in predictions.values() if row["status"] == "complete"]
        if not completed:
            raise ValueError(f"no complete predictions in {prediction_path}")
        model = completed[0]["model"]
        per_target = {}
        observations = {}
        for target in ("cdk2", "tyk2", "jnk1", "p38"):
            target_samples = [sample for sample in samples.values() if sample["target_id"] == target]
            if not all(sample["sample_id"] in predictions and predictions[sample["sample_id"]]["status"] == "complete" for sample in target_samples):
                continue
            observed_log = np.asarray([sample["measurement"]["log10_value_uM"] for sample in target_samples])
            predicted_log = np.asarray([float(predictions[sample["sample_id"]]["affinity_pred_value"]) for sample in target_samples])
            observed = log10_micromolar_to_kcal_per_mol(observed_log)
            predicted = log10_micromolar_to_kcal_per_mol(predicted_log)
            per_target[target] = assay_metrics(observed, predicted)
            observations[target] = (observed, predicted)

        weighted = {}
        confidence_intervals = {}
        for metric in ("pearson_r", "spearman_rho", "kendall_tau", "pairwise_mae_kcal_mol", "mae_kcal_mol", "centered_mae_kcal_mol"):
            weighted[metric] = compound_weighted_average(per_target, metric)
            confidence_intervals[metric] = bootstrap_weighted_metric(
                observations,
                metric,
                iterations=args.bootstrap_iterations,
                seed=args.bootstrap_seed,
            )
        results["models"][model] = {
            "complete_samples": len(completed),
            "targets_analyzed": sorted(per_target),
            "per_target": per_target,
            "compound_weighted": weighted,
            "bootstrap_95_percent": confidence_intervals,
        }
        observations_by_model[model] = observations

    if {"nesso1", "boltz2"}.issubset(observations_by_model):
        common_targets = sorted(
            set(observations_by_model["nesso1"])
            & set(observations_by_model["boltz2"])
        )
        paired = {
            target: (
                observations_by_model["nesso1"][target][0],
                observations_by_model["nesso1"][target][1],
                observations_by_model["boltz2"][target][1],
            )
            for target in common_targets
        }
        metric_names = (
            "pearson_r",
            "spearman_rho",
            "kendall_tau",
            "pairwise_mae_kcal_mol",
            "mae_kcal_mol",
            "centered_mae_kcal_mol",
        )
        results["paired_comparison"] = {
            "difference": "nesso1_minus_boltz2",
            "targets": common_targets,
            "point_estimate": {
                metric: results["models"]["nesso1"]["compound_weighted"][metric]
                - results["models"]["boltz2"]["compound_weighted"][metric]
                for metric in metric_names
            },
            "bootstrap_95_percent": {
                metric: paired_bootstrap_weighted_metric_difference(
                    paired,
                    metric,
                    iterations=args.bootstrap_iterations,
                    seed=args.bootstrap_seed,
                )
                for metric in metric_names
            },
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
