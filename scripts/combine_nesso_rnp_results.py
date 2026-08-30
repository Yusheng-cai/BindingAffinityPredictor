#!/usr/bin/env python3
"""Combine frozen Nesso Runs N' Poses manifests and scored records."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from affinity_benchmark.data.manifest import validate_manifest
from analyze_nesso_rnp_distograms import _write_csv, aggregate_results


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, action="append", required=True)
    parser.add_argument("--analysis", type=Path, action="append", required=True)
    parser.add_argument("--run-json", type=Path, action="append", required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    parser.add_argument("--results-output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not (len(args.manifest) == len(args.analysis) == len(args.run_json)):
        raise ValueError("provide the same number of manifests, analyses, and run records")

    manifests = [json.loads(path.read_text()) for path in args.manifest]
    analyses = [json.loads(path.read_text()) for path in args.analysis]
    runs = [json.loads(path.read_text()) for path in args.run_json]
    for manifest in manifests:
        validate_manifest(manifest)

    samples = [sample for manifest in manifests for sample in manifest["samples"]]
    sample_ids = [sample["sample_id"] for sample in samples]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("component manifests contain duplicate sample IDs")
    records = [record for analysis in analyses for record in analysis["records"]]
    record_ids = [record["sample_id"] for record in records]
    if set(record_ids) != set(sample_ids) or len(record_ids) != len(set(record_ids)):
        raise ValueError("combined analysis records do not map one-to-one to manifest samples")
    failures = [failure for analysis in analyses for failure in analysis["failures"]]

    bin_counts = Counter(
        sample["runs_n_poses"]["similarity_bin"] for sample in samples
    )
    combined_manifest = {
        "schema_version": 1,
        "manifest_id": "rnp_nesso1_combined150_sep2021",
        "created_on": "2026-08-20",
        "purpose": "Combined 150-system Runs N' Poses analysis for Nesso-1 structural familiarity",
        "source_dataset": manifests[0]["source_dataset"],
        "selection": {
            "status": "posthoc_union_of_two_manifests_frozen_before_their_respective_inference_runs",
            "structural_training_cutoff": "2021-09-30",
            "familiarity_metric": "sucos_shape_pocket_qcov",
            "component_manifests": [
                {
                    "path": str(path),
                    "manifest_id": manifest["manifest_id"],
                    "sha256": sha256(path),
                    "samples": len(manifest["samples"]),
                }
                for path, manifest in zip(args.manifest, manifests)
            ],
            "total_samples": len(samples),
            "bin_counts": dict(sorted(bin_counts.items())),
            "unique_sample_ids": len(set(sample_ids)),
            "unique_pdb_entries": len(
                {sample["structure_reference"]["pdb_id"] for sample in samples}
            ),
            "unique_runs_n_poses_clusters": len(
                {sample["runs_n_poses"]["cluster"] for sample in samples}
            ),
            "unique_canonical_ligand_smiles": len(
                {sample["ligand"]["input_smiles"] for sample in samples}
            ),
            "combined_uncertainty_resampling_unit": "Runs N' Poses cluster",
        },
        "samples": samples,
    }
    validate_manifest(combined_manifest)
    args.manifest_output.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_output.write_text(
        json.dumps(combined_manifest, indent=2) + "\n", encoding="utf-8"
    )

    aggregate = aggregate_results(records, failures, iterations=2000, seed=2026082051)
    combined_analysis = {
        "schema_version": 1,
        "experiment_id": "exp014_nesso1_rnp_combined150",
        "manifest": str(args.manifest_output),
        "component_analyses": [str(path) for path in args.analysis],
        "component_runs": runs,
        "metric_definitions": analyses[0]["metric_definitions"],
        "records": records,
        "failures": failures,
        "aggregate": aggregate,
    }
    args.results_output_dir.mkdir(parents=True, exist_ok=True)
    (args.results_output_dir / "analysis.json").write_text(
        json.dumps(combined_analysis, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    compact = {key: value for key, value in combined_analysis.items() if key != "records"}
    (args.results_output_dir / "summary.json").write_text(
        json.dumps(compact, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    _write_csv(records, args.results_output_dir / "per_system_metrics.csv")
    print(
        f"Combined {len(samples)} manifest samples and {len(records)} scored records "
        f"with {len(failures)} failures"
    )


if __name__ == "__main__":
    main()
