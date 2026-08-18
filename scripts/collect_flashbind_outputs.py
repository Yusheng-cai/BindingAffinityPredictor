#!/usr/bin/env python3
"""Convert a FlashBind ensemble JSON file to the canonical affinity CSV schema."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from affinity_benchmark.adapters.affinity_models import (
    flashbind_record_id,
    load_flashbind_predictions,
)
from affinity_benchmark.data.manifest import load_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--predictions-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    predictions = load_flashbind_predictions(args.predictions_json)
    rows = []
    expected_record_ids = set()
    for sample in manifest["samples"]:
        record_id = flashbind_record_id(sample)
        expected_record_ids.add(record_id)
        native = predictions.get(record_id)
        row = {
            "sample_id": sample["sample_id"],
            "target_id": sample["target_id"],
            "model": "flashbind",
            "flashbind_record_id": record_id,
            "native_output_path": str(args.predictions_json),
        }
        if native is None:
            row["status"] = "missing"
        elif native["status"] != "success":
            row["status"] = "failed"
            row["failure_reason"] = native.get("failure_reason", "upstream_failed")
        else:
            row.update(
                {
                    "status": "complete",
                    "affinity_pred_value": native["pred_value"],
                    "pred_value": native["pred_value"],
                    "pred_value_raw": native["pred_value_raw"],
                    "mw": native["mw"],
                    "n_models": native.get("n_models"),
                }
            )
        rows.append(row)

    unexpected = sorted(set(predictions) - expected_record_ids)
    if unexpected:
        raise ValueError(f"unexpected FlashBind records: {unexpected}")

    fieldnames = sorted({key for row in rows for key in row})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "model": "flashbind",
        "complete": sum(row["status"] == "complete" for row in rows),
        "failed": sum(row["status"] == "failed" for row in rows),
        "missing": sum(row["status"] == "missing" for row in rows),
        "expected": len(rows),
        "output": str(args.output),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
