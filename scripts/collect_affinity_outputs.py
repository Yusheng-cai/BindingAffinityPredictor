#!/usr/bin/env python3
"""Collect model-native affinity JSON files into one canonical CSV table."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from affinity_benchmark.adapters.affinity_models import load_native_affinity
from affinity_benchmark.data.manifest import load_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--model", choices=("boltz2", "nesso1"), required=True)
    parser.add_argument(
        "--raw-root",
        type=Path,
        action="append",
        required=True,
        help="Root containing native outputs; repeat for disjoint target runs.",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    samples = {sample["sample_id"]: sample for sample in manifest["samples"]}
    found: dict[str, Path] = {}
    pattern = "affinity_*.json" if args.model == "boltz2" else "affinity.json"
    for raw_root in args.raw_root:
        for path in raw_root.rglob(pattern):
            sample_id = path.parent.name
            if sample_id not in samples:
                continue
            if sample_id in found:
                raise ValueError(f"duplicate output for {sample_id}: {found[sample_id]} and {path}")
            found[sample_id] = path

    rows = []
    for sample_id, sample in samples.items():
        path = found.get(sample_id)
        if path is None:
            rows.append({"sample_id": sample_id, "target_id": sample["target_id"], "model": args.model, "status": "missing"})
            continue
        native = load_native_affinity(path, args.model)
        row = {
            "sample_id": sample_id,
            "target_id": sample["target_id"],
            "model": args.model,
            "status": "complete",
            "native_output_path": str(path),
        }
        row.update(native)
        rows.append(row)

    fieldnames = sorted({key for row in rows for key in row})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    summary = {"model": args.model, "complete": len(found), "expected": len(samples), "output": str(args.output)}
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
