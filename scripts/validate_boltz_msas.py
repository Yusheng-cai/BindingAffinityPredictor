#!/usr/bin/env python3
"""Validate and inventory Boltz CSV MSAs against a frozen manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_msa(path: Path, query_sequence: str) -> dict:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["key", "sequence"]:
            raise ValueError(f"{path}: expected key,sequence columns")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: empty MSA")
    observed_query = rows[0]["sequence"].replace("-", "").upper()
    if observed_query != query_sequence:
        raise ValueError(f"{path}: first MSA sequence does not match manifest query")
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "rows": len(rows),
        "query_length": len(query_sequence),
        "query_matches_manifest": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--msa-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    records = []
    for sample in manifest["samples"]:
        chain = sample["protein"]["chains"][0]
        target = sample["target_id"]
        path = args.msa_root / target / f"{target}_{chain['id']}.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        records.append(
            {
                "sample_id": sample["sample_id"],
                "target_id": target,
                "chain_id": chain["id"],
                **inspect_msa(path, chain["sequence"]),
            }
        )
    result = {
        "schema_version": 1,
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": sha256(args.manifest),
        "server": "https://api.colabfold.com",
        "pairing_strategy": "greedy",
        "systems": len(records),
        "minimum_rows": min(record["rows"] for record in records),
        "median_rows": statistics.median(record["rows"] for record in records),
        "maximum_rows": max(record["rows"] for record in records),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"Validated {len(records)} MSAs; depth range "
        f"{result['minimum_rows']}–{result['maximum_rows']} rows"
    )


if __name__ == "__main__":
    main()
