#!/usr/bin/env python3
"""Write one sequence-and-SMILES Nesso YAML input per manifest sample."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from affinity_benchmark.data.manifest import validate_manifest


def nesso_record(sample: dict) -> dict:
    chains = sample["protein"]["chains"]
    if len(chains) != 1:
        raise ValueError(f"{sample['sample_id']}: expected exactly one protein chain")
    return {
        "sequences": [
            {
                "protein": {
                    "id": chains[0]["id"],
                    "sequence": chains[0]["sequence"],
                }
            },
            {
                "ligand": {
                    "id": "L",
                    "smiles": sample["ligand"]["input_smiles"],
                }
            },
        ],
        "properties": [{"affinity": {"binder": "L"}}],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    validate_manifest(manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    written = 0
    for sample in manifest["samples"]:
        if not sample.get("included", True):
            continue
        output = args.output_dir / f"{sample['sample_id']}.yaml"
        output.write_text(
            yaml.safe_dump(nesso_record(sample), sort_keys=False),
            encoding="utf-8",
        )
        written += 1
    print(f"Wrote {written} Nesso inputs to {args.output_dir}")


if __name__ == "__main__":
    main()
