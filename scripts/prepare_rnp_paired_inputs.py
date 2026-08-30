#!/usr/bin/env python3
"""Prepare structure-only Nesso-1 or Boltz-2 inputs from a frozen RnP manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from affinity_benchmark.data.manifest import load_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--model", choices=("nesso1", "boltz2"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--msa-root", type=Path)
    parser.add_argument("--sample-id", action="append", dest="sample_ids")
    return parser.parse_args()


def document_for_sample(sample: dict, model: str, msa_root: Path | None) -> dict:
    chains = sample["protein"]["chains"]
    if len(chains) != 1:
        raise ValueError(f"{sample['sample_id']}: expected one protein chain")
    chain = chains[0]
    protein = {"id": "A", "sequence": chain["sequence"]}
    if model == "boltz2":
        if msa_root is None:
            raise ValueError("--msa-root is required for Boltz-2")
        target = sample["target_id"]
        msa_path = (
            msa_root / target / f"{target}_{chain['id']}.csv"
        ).resolve()
        if not msa_path.is_file():
            raise FileNotFoundError(msa_path)
        protein["msa"] = str(msa_path)

    document = {
        "sequences": [
            {"protein": protein},
            {
                "ligand": {
                    "id": "B",
                    "smiles": sample["ligand"]["input_smiles"],
                }
            },
        ]
    }
    if model == "boltz2":
        document = {"version": 1, **document}
    return document


def main() -> None:
    args = parse_args()
    manifest = load_manifest(args.manifest)
    selected_ids = set(args.sample_ids or [])
    known_ids = {sample["sample_id"] for sample in manifest["samples"]}
    unknown = sorted(selected_ids - known_ids)
    if unknown:
        raise ValueError(f"unknown requested sample IDs: {unknown}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    written = 0
    for sample in manifest["samples"]:
        if not sample.get("included", True):
            continue
        if selected_ids and sample["sample_id"] not in selected_ids:
            continue
        document = document_for_sample(sample, args.model, args.msa_root)
        output = args.output_dir / f"{sample['sample_id']}.yaml"
        output.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
        written += 1
    print(f"Wrote {written} structure-only {args.model} inputs to {args.output_dir}")


if __name__ == "__main__":
    main()
