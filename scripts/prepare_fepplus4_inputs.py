#!/usr/bin/env python3
"""Create model-native YAML inputs from the frozen FEP+4 manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from affinity_benchmark.adapters.affinity_models import affinity_input_document, protein_chains
from affinity_benchmark.data.manifest import load_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--model", choices=("boltz2", "nesso1"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--msa-root", type=Path)
    parser.add_argument("--target", action="append", dest="targets")
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    selected = set(args.targets or [])
    count = 0
    for sample in manifest["samples"]:
        target = sample["target_id"]
        if selected and target not in selected:
            continue
        msa_paths = None
        if args.model == "boltz2":
            if args.msa_root is None:
                raise ValueError("--msa-root is required for Boltz-2 inputs")
            msa_paths = {}
            for chain in protein_chains(sample):
                path = (args.msa_root / target / f"{target}_{chain['id']}.csv").resolve()
                if not path.exists():
                    raise FileNotFoundError(path)
                msa_paths[chain["id"]] = str(path)
        document = affinity_input_document(sample, args.model, msa_paths=msa_paths)
        target_dir = args.output_dir / target
        target_dir.mkdir(parents=True, exist_ok=True)
        path = target_dir / f"{sample['sample_id']}.yaml"
        path.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
        count += 1
    print(f"wrote {count} {args.model} inputs beneath {args.output_dir}")


if __name__ == "__main__":
    main()
