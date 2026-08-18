#!/usr/bin/env python3
"""Generate one reusable Boltz CSV MSA per unique target protein chain."""

from __future__ import annotations

import argparse
from pathlib import Path

from boltz.main import compute_msa

from affinity_benchmark.adapters.affinity_models import protein_chains
from affinity_benchmark.data.manifest import load_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target", action="append", dest="targets")
    parser.add_argument("--server-url", default="https://api.colabfold.com")
    parser.add_argument("--pairing-strategy", default="greedy")
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    selected = set(args.targets or [])
    representative = {}
    for sample in manifest["samples"]:
        representative.setdefault(sample["target_id"], sample)

    for target, sample in representative.items():
        if selected and target not in selected:
            continue
        target_dir = args.output_dir / target
        target_dir.mkdir(parents=True, exist_ok=True)
        chains = protein_chains(sample)
        expected = [target_dir / f"{target}_{chain['id']}.csv" for chain in chains]
        if all(path.exists() for path in expected):
            print(f"{target}: reusing {len(expected)} existing MSA file(s)")
            continue
        data = {f"{target}_{chain['id']}": chain["sequence"] for chain in chains}
        compute_msa(
            data=data,
            target_id=target,
            msa_dir=target_dir,
            msa_server_url=args.server_url,
            msa_pairing_strategy=args.pairing_strategy,
        )
        missing = [str(path) for path in expected if not path.exists()]
        if missing:
            raise FileNotFoundError(f"MSA generation did not produce {missing}")


if __name__ == "__main__":
    main()
