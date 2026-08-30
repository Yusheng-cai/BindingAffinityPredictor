#!/usr/bin/env python3
"""Generate one reusable Boltz CSV MSA per unique target protein chain."""

from __future__ import annotations

import argparse
import tarfile
import time
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
    parser.add_argument("--retries", type=int, default=3)
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
        if args.retries < 1:
            raise ValueError("--retries must be at least 1")
        for attempt in range(1, args.retries + 1):
            try:
                compute_msa(
                    data=data,
                    target_id=target,
                    msa_dir=target_dir,
                    msa_server_url=args.server_url,
                    msa_pairing_strategy=args.pairing_strategy,
                )
                break
            except Exception as error:
                if isinstance(error, tarfile.ReadError):
                    for archive in target_dir.glob("*_tmp*/out.tar.gz"):
                        quarantined = archive.with_name(
                            f"{archive.name}.invalid_attempt{attempt}"
                        )
                        archive.replace(quarantined)
                        print(
                            f"{target}: quarantined malformed server response "
                            f"as {quarantined}"
                        )
                if attempt == args.retries:
                    raise
                print(
                    f"{target}: MSA attempt {attempt}/{args.retries} failed "
                    f"({type(error).__name__}: {error}); retrying"
                )
                time.sleep(2 * attempt)
        missing = [str(path) for path in expected if not path.exists()]
        if missing:
            raise FileNotFoundError(f"MSA generation did not produce {missing}")


if __name__ == "__main__":
    main()
