#!/usr/bin/env python3
"""Extract each Runs N' Poses system's closest pre-cutoff structural neighbor."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


DEFAULT_INPUT = Path("data/raw/runs_n_poses_metadata/all_similarity_scores.parquet")
DEFAULT_OUTPUT = Path("data/processed/runs_n_poses/boltz2_2023_similarity.csv")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cutoff", default="2023-06-01")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cutoff = pd.Timestamp(args.cutoff)
    columns = [
        "group_key",
        "target_system",
        "target_release_date",
        "sucos_shape_pocket_qcov",
    ]
    similarities = pd.read_parquet(args.input, columns=columns)
    similarities = similarities[
        similarities["target_release_date"] < cutoff
    ].copy()
    closest = (
        similarities.sort_values(
            ["group_key", "sucos_shape_pocket_qcov", "target_system"],
            ascending=[True, False, True],
        )
        .groupby("group_key", sort=True, as_index=False)
        .head(1)
        .sort_values("group_key")
        .reset_index(drop=True)
    )
    closest = closest.rename(
        columns={
            "target_system": "closest_pre_cutoff_system",
            "target_release_date": "closest_system_release_date",
            "sucos_shape_pocket_qcov": "similarity_score_0_to_100",
        }
    )
    closest.insert(1, "structural_cutoff", args.cutoff)
    closest.insert(2, "similarity_metric", "sucos_shape_pocket_qcov")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    closest.to_csv(args.output, index=False)
    print(f"Input SHA256: {sha256(args.input)}")
    print(f"Candidate target pairs before the cutoff: {len(similarities):,}")
    print(f"Unique query systems written: {len(closest):,}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
