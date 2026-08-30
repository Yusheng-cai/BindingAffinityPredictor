#!/usr/bin/env python3
"""Extract only manifest-selected Runs N' Poses ground-truth directories."""

from __future__ import annotations

import argparse
import json
import tarfile
from pathlib import Path, PurePosixPath


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text())
    system_ids = {
        sample["runs_n_poses"]["system_id"]
        for sample in manifest["samples"]
        if sample.get("included", True)
    }
    selected_members: list[tarfile.TarInfo] = []
    found_systems: set[str] = set()
    with tarfile.open(args.archive, "r:gz") as archive:
        for member in archive.getmembers():
            path = PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(f"unsafe archive member: {member.name}")
            if len(path.parts) >= 2 and path.parts[0] == "ground_truth":
                system_id = path.parts[1]
                if system_id in system_ids:
                    selected_members.append(member)
                    found_systems.add(system_id)
        missing = sorted(system_ids - found_systems)
        if missing:
            raise FileNotFoundError(f"archive is missing selected systems: {missing}")
        args.output_dir.mkdir(parents=True, exist_ok=True)
        archive.extractall(args.output_dir, members=selected_members)

    required_missing = []
    for system_id in sorted(system_ids):
        system_dir = args.output_dir / "ground_truth" / system_id
        for name in ("system.cif", "receptor.cif", "sequences.fasta"):
            if not (system_dir / name).is_file():
                required_missing.append(str(system_dir / name))
    if required_missing:
        raise FileNotFoundError(f"required extracted files are missing: {required_missing}")
    print(f"Extracted {len(found_systems)} systems to {args.output_dir / 'ground_truth'}")


if __name__ == "__main__":
    main()
