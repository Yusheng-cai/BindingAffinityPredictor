#!/usr/bin/env python3
"""Fetch the frozen RCSB mmCIF and CCD files for exp008."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPOSITORY_ROOT / "data/manifests/fepplus4_crystal16.json"
DEFAULT_OUTPUT = REPOSITORY_ROOT / "data/raw/exp008_boltz2_crystal_pose"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url, headers={"User-Agent": "BindingAffinityPredictor-exp008/1.0"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = response.read()
    if not payload:
        raise ValueError(f"empty response from {url}")
    destination.write_bytes(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    records: list[dict[str, object]] = []
    unique_pdb = sorted({sample["pdb_id"].upper() for sample in manifest["samples"]})
    unique_ccd = sorted({sample["ccd_id"].upper() for sample in manifest["samples"]})

    requests: list[tuple[str, Path, str, str]] = []
    for pdb_id in unique_pdb:
        requests.append(
            (
                f"https://files.rcsb.org/download/{pdb_id}.cif",
                args.output_dir / "references" / f"{pdb_id.lower()}.cif",
                "pdb_mmcif",
                pdb_id,
            )
        )
    for ccd_id in unique_ccd:
        requests.extend(
            [
                (
                    f"https://files.rcsb.org/ligands/download/{ccd_id}.cif",
                    args.output_dir / "ccd" / f"{ccd_id}.cif",
                    "ccd_cif",
                    ccd_id,
                ),
                (
                    f"https://files.rcsb.org/ligands/download/{ccd_id}_ideal.sdf",
                    args.output_dir / "ccd" / f"{ccd_id}_ideal.sdf",
                    "ccd_ideal_sdf",
                    ccd_id,
                ),
            ]
        )

    for url, destination, kind, identifier in requests:
        if args.refresh or not destination.exists():
            fetch(url, destination)
        records.append(
            {
                "kind": kind,
                "identifier": identifier,
                "url": url,
                "path": str(destination.relative_to(REPOSITORY_ROOT)),
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
            }
        )

    provenance = {
        "manifest_id": manifest["manifest_id"],
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "records": records,
    }
    provenance_path = args.output_dir / "download_manifest.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Verified {len(records)} files; provenance: {provenance_path}")


if __name__ == "__main__":
    main()
