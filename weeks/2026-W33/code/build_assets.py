#!/usr/bin/env python3
"""Snapshot the canonical bibliography for the 2026-W33 report."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parents[3]
    report_dir = root / "weeks/2026-W33/report"
    destination = report_dir / "figures"
    destination.mkdir(parents=True, exist_ok=True)
    records = []
    bibliography = root / "literature/references.bib"
    bibliography_snapshot = report_dir / "references.bib"
    shutil.copy2(bibliography, bibliography_snapshot)
    records.append(
        {
            "source": str(bibliography.relative_to(root)),
            "source_sha256": sha256(bibliography),
            "snapshot": str(bibliography_snapshot.relative_to(root)),
            "snapshot_sha256": sha256(bibliography_snapshot),
        }
    )
    (destination / "assets.json").write_text(json.dumps(records, indent=2) + "\n")
    print("Collected the canonical bibliography; this report has no figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
