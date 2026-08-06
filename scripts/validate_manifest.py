#!/usr/bin/env python3
"""Validate a canonical affinity manifest without downloading anything."""

from __future__ import annotations

import argparse
from pathlib import Path

from affinity_benchmark.data.manifest import load_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    manifest = load_manifest(args.manifest)
    print(f"valid: {manifest['manifest_id']} ({len(manifest['samples'])} sample(s))")


if __name__ == "__main__":
    main()
