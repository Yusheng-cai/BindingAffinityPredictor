"""Load and validate the small canonical JSON affinity manifest."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

from affinity_benchmark.data.affinity import (
    concentration_to_log10_micromolar,
    concentration_to_px,
)


_STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_MEASUREMENT_TYPES = frozenset({"ki", "kd", "ic50"})
_QUALIFIERS = frozenset({"=", "<", ">", "<=", ">="})
_PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def load_manifest(path: str | Path) -> dict[str, Any]:
    """Load JSON from *path*, validate it, and return its parsed object."""

    manifest_path = Path(path)
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    validate_manifest(manifest)
    return manifest


def validate_manifest(manifest: Any) -> None:
    """Raise ``ValueError`` with all dependency-free schema violations."""

    issues: list[str] = []
    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a JSON object")

    if manifest.get("schema_version") != 1:
        issues.append("schema_version must equal 1")
    if not _nonempty_string(manifest.get("manifest_id")):
        issues.append("manifest_id must be a non-empty string")

    samples = manifest.get("samples")
    if not isinstance(samples, list) or not samples:
        issues.append("samples must be a non-empty list")
        _raise_issues(issues)
        return

    seen_ids: set[str] = set()
    for index, sample in enumerate(samples):
        prefix = f"samples[{index}]"
        if not isinstance(sample, dict):
            issues.append(f"{prefix} must be an object")
            continue

        sample_id = sample.get("sample_id")
        if not _nonempty_string(sample_id):
            issues.append(f"{prefix}.sample_id must be a non-empty string")
        elif sample_id in seen_ids:
            issues.append(f"{prefix}.sample_id duplicates {sample_id!r}")
        else:
            seen_ids.add(sample_id)

        _validate_protein(sample.get("protein"), prefix, issues)
        _validate_ligand(sample.get("ligand"), prefix, issues)
        _validate_measurement(sample.get("measurement"), prefix, issues)
        _validate_structure_reference(sample.get("structure_reference"), prefix, issues)

    _raise_issues(issues)


def _validate_protein(protein: Any, prefix: str, issues: list[str]) -> None:
    path = f"{prefix}.protein"
    if not isinstance(protein, dict):
        issues.append(f"{path} must be an object")
        return

    sequence = protein.get("sequence")
    if not _nonempty_string(sequence):
        issues.append(f"{path}.sequence must be a non-empty string")
        return
    if sequence != sequence.upper() or set(sequence) - _STANDARD_AMINO_ACIDS:
        invalid = "".join(sorted(set(sequence) - _STANDARD_AMINO_ACIDS))
        issues.append(f"{path}.sequence must contain uppercase standard amino acids; invalid={invalid!r}")

    sequence_length = protein.get("sequence_length")
    if isinstance(sequence_length, bool) or not isinstance(sequence_length, int):
        issues.append(f"{path}.sequence_length must be an integer")
    elif sequence_length != len(sequence):
        issues.append(
            f"{path}.sequence_length is {sequence_length}, but the sequence has {len(sequence)} residues"
        )


def _validate_ligand(ligand: Any, prefix: str, issues: list[str]) -> None:
    path = f"{prefix}.ligand"
    if not isinstance(ligand, dict):
        issues.append(f"{path} must be an object")
        return
    if not _nonempty_string(ligand.get("input_smiles")):
        issues.append(f"{path}.input_smiles must be a non-empty string")


def _validate_measurement(measurement: Any, prefix: str, issues: list[str]) -> None:
    path = f"{prefix}.measurement"
    if not isinstance(measurement, dict):
        issues.append(f"{path} must be an object")
        return

    measurement_type = measurement.get("type")
    if measurement_type not in _MEASUREMENT_TYPES:
        issues.append(f"{path}.type must be one of {sorted(_MEASUREMENT_TYPES)}")
    if measurement.get("qualifier") not in _QUALIFIERS:
        issues.append(f"{path}.qualifier must be one of {sorted(_QUALIFIERS)}")

    value = measurement.get("value")
    unit = measurement.get("unit")
    try:
        expected_log = concentration_to_log10_micromolar(value, unit)
        expected_px = concentration_to_px(value, unit)
    except (TypeError, ValueError) as error:
        issues.append(f"{path} has invalid value/unit: {error}")
        return

    error_value = measurement.get("error")
    if error_value is not None:
        if isinstance(error_value, bool) or not isinstance(error_value, (int, float)):
            issues.append(f"{path}.error must be a non-negative finite number or null")
        elif not math.isfinite(float(error_value)) or error_value < 0:
            issues.append(f"{path}.error must be a non-negative finite number or null")

    _check_close(measurement.get("log10_value_uM"), expected_log, f"{path}.log10_value_uM", issues)
    _check_close(measurement.get("pX_M"), expected_px, f"{path}.pX_M", issues)


def _validate_structure_reference(reference: Any, prefix: str, issues: list[str]) -> None:
    path = f"{prefix}.structure_reference"
    if not isinstance(reference, dict):
        issues.append(f"{path} must be an object")
        return
    pdb_id = reference.get("pdb_id")
    if not _nonempty_string(pdb_id) or _PDB_ID.fullmatch(pdb_id) is None:
        issues.append(f"{path}.pdb_id must be a four-character PDB identifier")
    if reference.get("use_as_model_input") is not False:
        issues.append(f"{path}.use_as_model_input must be false for this blind-pose protocol")


def _check_close(actual: Any, expected: float, path: str, issues: list[str]) -> None:
    if isinstance(actual, bool) or not isinstance(actual, (int, float)):
        issues.append(f"{path} must be a finite number")
        return
    if not math.isfinite(float(actual)) or not math.isclose(float(actual), expected, abs_tol=1e-12):
        issues.append(f"{path}={actual!r}, expected {expected!r} from value and unit")


def _nonempty_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _raise_issues(issues: list[str]) -> None:
    if issues:
        formatted = "\n".join(f"- {issue}" for issue in issues)
        raise ValueError(f"manifest validation failed:\n{formatted}")
