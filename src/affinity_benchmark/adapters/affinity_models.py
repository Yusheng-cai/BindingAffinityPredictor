"""Shared input and native-output helpers for affinity-model inference."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


def protein_chains(sample: dict[str, Any]) -> list[dict[str, Any]]:
    """Return explicit protein chains while supporting legacy one-chain manifests."""

    protein = sample["protein"]
    if "chains" in protein:
        return protein["chains"]
    return [{"id": "A", "sequence": protein["sequence"]}]


def affinity_input_document(
    sample: dict[str, Any], model: str, msa_paths: dict[str, str] | None = None
) -> dict[str, Any]:
    """Build a Boltz-2 or Nesso-1 sequence-and-SMILES input document."""

    if model not in {"boltz2", "nesso1"}:
        raise ValueError(f"unsupported model {model!r}")

    sequences: list[dict[str, Any]] = []
    for chain in protein_chains(sample):
        protein = {"id": chain["id"], "sequence": chain["sequence"]}
        if model == "boltz2" and msa_paths is not None:
            protein["msa"] = msa_paths[chain["id"]]
        sequences.append({"protein": protein})

    ligand_id = "L"
    sequences.append(
        {
            "ligand": {
                "id": ligand_id,
                "smiles": sample["ligand"]["input_smiles"],
            }
        }
    )
    document: dict[str, Any] = {
        "sequences": sequences,
        "properties": [{"affinity": {"binder": ligand_id}}],
    }
    if model == "boltz2":
        document = {"version": 1, **document}
    return document


def load_native_affinity(path: str | Path, model: str) -> dict[str, Any]:
    """Parse and validate the scalar affinity fields emitted by a model."""

    with Path(path).open(encoding="utf-8") as handle:
        native = json.load(handle)
    required = {"affinity_pred_value", "affinity_probability_binary"}
    missing = sorted(required - native.keys())
    if missing:
        raise ValueError(f"{model} output {path} is missing {missing}")
    return native


def flashbind_record_id(sample: dict[str, Any]) -> str:
    """Return the upstream FEP4 record ID for one canonical manifest sample."""

    ligand_name = sample["ligand"]["name"].replace("_", "-")
    return f'{sample["target_id"]}_{ligand_name}'


def load_flashbind_predictions(path: str | Path) -> dict[str, dict[str, Any]]:
    """Load and validate a released FlashBind prediction JSON document."""

    with Path(path).open(encoding="utf-8") as handle:
        document = json.load(handle)
    if not isinstance(document, dict):
        raise ValueError("FlashBind predictions must be a JSON object keyed by record ID")

    validated: dict[str, dict[str, Any]] = {}
    for record_id, prediction in document.items():
        if not isinstance(prediction, dict):
            raise ValueError(f"FlashBind prediction {record_id!r} must be an object")
        status = prediction.get("status")
        if status not in {"success", "failed"}:
            raise ValueError(f"FlashBind prediction {record_id!r} has invalid status {status!r}")
        if status == "success":
            for field in ("pred_value", "pred_value_raw", "mw"):
                value = prediction.get(field)
                if not isinstance(value, (int, float)) or not math.isfinite(value):
                    raise ValueError(
                        f"FlashBind prediction {record_id!r} has invalid {field}: {value!r}"
                    )
        validated[record_id] = prediction
    return validated
