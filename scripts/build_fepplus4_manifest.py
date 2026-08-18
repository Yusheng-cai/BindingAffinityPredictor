#!/usr/bin/env python3
"""Build the frozen 87-compound FEP+ four-kinase benchmark manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, inchi

from affinity_benchmark.data.affinity import (
    concentration_to_log10_micromolar,
    concentration_to_px,
)
from affinity_benchmark.data.manifest import validate_manifest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = REPOSITORY_ROOT / "data" / "raw" / "fepplus4"
OUTPUT = REPOSITORY_ROOT / "data" / "manifests" / "fepplus4_87.json"
SOURCE_REVISION = "da7c3372256446222e424368be38ef3d2b55a67b"

TARGETS = {
    "cdk2": {
        "pdb_id": "1H1Q",
        "chains": [("A", "CELL DIVISION PROTEIN KINASE 2"), ("B", "CYCLIN A2")],
    },
    "tyk2": {
        "pdb_id": "4GIH",
        "chains": [("A", "Non-receptor tyrosine-protein kinase TYK2")],
    },
    "jnk1": {
        "pdb_id": "2GMX",
        "chains": [
            ("A", "Mitogen-activated protein kinase 8"),
            ("F", "C-jun-amino-terminal kinase-interacting protein 1"),
        ],
    },
    "p38": {
        "pdb_id": "3FLY",
        "chains": [("A", "Mitogen-activated protein kinase 14")],
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        else:
            sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def chains_for_target(target: str) -> list[dict]:
    fasta_path = RAW_ROOT / "rcsb" / f"{TARGETS[target]['pdb_id']}.fasta"
    records = read_fasta(fasta_path)
    chains: list[dict] = []
    for chain_id, entity_name in TARGETS[target]["chains"]:
        matches = [(header, sequence) for header, sequence in records if entity_name.lower() in header.lower()]
        if len(matches) != 1:
            raise ValueError(f"expected one {target} FASTA record matching {entity_name!r}, found {len(matches)}")
        header, sequence = matches[0]
        chains.append(
            {
                "id": chain_id,
                "name": entity_name,
                "sequence": sequence,
                "sequence_length": len(sequence),
                "sequence_source": f"https://www.rcsb.org/fasta/entry/{TARGETS[target]['pdb_id']}/display",
                "source_header": header,
            }
        )
    return chains


def normalize_measurement(entry: dict) -> dict:
    measurement = entry["measurement"]
    measurement_type = next(key for key in measurement if key not in {"comment", "doi"})
    value, error, unit = measurement[measurement_type]
    unit_scale = {"nM": 1e-3, "uM": 1.0}[unit]
    value_um = float(value) * unit_scale
    error_um = None if error is None or float(error) < 0 else float(error) * unit_scale
    return {
        "type": measurement_type,
        "value": value_um,
        "error": error_um,
        "unit": "uM",
        "source_value": value,
        "source_error": error,
        "source_unit": unit,
        "qualifier": "=",
        "log10_value_uM": concentration_to_log10_micromolar(value_um, "uM"),
        "pX_M": concentration_to_px(value_um, "uM"),
        "source_doi": measurement.get("doi"),
        "source_location": measurement.get("comment"),
    }


def build_manifest() -> dict:
    samples: list[dict] = []
    source_files: list[dict] = []
    for target in ("cdk2", "tyk2", "jnk1", "p38"):
        ligand_path = RAW_ROOT / f"{target}_ligands_da7c337.yml"
        target_path = RAW_ROOT / f"{target}_target_da7c337.yml"
        source_files.extend(
            [
                {"path": str(ligand_path.relative_to(REPOSITORY_ROOT)), "sha256": sha256(ligand_path)},
                {"path": str(target_path.relative_to(REPOSITORY_ROOT)), "sha256": sha256(target_path)},
            ]
        )
        entries = list(yaml.safe_load_all(ligand_path.read_text(encoding="utf-8")))
        chains = chains_for_target(target)
        for entry in entries:
            molecule = Chem.MolFromSmiles(entry["smiles"])
            if molecule is None:
                raise ValueError(f"RDKit could not parse {target}/{entry['name']}")
            canonical_smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
            sample_id = f"{target}_{entry['name'].removeprefix('lig_')}"
            samples.append(
                {
                    "sample_id": sample_id,
                    "target_id": target,
                    "assay_id": f"fepplus4_{target}",
                    "protein": {
                        "name": f"{target.upper()} deposited benchmark complex construct",
                        "chains": chains,
                        "construct_policy": "All deposited protein/peptide entities present in the curated benchmark complex are included by sequence.",
                    },
                    "ligand": {
                        "name": entry["name"],
                        "input_smiles": canonical_smiles,
                        "source_smiles": entry["smiles"],
                        "inchi_key": inchi.MolToInchiKey(molecule),
                        "formal_charge": Chem.GetFormalCharge(molecule),
                        "molecular_weight_da": Descriptors.MolWt(molecule),
                        "clogp": Crippen.MolLogP(molecule),
                    },
                    "measurement": normalize_measurement(entry),
                    "structure_reference": {
                        "pdb_id": TARGETS[target]["pdb_id"],
                        "use_as_model_input": False,
                        "use_for_pose_scoring": False,
                        "role": "construct_identity_and_benchmark_provenance_only",
                    },
                    "split": "fepplus4_public_benchmark",
                    "included": True,
                }
            )

    manifest = {
        "schema_version": 1,
        "manifest_id": "fepplus4_87_da7c337",
        "created_on": "2026-08-17",
        "purpose": "Released-model reproduction of the 87-neutral-compound FEP+ four-kinase affinity benchmark",
        "source_dataset": {
            "name": "OpenFF Protein-Ligand Benchmark historical four-target set",
            "source_url": "https://github.com/openforcefield/protein-ligand-benchmark",
            "source_revision": SOURCE_REVISION,
            "license": "CC-BY-4.0",
            "source_files": source_files,
            "revision_note": "The current main branch contains a reduced 73-compound form for these targets; the historical revision recovers the 16+16+21+34=87 neutral compounds reported by Boltz-2 and Nesso-1.",
        },
        "aggregation_unit": "target",
        "samples": samples,
    }
    validate_manifest(manifest)
    counts = {target: sum(sample["target_id"] == target for sample in samples) for target in TARGETS}
    if counts != {"cdk2": 16, "tyk2": 16, "jnk1": 21, "p38": 34}:
        raise ValueError(f"unexpected target counts: {counts}")
    if any(sample["ligand"]["formal_charge"] != 0 for sample in samples):
        raise ValueError("the frozen 87-compound subset must be neutral")
    return manifest


def main() -> None:
    manifest = build_manifest()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(manifest['samples'])} samples to {OUTPUT}")


if __name__ == "__main__":
    main()
