#!/usr/bin/env python3
"""Audit meaningful Boltz-2 intermediate and output artifacts for one run.

Run this with the isolated Boltz environment, not the repository's base Python.
The script reconstructs model-ready feature tensors deterministically from
Boltz's processed artifacts, but does not execute the neural network again.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
import torch

from boltz.data.module.inferencev2 import PredictionDataset
from boltz.data.types import MSA, Manifest, StructureV2


def array_summary(array: np.ndarray) -> dict[str, Any]:
    """Return shape, dtype, and finite numeric statistics for an array."""

    summary: dict[str, Any] = {
        "shape": list(array.shape),
        "dtype": str(array.dtype),
    }
    if array.size and np.issubdtype(array.dtype, np.number):
        summary.update(
            {
                "all_finite": bool(np.isfinite(array).all()),
                "minimum": float(np.min(array)),
                "maximum": float(np.max(array)),
                "mean": float(np.mean(array)),
            }
        )
    return summary


def feature_shapes(features: dict[str, Any]) -> dict[str, Any]:
    """Describe model-ready, unbatched features without serializing tensors."""

    result: dict[str, Any] = {}
    for name, value in features.items():
        if name == "record":
            result[name] = {"record_id": value.id}
        elif isinstance(value, torch.Tensor):
            result[name] = {"shape": list(value.shape), "dtype": str(value.dtype)}
        elif isinstance(value, list):
            result[name] = {"type": "list", "length": len(value)}
        else:
            result[name] = {"type": type(value).__name__}
    return result


def sha256(path: Path) -> str:
    """Hash one artifact in bounded chunks."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def contacts_from_cif(path: Path, cutoff: float = 4.0) -> list[dict[str, Any]]:
    """List protein residues with a heavy-atom contact to ligand chain B."""

    structure = gemmi.make_structure_from_block(gemmi.cif.read(str(path)).sole_block())
    model = structure[0]
    protein = model["A"]
    ligand = model["B"]
    ligand_atoms = [
        atom
        for residue in ligand
        for atom in residue
        if atom.element.name != "H"
    ]
    contacts: list[dict[str, Any]] = []
    for residue in protein:
        distances = [
            atom.pos.dist(ligand_atom.pos)
            for atom in residue
            if atom.element.name != "H"
            for ligand_atom in ligand_atoms
        ]
        if distances and min(distances) <= cutoff:
            contacts.append(
                {
                    "chain": protein.name,
                    "residue_name": residue.name,
                    "residue_number": residue.seqid.num,
                    "minimum_heavy_atom_distance_angstrom": min(distances),
                }
            )
    return contacts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--mol-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    base = args.results_dir.resolve()
    processed = base / "processed"
    prediction_root = base / "predictions"
    manifest = Manifest.load(processed / "manifest.json")
    if len(manifest.records) != 1:
        raise ValueError("this focused audit expects exactly one record")
    record = manifest.records[0]
    record_id = record.id
    prediction = prediction_root / record_id

    msa = MSA.load(processed / "msa" / f"{record_id}_0.npz")
    molecular_system = StructureV2.load(processed / "structures" / f"{record_id}.npz")

    structure_dataset = PredictionDataset(
        manifest=manifest,
        target_dir=processed / "structures",
        msa_dir=processed / "msa",
        mol_dir=args.mol_dir.resolve(),
        constraints_dir=processed / "constraints",
        template_dir=processed / "templates",
        extra_mols_dir=processed / "mols",
        affinity=False,
    )
    affinity_dataset = PredictionDataset(
        manifest=manifest,
        target_dir=prediction_root,
        msa_dir=processed / "msa",
        mol_dir=args.mol_dir.resolve(),
        constraints_dir=processed / "constraints",
        template_dir=processed / "templates",
        extra_mols_dir=processed / "mols",
        override_method="other",
        affinity=True,
    )

    structure_features = structure_dataset[0]
    affinity_features = affinity_dataset[0]

    embeddings_path = prediction / f"embeddings_{record_id}.npz"
    plddt_path = prediction / f"plddt_{record_id}_model_0.npz"
    pae_path = prediction / f"pae_{record_id}_model_0.npz"
    pde_path = prediction / f"pde_{record_id}_model_0.npz"
    cif_path = prediction / f"{record_id}_model_0.cif"
    confidence_path = prediction / f"confidence_{record_id}_model_0.json"
    affinity_path = prediction / f"affinity_{record_id}.json"

    with np.load(embeddings_path) as arrays:
        embedding_summary = {name: array_summary(arrays[name]) for name in arrays.files}
    with np.load(plddt_path) as arrays:
        plddt = arrays["plddt"]
    with np.load(pae_path) as arrays:
        pae = arrays["pae"]
    with np.load(pde_path) as arrays:
        pde = arrays["pde"]

    protein_tokens = record.chains[0].num_residues
    total_tokens = int(plddt.shape[0])
    ligand_slice = slice(protein_tokens, total_tokens)
    protein_slice = slice(0, protein_tokens)

    artifacts = [
        embeddings_path,
        plddt_path,
        pae_path,
        pde_path,
        cif_path,
        confidence_path,
        affinity_path,
        processed / "manifest.json",
        processed / "msa" / f"{record_id}_0.npz",
        processed / "structures" / f"{record_id}.npz",
    ]

    report = {
        "record_id": record_id,
        "preprocessing": {
            "chains": len(molecular_system.chains),
            "residues": len(molecular_system.residues),
            "atoms": len(molecular_system.atoms),
            "bonds": len(molecular_system.bonds),
            "protein_residues": protein_tokens,
            "ligand_atom_tokens": total_tokens - protein_tokens,
            "msa_sequences": len(msa.sequences),
            "msa_aligned_residue_symbols": len(msa.residues),
            "msa_deletion_records": len(msa.deletions),
        },
        "model_ready_structure_feature_shapes_unbatched": feature_shapes(structure_features),
        "model_ready_affinity_feature_shapes_unbatched": feature_shapes(affinity_features),
        "trunk_embeddings": embedding_summary,
        "confidence": json.loads(confidence_path.read_text()),
        "confidence_arrays": {
            "plddt": {
                **array_summary(plddt),
                "protein_mean": float(np.mean(plddt[protein_slice])),
                "ligand_mean": float(np.mean(plddt[ligand_slice])),
            },
            "pae": {
                **array_summary(pae),
                "protein_protein_mean": float(np.mean(pae[protein_slice, protein_slice])),
                "protein_to_ligand_mean": float(np.mean(pae[protein_slice, ligand_slice])),
                "ligand_to_protein_mean": float(np.mean(pae[ligand_slice, protein_slice])),
                "ligand_ligand_mean": float(np.mean(pae[ligand_slice, ligand_slice])),
            },
            "pde": {
                **array_summary(pde),
                "protein_protein_mean": float(np.mean(pde[protein_slice, protein_slice])),
                "protein_ligand_mean": float(np.mean(pde[protein_slice, ligand_slice])),
                "ligand_ligand_mean": float(np.mean(pde[ligand_slice, ligand_slice])),
            },
        },
        "affinity_native": json.loads(affinity_path.read_text()),
        "predicted_contacts_at_4_angstrom": contacts_from_cif(cif_path),
        "artifact_integrity": {
            str(path.relative_to(base)): {
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in artifacts
        },
        "interpretation_limits": [
            "Feature shapes are reconstructed deterministically from saved processed inputs; the neural network is not rerun.",
            "Saved embeddings are final trunk outputs, not per-layer activations.",
            "Transient coordinates from all diffusion steps are not emitted by the stock CLI.",
            "Native affinity fields are recorded without claiming accuracy when no matched quantitative assay is available.",
        ],
    }

    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
