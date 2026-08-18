#!/usr/bin/env python3
"""Export a compact, browser-native view of one RFdiffusion/ProteinMPNN run.

The output is JavaScript rather than JSON so that the tutorial continues to work
when index.html is opened directly from disk (where fetch() is commonly blocked).
C-alpha trajectories, final N/C-alpha/C/O backbone atoms, and compact summaries
are exported; raw model outputs remain under the git-ignored runs/ tree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import re
from pathlib import Path
from typing import Any

import numpy as np


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA1_TO_3 = {one: three for three, one in AA3_TO_1.items()}
BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_ca_models(path: Path) -> list[list[dict[str, Any]]]:
    """Read C-alpha atoms from a single- or multi-model PDB.

    RFdiffusion's trajectory writer separates models with ENDMDL but does not
    write MODEL records, so ENDMDL is the authoritative delimiter here.
    """
    models: list[list[dict[str, Any]]] = []
    current: list[dict[str, Any]] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                current.append(
                    {
                        "chain": line[21].strip(),
                        "residue": int(line[22:26]),
                        "name": line[17:20].strip(),
                        "coord": [
                            round(float(line[30:38]), 3),
                            round(float(line[38:46]), 3),
                            round(float(line[46:54]), 3),
                        ],
                    }
                )
            elif line.startswith("ENDMDL") and current:
                models.append(current)
                current = []
    if current:
        models.append(current)
    if not models:
        raise ValueError(f"No C-alpha coordinates found in {path}")
    return models


def parse_final_backbone(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name not in BACKBONE_ATOMS:
                continue
            atoms.append(
                {
                    "chain": line[21].strip(),
                    "residue": int(line[22:26]),
                    "atom": atom_name,
                    "element": atom_name[0],
                    "coord": [
                        round(float(line[30:38]), 3),
                        round(float(line[38:46]), 3),
                        round(float(line[46:54]), 3),
                    ],
                }
            )
    if not atoms:
        raise ValueError(f"No final backbone atoms found in {path}")
    return atoms


def chain_atoms(model: list[dict[str, Any]], chain: str) -> list[dict[str, Any]]:
    atoms = [atom for atom in model if atom["chain"] == chain]
    if not atoms:
        raise ValueError(f"Chain {chain!r} absent from parsed PDB model")
    return atoms


def parse_mpnn_fasta(path: Path) -> list[dict[str, Any]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    candidates: list[dict[str, Any]] = []
    for index in range(0, len(lines), 2):
        header, sequence = lines[index], lines[index + 1]
        if "sample=" not in header:
            continue
        fields = dict(re.findall(r"(sample|score|global_score)=([^,]+)", header))
        candidates.append(
            {
                "sample": int(fields["sample"]),
                "score": round(float(fields["score"]), 4),
                "globalScore": round(float(fields["global_score"]), 4),
                "sequence": sequence,
            }
        )
    if not candidates:
        raise ValueError(f"No ProteinMPNN samples found in {path}")
    return candidates


def nearest_contacts(
    target: list[dict[str, Any]], binder: list[dict[str, Any]], limit: int = 14
) -> list[dict[str, Any]]:
    target_xyz = np.asarray([atom["coord"] for atom in target], dtype=float)
    binder_xyz = np.asarray([atom["coord"] for atom in binder], dtype=float)
    distances = np.linalg.norm(
        binder_xyz[:, None, :] - target_xyz[None, :, :], axis=-1
    )
    pairs: list[tuple[float, int, int]] = []
    used_binder: set[int] = set()
    for flat_index in np.argsort(distances, axis=None):
        binder_index, target_index = np.unravel_index(flat_index, distances.shape)
        distance = float(distances[binder_index, target_index])
        if distance > 9.0:
            break
        if int(binder_index) in used_binder:
            continue
        pairs.append((distance, int(binder_index), int(target_index)))
        used_binder.add(int(binder_index))
        if len(pairs) == limit:
            break
    return [
        {
            "binderIndex": binder_index,
            "targetIndex": target_index,
            "binderResidue": binder[binder_index]["residue"],
            "targetResidue": target[target_index]["residue"],
            "distanceAngstrom": round(distance, 2),
        }
        for distance, binder_index, target_index in pairs
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target-chain", default="A")
    parser.add_argument("--binder-chain", default="B")
    args = parser.parse_args()

    rf_root = args.run_root / "rfdiffusion/seed42/raw"
    mpnn_root = args.run_root / "proteinmpnn/seed42/raw"
    final_pdb = rf_root / "design_ppi_42.pdb"
    xt_path = rf_root / "traj/design_ppi_42_Xt-1_traj.pdb"
    px0_path = rf_root / "traj/design_ppi_42_pX0_traj.pdb"
    trb_path = rf_root / "design_ppi_42.trb"
    fasta_path = mpnn_root / "seqs/design_ppi_42.fa"
    scores_path = mpnn_root / "scores/design_ppi_42.npz"

    final_model = parse_ca_models(final_pdb)[0]
    final_backbone = parse_final_backbone(final_pdb)
    target = chain_atoms(final_model, args.target_chain)
    binder = chain_atoms(final_model, args.binder_chain)
    xt_models = parse_ca_models(xt_path)
    px0_models = parse_ca_models(px0_path)
    if len(xt_models) != 50 or len(px0_models) != 50:
        raise ValueError(
            f"Expected 50 models per trajectory; got {len(xt_models)} and {len(px0_models)}"
        )

    # RFdiffusion writes trajectories final-to-initial for PyMOL playback. The
    # browser presents the generative direction, t=50 down to t=1, so reverse.
    xt_frames = [
        [atom["coord"] for atom in chain_atoms(model, args.binder_chain)]
        for model in reversed(xt_models)
    ]
    px0_frames = [
        [atom["coord"] for atom in chain_atoms(model, args.binder_chain)]
        for model in reversed(px0_models)
    ]
    if any(len(frame) != len(binder) for frame in xt_frames + px0_frames):
        raise ValueError("Binder length changed within a trajectory")

    # The pickle is a trusted local artifact generated by the pinned RFdiffusion
    # checkout. Never unpickle untrusted external files.
    with trb_path.open("rb") as handle:
        trb = pickle.load(handle)
    plddt = np.asarray(trb["plddt"], dtype=float)
    binder_length = len(binder)
    if plddt.shape != (50, len(final_model)):
        raise ValueError(f"Unexpected confidence array shape: {plddt.shape}")
    confidence = [
        {
            "timestep": 50 - index,
            "binderMean": round(float(row[:binder_length].mean()), 4),
            "targetMean": round(float(row[binder_length:].mean()), 4),
            "complexMean": round(float(row.mean()), 4),
        }
        for index, row in enumerate(plddt)
    ]

    candidates = parse_mpnn_fasta(fasta_path)
    selected = min(candidates, key=lambda candidate: candidate["score"])
    if len(selected["sequence"]) != binder_length:
        raise ValueError(
            f"Selected sequence length {len(selected['sequence'])} does not match "
            f"binder length {binder_length}"
        )
    with np.load(scores_path) as score_archive:
        archived_scores = score_archive["score"].astype(float).tolist()
    if not np.allclose(
        archived_scores, [candidate["score"] for candidate in candidates], atol=5e-5
    ):
        raise ValueError("FASTA and ProteinMPNN score archive disagree")

    hotspot_residues = [59, 83, 91]
    hotspot_atoms = [atom for atom in target if atom["residue"] in hotspot_residues]
    if len(hotspot_atoms) != len(hotspot_residues):
        raise ValueError("Not all configured hotspot residues were found")

    payload = {
        "schemaVersion": 2,
        "status": "RFdiffusion and ProteinMPNN completed; independent structure validation not run",
        "metadata": {
            "experiment": "exp003_rfdiffusion_insr_binder_smoke",
            "target": "insulin-receptor example, residues A1-A150",
            "seed": 42,
            "contig": "[A1-150/0 70-100]",
            "sampledBinderLength": binder_length,
            "hotspots": ["A59", "A83", "A91"],
            "diffusionSteps": 50,
            "gpu": trb["device"],
            "modelReportedSeconds": round(float(trb["time"]), 2),
            "rfdiffusionRevision": "86507b6538f51fce57b5a72477165f03999ed7ae",
            "proteinmpnnRevision": "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57",
        },
        "target": {
            "chain": args.target_chain,
            "residues": [atom["residue"] for atom in target],
            "sequence": "".join(AA3_TO_1.get(atom["name"], "X") for atom in target),
            "ca": [atom["coord"] for atom in target],
            "backboneAtoms": [
                atom for atom in final_backbone if atom["chain"] == args.target_chain
            ],
            "hotspots": [
                {"label": f"A{atom['residue']}", "residue": atom["residue"], "coord": atom["coord"]}
                for atom in hotspot_atoms
            ],
        },
        "binder": {
            "chain": args.binder_chain,
            "length": binder_length,
            "finalCa": [atom["coord"] for atom in binder],
            "backboneAtoms": [
                {
                    **atom,
                    "oneLetter": selected["sequence"][atom["residue"] - 1],
                    "residueName": AA1_TO_3[selected["sequence"][atom["residue"] - 1]],
                }
                for atom in final_backbone
                if atom["chain"] == args.binder_chain
            ],
            "assignedResidues": [
                {
                    "residue": index,
                    "oneLetter": one_letter,
                    "residueName": AA1_TO_3[one_letter],
                }
                for index, one_letter in enumerate(selected["sequence"], start=1)
            ],
            "coordinateCompleteness": "backbone_N_CA_C_O_only",
            "sidechainsPacked": False,
            "timesteps": list(range(50, 0, -1)),
            "trajectories": {"xtMinus1": xt_frames, "pX0": px0_frames},
            "confidence": confidence,
        },
        "contacts": nearest_contacts(target, binder),
        "proteinMpnn": {
            "model": "v_48_020",
            "temperature": 0.1,
            "scoreMeaning": "mean negative log-probability over designed chain; lower is better within this run",
            "candidates": candidates,
            "selectedSample": selected["sample"],
        },
        "provenance": {
            "finalPdbSha256": sha256(final_pdb),
            "xtTrajectorySha256": sha256(xt_path),
            "px0TrajectorySha256": sha256(px0_path),
            "trbSha256": sha256(trb_path),
            "fastaSha256": sha256(fasta_path),
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(payload, separators=(",", ":"), allow_nan=False)
    args.output.write_text(
        "// Generated by scripts/export_rfdiffusion_tutorial_data.py\n"
        f"window.RFDIFFUSION_REAL_RUN={serialized};\n"
    )
    print(f"Wrote {args.output} ({args.output.stat().st_size:,} bytes)")
    print(
        f"Validated target={len(target)}, binder={binder_length}, "
        f"trajectories={len(xt_frames)}/{len(px0_frames)}, sequences={len(candidates)}"
    )


if __name__ == "__main__":
    main()
