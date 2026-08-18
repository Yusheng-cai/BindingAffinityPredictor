#!/usr/bin/env python3
"""Run and audit the exp006 fixed-backbone ProteinMPNN ensemble."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


EXPERIMENT_ID = "exp006_proteinmpnn_fixed_backbone_ensemble"
MODEL_SOURCE_URL = "https://github.com/dauparas/ProteinMPNN"
MODEL_REVISION = "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57"
EXPECTED_INPUT_SHA256 = (
    "2e68ae2e859b660ea2cf158be121dad6622e34ae87991eed3b286cd64b567b44"
)
EXPECTED_CHECKPOINT_SHA256 = (
    "c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd"
)
FROZEN_SEEDS = tuple(range(6000, 6010))
FROZEN_PYTHON = Path("/home/yusheng/anaconda3/bin/python")
FROZEN_TORCH_VERSION = "2.5.1"
SAMPLE_HEADER = re.compile(
    r"^>T=(?P<temperature>[^,]+), sample=(?P<sample>\d+), "
    r"score=(?P<score>[^,]+), global_score=(?P<global_score>[^,]+),"
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def command_text(command: list[str]) -> str:
    return shlex.join(command)


def git_revision(root: Path) -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()


def snapshot(command: list[str], destination: Path) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    destination.write_text(
        f"$ {command_text(command)}\n"
        f"exit_status={result.returncode}\n\n"
        f"{result.stdout}{result.stderr}"
    )


def parse_samples(fasta: Path, seed: int) -> list[dict[str, object]]:
    lines = [line.strip() for line in fasta.read_text().splitlines() if line.strip()]
    if len(lines) % 2:
        raise ValueError(f"Expected two-line FASTA records in {fasta}")
    samples: list[dict[str, object]] = []
    for header, sequence in zip(lines[0::2], lines[1::2]):
        match = SAMPLE_HEADER.match(header)
        if not match:
            continue
        samples.append(
            {
                "batch_seed": seed,
                "sample": int(match.group("sample")),
                "temperature": float(match.group("temperature")),
                "designed_chain_score": float(match.group("score")),
                "global_score": float(match.group("global_score")),
                "sequence": sequence,
                "sequence_length": len(sequence),
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "native_header": header,
            }
        )
    return samples


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--seeds", nargs="+", type=int, default=list(FROZEN_SEEDS))
    parser.add_argument("--sequences-per-batch", type=int, default=10)
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--backbone-noise", type=float, default=0.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    python = FROZEN_PYTHON
    mpnn_root = root / "external/src/ProteinMPNN"
    mpnn_script = mpnn_root / "protein_mpnn_run.py"
    weights = mpnn_root / "vanilla_model_weights"
    checkpoint = weights / "v_48_020.pt"
    input_pdb = (
        root
        / "runs/exp003_rfdiffusion_insr_binder_smoke"
        / "rfdiffusion/seed42/raw/design_ppi_42.pdb"
    )
    run_root = (
        root / "runs" / EXPERIMENT_ID / "proteinmpnn" / "fixed_backbone_cpu"
    )

    missing = [
        str(path)
        for path in (python, mpnn_script, checkpoint, input_pdb)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(f"Missing required paths: {missing}")
    torch_version = importlib.metadata.version("torch")
    if Path(sys.executable).resolve() != python.resolve():
        raise RuntimeError(f"Run this wrapper with the frozen Python: {python}")
    if torch_version != FROZEN_TORCH_VERSION:
        raise RuntimeError(
            f"Expected torch {FROZEN_TORCH_VERSION}, found {torch_version}"
        )
    if run_root.exists():
        raise FileExistsError(f"Refusing to overwrite existing run: {run_root}")
    if tuple(args.seeds) != FROZEN_SEEDS:
        raise ValueError(f"Frozen seeds are {list(FROZEN_SEEDS)}")
    if (
        args.sequences_per_batch != 10
        or args.temperature != 0.1
        or args.backbone_noise != 0.0
    ):
        raise ValueError("Frozen protocol requires 10 samples, T=0.1, and noise=0")

    input_hash = file_sha256(input_pdb)
    checkpoint_hash = file_sha256(checkpoint)
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"Input checksum mismatch: {input_hash}")
    if checkpoint_hash != EXPECTED_CHECKPOINT_SHA256:
        raise ValueError(f"Checkpoint checksum mismatch: {checkpoint_hash}")

    run_root.mkdir(parents=True)
    snapshot(["nvidia-smi"], run_root / "nvidia_smi_before.txt")
    snapshot([str(python), "-m", "pip", "freeze"], run_root / "environment.txt")

    metadata: dict[str, object] = {
        "experiment_id": EXPERIMENT_ID,
        "status": "running",
        "started_at_utc": now_utc(),
        "git_revision": git_revision(root),
        "model": {
            "source_url": MODEL_SOURCE_URL,
            "source_revision": MODEL_REVISION,
            "checkpoint": "v_48_020.pt",
            "checkpoint_sha256": checkpoint_hash,
        },
        "input": {
            "path": str(input_pdb.relative_to(root)),
            "sha256": input_hash,
            "designed_chain": "B",
            "fixed_chain": "A",
        },
        "sampling": {
            "batch_seeds": list(args.seeds),
            "sequences_per_batch": args.sequences_per_batch,
            "temperature": args.temperature,
            "backbone_noise_angstrom": args.backbone_noise,
            "batch_size": 1,
        },
        "execution": {
            "device": "cpu",
            "cuda_visible_devices": "",
            "thread_limits": {
                "OMP_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1",
            },
            "python": str(python),
            "python_version": platform.python_version(),
            "torch_version": torch_version,
            "platform": platform.platform(),
        },
        "batches": [],
    }
    write_json(run_root / "run.json", metadata)

    child_env = os.environ.copy()
    child_env.update(
        {
            "CUDA_VISIBLE_DEVICES": "",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    all_samples: list[dict[str, object]] = []
    run_started = time.monotonic()

    try:
        for seed in args.seeds:
            batch_root = run_root / f"batch_seed{seed}"
            raw_root = batch_root / "raw"
            logs_root = batch_root / "logs"
            raw_root.mkdir(parents=True)
            logs_root.mkdir()

            model_command = [
                str(python),
                str(mpnn_script),
                "--pdb_path",
                str(input_pdb),
                "--pdb_path_chains",
                "B",
                "--out_folder",
                str(raw_root),
                "--path_to_model_weights",
                str(weights),
                "--model_name",
                "v_48_020",
                "--num_seq_per_target",
                str(args.sequences_per_batch),
                "--sampling_temp",
                str(args.temperature),
                "--seed",
                str(seed),
                "--backbone_noise",
                str(args.backbone_noise),
                "--save_probs",
                "1",
                "--save_score",
                "1",
                "--batch_size",
                "1",
            ]
            timed_command = [
                "/usr/bin/time",
                "-v",
                "-o",
                str(logs_root / "resource_usage.txt"),
                *model_command,
            ]
            record: dict[str, object] = {
                "seed": seed,
                "status": "running",
                "started_at_utc": now_utc(),
                "command": command_text(model_command),
            }
            metadata["batches"].append(record)
            write_json(run_root / "run.json", metadata)

            started = time.monotonic()
            with (logs_root / "stdout.log").open("w") as stdout, (
                logs_root / "stderr.log"
            ).open("w") as stderr:
                result = subprocess.run(
                    timed_command,
                    cwd=root,
                    env=child_env,
                    stdout=stdout,
                    stderr=stderr,
                )
            record.update(
                {
                    "finished_at_utc": now_utc(),
                    "wall_seconds": round(time.monotonic() - started, 3),
                    "exit_status": result.returncode,
                }
            )
            if result.returncode:
                record["status"] = "failed"
                metadata["status"] = "failed"
                metadata["finished_at_utc"] = now_utc()
                write_json(run_root / "run.json", metadata)
                return result.returncode

            fasta = raw_root / "seqs/design_ppi_42.fa"
            probabilities = raw_root / "probs/design_ppi_42.npz"
            scores = raw_root / "scores/design_ppi_42.npz"
            if not all(path.exists() for path in (fasta, probabilities, scores)):
                raise FileNotFoundError(f"Missing outputs for seed {seed}")
            samples = parse_samples(fasta, seed)
            if len(samples) != args.sequences_per_batch:
                raise ValueError(f"Seed {seed} produced {len(samples)} samples")
            if any(row["sequence_length"] != 90 for row in samples):
                raise ValueError(f"Seed {seed} produced a non-90-residue sequence")
            all_samples.extend(samples)
            record.update(
                {
                    "status": "completed",
                    "sample_count": len(samples),
                    "fasta_sha256": file_sha256(fasta),
                    "probabilities_sha256": file_sha256(probabilities),
                    "scores_sha256": file_sha256(scores),
                }
            )
            write_json(run_root / "run.json", metadata)
    except BaseException as error:
        metadata["status"] = (
            "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
        )
        metadata["failure"] = f"{type(error).__name__}: {error}"
        metadata["finished_at_utc"] = now_utc()
        metadata["wall_seconds"] = round(time.monotonic() - run_started, 3)
        write_json(run_root / "run.json", metadata)
        raise

    fields = [
        "batch_seed",
        "sample",
        "temperature",
        "designed_chain_score",
        "global_score",
        "sequence",
        "sequence_length",
        "sequence_sha256",
        "native_header",
    ]
    predictions = run_root / "sequences.csv"
    with predictions.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_samples)

    metadata.update(
        {
            "status": "completed",
            "finished_at_utc": now_utc(),
            "wall_seconds": round(time.monotonic() - run_started, 3),
            "sample_count": len(all_samples),
            "unique_sequence_count": len({row["sequence"] for row in all_samples}),
            "aggregate_predictions": "sequences.csv",
            "aggregate_predictions_sha256": file_sha256(predictions),
        }
    )
    write_json(run_root / "run.json", metadata)
    snapshot(["nvidia-smi"], run_root / "nvidia_smi_after.txt")
    print(
        f"Generated {metadata['sample_count']} samples "
        f"({metadata['unique_sequence_count']} unique) in "
        f"{metadata['wall_seconds']} s."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
