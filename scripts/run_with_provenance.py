#!/usr/bin/env python3
"""Run a model command while preserving logs, timing, Git state, and GPU usage."""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path


def gpu_snapshot() -> dict:
    command = [
        "nvidia-smi",
        "--query-gpu=index,name,memory.total,memory.used,utilization.gpu,temperature.gpu",
        "--format=csv,noheader,nounits",
    ]
    try:
        output = subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).strip()
        fields = [field.strip() for field in output.splitlines()[0].split(",")]
        return {
            "index": int(fields[0]),
            "name": fields[1],
            "memory_total_mib": int(fields[2]),
            "memory_used_mib": int(fields[3]),
            "utilization_percent": int(fields[4]),
            "temperature_c": int(fields[5]),
        }
    except Exception as error:  # provenance capture must not abort inference
        return {"error": repr(error)}


def stream(pipe, destination, log_handle) -> None:
    for line in iter(pipe.readline, ""):
        destination.write(line)
        destination.flush()
        log_handle.write(line)
        log_handle.flush()
    pipe.close()


def git_metadata(cwd: Path) -> dict:
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=cwd, text=True).strip()
        status = subprocess.check_output(["git", "status", "--short"], cwd=cwd, text=True)
        return {"revision": revision, "working_tree_dirty": bool(status.strip())}
    except Exception as error:
        return {"error": repr(error)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--cwd", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command and args.command[0] == "--" else args.command
    if not command:
        raise ValueError("a command is required after --")

    logs = args.run_dir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    before = gpu_snapshot()
    started = datetime.now(timezone.utc)
    start_monotonic = time.monotonic()
    peak_gpu = before.get("memory_used_mib")
    samples = []

    with (logs / "stdout.log").open("w", encoding="utf-8") as stdout_log, (
        logs / "stderr.log"
    ).open("w", encoding="utf-8") as stderr_log:
        process = subprocess.Popen(
            command,
            cwd=args.cwd,
            env=os.environ.copy(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        stdout_thread = threading.Thread(
            target=stream, args=(process.stdout, sys.stdout, stdout_log), daemon=True
        )
        stderr_thread = threading.Thread(
            target=stream, args=(process.stderr, sys.stderr, stderr_log), daemon=True
        )
        stdout_thread.start()
        stderr_thread.start()
        while process.poll() is None:
            snapshot = gpu_snapshot()
            snapshot["elapsed_seconds"] = time.monotonic() - start_monotonic
            samples.append(snapshot)
            used = snapshot.get("memory_used_mib")
            if isinstance(used, int):
                peak_gpu = used if peak_gpu is None else max(peak_gpu, used)
            time.sleep(0.5)
        stdout_thread.join()
        stderr_thread.join()
        return_code = process.returncode

    ended = datetime.now(timezone.utc)
    after = gpu_snapshot()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    baseline = before.get("memory_used_mib")
    record = {
        "schema_version": 1,
        "model": args.model,
        "status": "complete" if return_code == 0 else "failed",
        "return_code": return_code,
        "command": command,
        "cwd": str(args.cwd.resolve()),
        "started_utc": started.isoformat(),
        "ended_utc": ended.isoformat(),
        "wall_time_seconds": time.monotonic() - start_monotonic,
        "max_host_rss_kib": usage.ru_maxrss,
        "gpu_before": before,
        "gpu_after": after,
        "peak_total_gpu_memory_mib": peak_gpu,
        "peak_gpu_memory_above_baseline_mib": (
            peak_gpu - baseline if isinstance(peak_gpu, int) and isinstance(baseline, int) else None
        ),
        "gpu_samples": len(samples),
        "python": sys.version,
        "platform": platform.platform(),
        "project_git": git_metadata(args.cwd),
        "selected_environment": {
            key: os.environ.get(key)
            for key in (
                "CUDA_VISIBLE_DEVICES",
                "BOLTZ_CACHE",
                "NESSO_CACHE",
                "NESSO_MODEL_REVISION",
                "FLASHBIND_CACHE",
                "FLASHBIND_SOURCE_REVISION",
                "FLASHBIND_DATA_REVISION",
                "FLASHBIND_MODEL_REVISION",
            )
            if os.environ.get(key) is not None
        },
    }
    (args.run_dir / "run.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    if return_code != 0:
        raise SystemExit(return_code)


if __name__ == "__main__":
    main()
