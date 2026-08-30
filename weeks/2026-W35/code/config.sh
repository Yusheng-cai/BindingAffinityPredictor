#!/usr/bin/env bash
# Shared paths and frozen settings for the 2026-W35 reproduction workflow.

BAP_CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAP_WEEK_DIR="$(cd "$BAP_CODE_DIR/.." && pwd)"
BAP_REPO_ROOT="$(cd "$BAP_WEEK_DIR/../.." && pwd)"

: "${BAP_PROGRAM_ROOT:=/home/yusheng/programs}"
: "${BAP_CACHE_ROOT:=$BAP_REPO_ROOT/cache}"
: "${BAP_DATA_PYTHON:=/home/yusheng/anaconda3/bin/python}"

BAP_EXP_ID="exp015_boltz2_nesso1_rnp_postcutoff50"
BAP_EXPERIMENT="$BAP_REPO_ROOT/experiments/$BAP_EXP_ID"
BAP_BENCHMARK="$BAP_REPO_ROOT/configs/benchmarks/rnp_boltz2_nesso1_postcutoff50.yaml"
BAP_MANIFEST="$BAP_REPO_ROOT/data/manifests/rnp_boltz2_nesso1_postcutoff50.json"
BAP_GROUND_TRUTH_ARCHIVE="$BAP_REPO_ROOT/data/raw/runs_n_poses_metadata/ground_truth.tar.gz"
BAP_GROUND_TRUTH_ROOT="$BAP_REPO_ROOT/data/raw/runs_n_poses_ground_truth_exp015"
BAP_RUN_ROOT="$BAP_REPO_ROOT/runs/$BAP_EXP_ID"
BAP_RESULT_ROOT="$BAP_REPO_ROOT/reports/$BAP_EXP_ID/results"

BAP_NESSO_ENV="$BAP_PROGRAM_ROOT/nesso-1.0.0"
BAP_BOLTZ_ENV="$BAP_PROGRAM_ROOT/boltz-2.2.1"
BAP_NESSO_CACHE="$BAP_CACHE_ROOT/nesso"
BAP_BOLTZ_CACHE="$BAP_CACHE_ROOT/boltz"
BAP_MSA_ROOT="$BAP_RUN_ROOT/boltz2/msas"
BAP_SEED=42

export BAP_CODE_DIR BAP_WEEK_DIR BAP_REPO_ROOT BAP_PROGRAM_ROOT BAP_CACHE_ROOT
export BAP_DATA_PYTHON
export BAP_EXP_ID BAP_EXPERIMENT BAP_BENCHMARK BAP_MANIFEST
export BAP_GROUND_TRUTH_ARCHIVE BAP_GROUND_TRUTH_ROOT BAP_RUN_ROOT BAP_RESULT_ROOT
export BAP_NESSO_ENV BAP_BOLTZ_ENV BAP_NESSO_CACHE BAP_BOLTZ_CACHE BAP_MSA_ROOT BAP_SEED
