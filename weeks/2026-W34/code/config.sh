#!/usr/bin/env bash
# Shared paths and pinned scientific settings for the 2026-W34 workflow.

# This file is sourced by the numbered workflow scripts. Machine-specific root
# paths can be overridden in the environment. Scientific settings remain fixed
# here; changing them defines a new protocol rather than reproducing this week.

BAP_CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAP_WEEK_DIR="$(cd "$BAP_CODE_DIR/.." && pwd)"
BAP_REPO_ROOT="$(cd "$BAP_WEEK_DIR/../.." && pwd)"

: "${BAP_SOURCE_ROOT:=/home/yusheng/source}"
: "${BAP_PROGRAM_ROOT:=/home/yusheng/programs}"
: "${BAP_CACHE_ROOT:=$BAP_REPO_ROOT/cache}"
: "${BAP_RUN_ROOT:=$BAP_REPO_ROOT/runs}"
: "${BAP_CONDA_EXE:=/home/yusheng/anaconda3/bin/conda}"

BAP_SEED=42
BAP_TARGETS="cdk2 tyk2 jnk1 p38"

BAP_BOLTZ_SOURCE="$BAP_SOURCE_ROOT/boltz"
BAP_NESSO_SOURCE="$BAP_SOURCE_ROOT/nesso"
BAP_FLASHBIND_SOURCE="$BAP_SOURCE_ROOT/FlashBind"
BAP_OPENFF_SOURCE="$BAP_SOURCE_ROOT/protein-ligand-benchmark"

BAP_BOLTZ_ENV="$BAP_PROGRAM_ROOT/boltz-2.2.1"
BAP_NESSO_ENV="$BAP_PROGRAM_ROOT/nesso-1.0.0"
BAP_FLASHBIND_ENV="$BAP_PROGRAM_ROOT/flashbind-f161268"

BAP_BOLTZ_CACHE="$BAP_CACHE_ROOT/boltz"
BAP_NESSO_CACHE="$BAP_CACHE_ROOT/nesso"
BAP_FLASHBIND_CACHE="$BAP_CACHE_ROOT/flashbind"

BAP_BOLTZ_REVISION="b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc"
BAP_NESSO_REVISION="f0156e9a22326448684bae09ee96f73415902dcd"
BAP_FLASHBIND_REVISION="f161268176237ab6ce5757031a8c1b93937e0d37"
BAP_OPENFF_REVISION="da7c3372256446222e424368be38ef3d2b55a67b"
BAP_NESSO_MODEL_REVISION="1896c84c7186c506c7efd79051480809d51098bf"
BAP_NESSO_ESM_REVISION="08e4846e537177426273712802403f7ba8261b6c"
BAP_FLASHBIND_MODEL_REVISION="fa6362c85b3350109af4634e3e2c364644338b3c"
BAP_FLASHBIND_DATA_REVISION="50b1511e080236d80f9b8ff1e4d0cd38bb2480b9"

BAP_MANIFEST="$BAP_REPO_ROOT/data/manifests/fepplus4_87.json"
BAP_POSE_MANIFEST="$BAP_REPO_ROOT/data/manifests/fepplus4_crystal16.json"
BAP_EXP007="$BAP_REPO_ROOT/experiments/exp007_fepplus4_boltz2_nesso1"
BAP_EXP009="$BAP_REPO_ROOT/experiments/exp009_flashbind_fepplus4_released_poses"
BAP_EXP007_RUN="$BAP_RUN_ROOT/exp007_fepplus4_boltz2_nesso1"
BAP_EXP009_RUN="$BAP_RUN_ROOT/exp009_flashbind_fepplus4_released_poses"

BAP_BOLTZ_STRUCTURE_CHECKPOINT="$BAP_BOLTZ_ENV/checkpoints/boltz2_conf.ckpt"
BAP_BOLTZ_AFFINITY_CHECKPOINT="$BAP_BOLTZ_ENV/checkpoints/boltz2_aff.ckpt"
BAP_FLASHBIND_DATA="$BAP_FLASHBIND_CACHE/datasets/fep4/fep4"
BAP_FLASHBIND_CHECKPOINT_1="$BAP_FLASHBIND_CACHE/checkpoints/value_1.ckpt"
BAP_FLASHBIND_CHECKPOINT_2="$BAP_FLASHBIND_CACHE/checkpoints/value_2.ckpt"

export BAP_CODE_DIR BAP_WEEK_DIR BAP_REPO_ROOT
export BAP_SOURCE_ROOT BAP_PROGRAM_ROOT BAP_CACHE_ROOT BAP_RUN_ROOT
export BAP_SEED BAP_TARGETS BAP_CONDA_EXE
export BAP_BOLTZ_CACHE BAP_NESSO_CACHE BAP_FLASHBIND_CACHE
