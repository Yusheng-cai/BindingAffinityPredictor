#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 04_run_inference.sh --model MODEL [--target TARGET] [--execute]

MODEL is boltz2, nesso1, flashbind, or all.
TARGET is cdk2, tyk2, jnk1, p38, or all (default).

Without --execute, the script prints the resolved commands. Completed run
records are preserved and skipped. Boltz-2 uses the approved MSA-1024 variant.
EOF
}

model=""
target_selection="all"
while (($#)); do
  case "$1" in
    --model=*) model="${1#*=}"; shift ;;
    --model) [[ $# -ge 2 ]] || bap_die "--model requires a value"; model="$2"; shift 2 ;;
    --target=*) target_selection="${1#*=}"; shift ;;
    --target) [[ $# -ge 2 ]] || bap_die "--target requires a value"; target_selection="$2"; shift 2 ;;
    --execute) BAP_EXECUTE=1; shift ;;
    --check) BAP_EXECUTE=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) bap_die "Unknown argument: $1" ;;
  esac
done
[[ -n "$model" ]] || { usage; bap_die "--model is required"; }
[[ "$model" =~ ^(boltz2|nesso1|flashbind|all)$ ]] || bap_die "Invalid model: $model"
[[ "$target_selection" =~ ^(cdk2|tyk2|jnk1|p38|all)$ ]] || \
  bap_die "Invalid target: $target_selection"

selected_targets="$BAP_TARGETS"
[[ "$target_selection" == all ]] || selected_targets="$target_selection"
runner="$BAP_REPO_ROOT/scripts/run_with_provenance.py"
bap_require_file "$runner"
bap_require_file "$BAP_MANIFEST"

if [[ "$model" == nesso1 || "$model" == all ]]; then
  bap_require_executable "$BAP_NESSO_ENV/bin/nesso"
  for target in $selected_targets; do bap_require_dir "$BAP_EXP007/inputs/nesso1/$target"; done
fi
if [[ "$model" == boltz2 || "$model" == all ]]; then
  bap_require_executable "$BAP_BOLTZ_ENV/bin/boltz"
  bap_require_file "$BAP_BOLTZ_STRUCTURE_CHECKPOINT"
  bap_require_file "$BAP_BOLTZ_AFFINITY_CHECKPOINT"
  for target in $selected_targets; do bap_require_dir "$BAP_EXP007/inputs/boltz2/$target"; done
fi
if [[ "$model" == flashbind || "$model" == all ]]; then
  bap_require_file "$BAP_FLASHBIND_SOURCE/scripts/predict.py"
  bap_require_file "$BAP_FLASHBIND_DATA/id.json"
  bap_require_file "$BAP_FLASHBIND_CHECKPOINT_1"
  bap_require_file "$BAP_FLASHBIND_CHECKPOINT_2"
fi

run_nesso_target() {
  local target="$1"
  local run_dir="$BAP_EXP007_RUN/nesso1/seed${BAP_SEED}/$target"
  if ((BAP_EXECUTE)) && bap_skip_complete_run "$run_dir"; then return; fi
  bap_run env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" "$runner" \
    --run-dir "$run_dir" --model nesso1 --cwd "$BAP_REPO_ROOT" -- \
    env -u PYTHONPATH "$BAP_NESSO_ENV/bin/nesso" predict \
    "$BAP_EXP007/inputs/nesso1/$target" \
    --out_dir "$run_dir/raw" --cache "$BAP_NESSO_CACHE" \
    --model_revision v1.0.0 --accelerator gpu --devices 1 \
    --precision bf16-mixed --recycling_steps 5 --num_workers 1 \
    --require_affinity --refine_protein_inference \
    --refine_protein_cutoff 22.0 --refine_protein_tokens_budget 256 \
    --no_kernels --seed "$BAP_SEED"
}

run_boltz_target() {
  local target="$1"
  local run_dir="$BAP_EXP007_RUN/boltz2_msa1024/seed${BAP_SEED}/$target"
  if ((BAP_EXECUTE)) && bap_skip_complete_run "$run_dir"; then return; fi
  bap_run env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/python" "$runner" \
    --run-dir "$run_dir" --model boltz2_msa1024 --cwd "$BAP_REPO_ROOT" -- \
    env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/boltz" predict \
    "$BAP_EXP007/inputs/boltz2/$target" \
    --out_dir "$run_dir/raw" --cache "$BAP_BOLTZ_CACHE" --model boltz2 \
    --checkpoint "$BAP_BOLTZ_STRUCTURE_CHECKPOINT" \
    --affinity_checkpoint "$BAP_BOLTZ_AFFINITY_CHECKPOINT" \
    --accelerator gpu --devices 1 --recycling_steps 3 --sampling_steps 200 \
    --diffusion_samples 1 --max_parallel_samples 1 --step_scale 1.5 \
    --output_format mmcif --num_workers 1 --seed "$BAP_SEED" \
    --sampling_steps_affinity 200 --diffusion_samples_affinity 5 \
    --max_msa_seqs 8192 --subsample_msa --num_subsampled_msa 1024 \
    --no_kernels
}

run_flashbind() {
  local run_dir="$BAP_EXP009_RUN/flashbind_released_poses/seed${BAP_SEED}/full87"
  if ((BAP_EXECUTE)) && bap_skip_complete_run "$run_dir"; then return; fi
  bap_run env -u PYTHONPATH "$BAP_FLASHBIND_ENV/bin/python" "$runner" \
    --run-dir "$run_dir" --model flashbind --cwd "$BAP_REPO_ROOT" -- \
    env PYTHONPATH="$BAP_FLASHBIND_SOURCE/src" \
    "$BAP_FLASHBIND_ENV/bin/python" "$BAP_FLASHBIND_SOURCE/scripts/predict.py" \
    --data "$BAP_FLASHBIND_DATA/id.json" --task value \
    --structure "$BAP_FLASHBIND_DATA/pdb" --structure_type pdb \
    --ligand "$BAP_FLASHBIND_DATA/ligand_sdf.lmdb" --ligand_type sdf \
    --protein_repr "$BAP_FLASHBIND_DATA/repr/esm3.lmdb" \
    --ligand_repr "$BAP_FLASHBIND_DATA/repr/torchdrug.lmdb" \
    --distance_threshold 20.0 --out_dir "$run_dir/raw" \
    --devices 1 --accelerator gpu --num_workers 0 --seed "$BAP_SEED" \
    --affinity_checkpoint "$BAP_FLASHBIND_CHECKPOINT_1" "$BAP_FLASHBIND_CHECKPOINT_2"
}

if [[ "$model" == nesso1 || "$model" == all ]]; then
  for target in $selected_targets; do run_nesso_target "$target"; done
fi
if [[ "$model" == boltz2 || "$model" == all ]]; then
  for target in $selected_targets; do run_boltz_target "$target"; done
fi
if [[ "$model" == flashbind || "$model" == all ]]; then
  [[ "$target_selection" == all ]] || bap_die "FlashBind is run as one released 87-record archive"
  run_flashbind
fi

bap_log "Requested inference commands are complete or already preserved."
