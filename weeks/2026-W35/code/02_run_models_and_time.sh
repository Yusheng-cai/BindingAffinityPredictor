#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 02_run_models_and_time.sh [--model nesso1|boltz2|all] [--check|--execute]

The model commands are frozen in the exp015 gate/full scripts. Each command is
wrapped by scripts/run_with_provenance.py, which records wall time, peak memory,
hardware, Git revision, resolved command, settings checksums, and failure status.
Existing complete run directories are preserved and never overwritten.
EOF
}

model=all
while (($#)); do
  case "$1" in
    --model=*) model="${1#*=}"; shift ;;
    --model) [[ $# -ge 2 ]] || bap_die "--model requires a value"; model="$2"; shift 2 ;;
    --check) BAP_EXECUTE=0; shift ;;
    --execute) BAP_EXECUTE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) bap_die "Unknown argument: $1" ;;
  esac
done
[[ "$model" =~ ^(nesso1|boltz2|all)$ ]] || bap_die "Invalid model: $model"

bap_require_file "$BAP_EXPERIMENT/run_gate.sh"
bap_require_file "$BAP_EXPERIMENT/run_full.sh"
bap_require_file "$BAP_REPO_ROOT/scripts/run_with_provenance.py"

models=(nesso1 boltz2)
[[ "$model" == all ]] || models=("$model")
for current in "${models[@]}"; do
  if [[ "$current" == nesso1 ]]; then
    gate_dir="$BAP_RUN_ROOT/nesso1/seed42/gate5"
    full_dir="$BAP_RUN_ROOT/nesso1/seed42/full50"
  else
    gate_dir="$BAP_RUN_ROOT/boltz2_msa1024/seed42/gate5"
    full_dir="$BAP_RUN_ROOT/boltz2_msa1024/seed42/full50"
  fi

  if bap_complete_run "$gate_dir"; then
    bap_log "Preserving completed $current gate: $gate_dir/run.json"
  elif ((BAP_EXECUTE)); then
    bap_refuse_incomplete_run "$gate_dir"
    bash "$BAP_EXPERIMENT/run_gate.sh" "$current"
  else
    bap_die "$current gate is not complete: $gate_dir"
  fi

  if bap_complete_run "$full_dir"; then
    bap_log "Preserving completed $current full run: $full_dir/run.json"
  elif ((BAP_EXECUTE)); then
    bap_refuse_incomplete_run "$full_dir"
    bash "$BAP_EXPERIMENT/run_full.sh" "$current"
  else
    bap_die "$current full run is not complete: $full_dir"
  fi

  "$BAP_NESSO_ENV/bin/python" - "$full_dir/run.json" <<'PY'
import json
import sys

record = json.load(open(sys.argv[1], encoding="utf-8"))
print(
    f"{record['model']}: status={record['status']}, "
    f"wall={record['wall_time_seconds']:.3f} s, "
    f"peak GPU={record.get('peak_total_gpu_memory_mib')} MiB"
)
PY
done

bap_log "The full-run wall times above are the raw timing measurements used in the report."
