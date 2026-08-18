#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 05_collect_outputs.sh [--check|--execute]

Convert unmodified model-native JSON outputs into the common prediction-table
schema. The native outputs remain in runs/*/raw/.
EOF
}

if ! bap_parse_execute_only "${1:---check}"; then
  usage
  exit 0
fi

nesso_run="$BAP_EXP007_RUN/nesso1/seed${BAP_SEED}"
boltz_run="$BAP_EXP007_RUN/boltz2_msa1024/seed${BAP_SEED}"
flashbind_run="$BAP_EXP009_RUN/flashbind_released_poses/seed${BAP_SEED}/full87"

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
  "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/collect_affinity_outputs.py" \
  --manifest "$BAP_MANIFEST" --model nesso1 --raw-root "$nesso_run" \
  --output "$nesso_run/predictions.csv"

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
  "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/collect_affinity_outputs.py" \
  --manifest "$BAP_MANIFEST" --model boltz2 --raw-root "$boltz_run" \
  --output "$boltz_run/predictions.csv"

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
  "$BAP_FLASHBIND_ENV/bin/python" "$BAP_REPO_ROOT/scripts/collect_flashbind_outputs.py" \
  --manifest "$BAP_MANIFEST" \
  --predictions-json "$flashbind_run/raw/affinity_predictions_ensemble.json" \
  --output "$flashbind_run/predictions.csv"

if ((BAP_EXECUTE == 0)); then
  bap_require_file "$nesso_run/predictions.csv"
  bap_require_file "$boltz_run/predictions.csv"
  bap_require_file "$flashbind_run/predictions.csv"
fi

bap_log "Canonical prediction tables are ready."
