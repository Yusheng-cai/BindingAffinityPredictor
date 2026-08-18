#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 06_analyze_affinity.sh [--check|--execute]

Calculate the predefined per-target and compound-weighted Pearson, Spearman,
Kendall, pairwise-error, centered-MAE, and bootstrap statistics.
EOF
}

if ! bap_parse_execute_only "${1:---check}"; then
  usage
  exit 0
fi

nesso_predictions="$BAP_EXP007_RUN/nesso1/seed${BAP_SEED}/predictions.csv"
boltz_predictions="$BAP_EXP007_RUN/boltz2_msa1024/seed${BAP_SEED}/predictions.csv"
flashbind_predictions="$BAP_EXP009_RUN/flashbind_released_poses/seed${BAP_SEED}/full87/predictions.csv"
analysis="$BAP_REPO_ROOT/scripts/analyze_fepplus4.py"
report007="$BAP_REPO_ROOT/reports/exp007_fepplus4_boltz2_nesso1"
report009="$BAP_REPO_ROOT/reports/exp009_flashbind_fepplus4_released_poses"

for path in "$nesso_predictions" "$boltz_predictions" "$flashbind_predictions"; do
  bap_require_file "$path"
done

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
  "$BAP_NESSO_ENV/bin/python" "$analysis" --manifest "$BAP_MANIFEST" \
  --predictions "$nesso_predictions" \
  --output "$report007/nesso1_metrics.json" \
  --bootstrap-iterations 2000 --bootstrap-seed 20260817

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
  "$BAP_NESSO_ENV/bin/python" "$analysis" --manifest "$BAP_MANIFEST" \
  --predictions "$nesso_predictions" --predictions "$boltz_predictions" \
  --output "$report007/boltz2_msa1024_vs_nesso1_metrics.json" \
  --bootstrap-iterations 2000 --bootstrap-seed 20260817

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
  "$BAP_NESSO_ENV/bin/python" "$analysis" --manifest "$BAP_MANIFEST" \
  --predictions "$flashbind_predictions" \
  --output "$report009/flashbind_metrics.json" \
  --bootstrap-iterations 2000 --bootstrap-seed 20260817

bap_log "Affinity benchmark metrics are ready."
