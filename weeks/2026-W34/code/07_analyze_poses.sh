#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 07_analyze_poses.sh [--check|--execute] [--allow-rcsb-download]

Download the 16 frozen crystal references when explicitly permitted, align
predicted and experimental protein pockets, calculate symmetry-corrected ligand
heavy-atom RMSD, and build the paired Boltz-2/FlashBind table.
EOF
}

BAP_EXECUTE=0
allow_rcsb=0
for argument in "$@"; do
  case "$argument" in
    --check) BAP_EXECUTE=0 ;;
    --execute) BAP_EXECUTE=1 ;;
    --allow-rcsb-download) allow_rcsb=1 ;;
    -h|--help) usage; exit 0 ;;
    *) bap_die "Unknown argument: $argument" ;;
  esac
done

raw_pose="$BAP_REPO_ROOT/data/raw/exp008_boltz2_crystal_pose"
boltz_pose_run="$BAP_RUN_ROOT/exp008_boltz2_crystal_pose/seed42"
flash_pose_run="$BAP_RUN_ROOT/exp010_flashbind_crystal_pose/released_poses"
paired_report="$BAP_REPO_ROOT/reports/exp010_flashbind_crystal_pose"

if [[ ! -f "$raw_pose/download_manifest.json" ]]; then
  ((BAP_EXECUTE && allow_rcsb)) || bap_die \
    "Crystal references are missing; use --execute --allow-rcsb-download"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_FLASHBIND_ENV/bin/python" \
    "$BAP_REPO_ROOT/scripts/fetch_fepplus4_crystal_references.py"
fi

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src:$BAP_REPO_ROOT/scripts" \
  "$BAP_FLASHBIND_ENV/bin/python" "$BAP_REPO_ROOT/scripts/analyze_boltz2_crystal_pose.py" \
  --manifest "$BAP_POSE_MANIFEST" --parent-manifest "$BAP_MANIFEST" \
  --raw-root "$raw_pose" --output-dir "$boltz_pose_run"

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src:$BAP_REPO_ROOT/scripts" \
  "$BAP_FLASHBIND_ENV/bin/python" "$BAP_REPO_ROOT/scripts/analyze_flashbind_crystal_pose.py" \
  --manifest "$BAP_POSE_MANIFEST" --parent-manifest "$BAP_MANIFEST" \
  --raw-root "$raw_pose" --flashbind-root "$BAP_FLASHBIND_DATA" \
  --output-dir "$flash_pose_run"

bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
  "$BAP_FLASHBIND_ENV/bin/python" "$BAP_REPO_ROOT/scripts/compare_pose_models.py" \
  --flashbind "$flash_pose_run/pose_metrics.csv" \
  --boltz2 "$boltz_pose_run/pose_metrics.csv" --output-dir "$paired_report"

if ((BAP_EXECUTE)); then
  install -m 0644 "$boltz_pose_run/pose_metrics.csv" \
    "$BAP_REPO_ROOT/reports/exp008_boltz2_crystal_pose/pose_metrics.csv"
  install -m 0644 "$flash_pose_run/pose_metrics.csv" \
    "$paired_report/flashbind_pose_metrics.csv"
else
  bap_require_file "$boltz_pose_run/pose_metrics.csv"
  bap_require_file "$flash_pose_run/pose_metrics.csv"
  bap_require_file "$paired_report/paired_pose_rmsd.csv"
  bap_require_file "$paired_report/paired_pose_comparison.json"
fi

bap_log "Crystal-pose analyses and paired RMSD comparison are ready."
