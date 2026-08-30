#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 03_analyze_runs_n_poses.sh [--check|--execute]

Score Nesso-1 distograms and Boltz-2 coordinates against the same 50 selected
experimental complexes, then build the 49-system paired/familiarity comparison.
Completed analysis records are preserved and skipped.
EOF
}

case "${1:---check}" in
  --check) BAP_EXECUTE=0 ;;
  --execute) BAP_EXECUTE=1 ;;
  -h|--help) usage; exit 0 ;;
  *) bap_die "Expected --check or --execute" ;;
esac

runner="$BAP_REPO_ROOT/scripts/run_with_provenance.py"
nesso_metrics="$BAP_RUN_ROOT/analysis/nesso1_full50/results/per_system_metrics.json"
boltz_metrics="$BAP_RUN_ROOT/analysis/boltz2_full50/results/per_system_metrics.json"
nesso_run="$BAP_RUN_ROOT/nesso1/seed42/full50/run.json"
boltz_run="$BAP_RUN_ROOT/boltz2_msa1024/seed42/full50/run.json"

bap_require_file "$runner"
bap_require_file "$nesso_run"
bap_require_file "$boltz_run"
bap_require_dir "$BAP_GROUND_TRUTH_ROOT/ground_truth"

run_nesso_analysis() {
  local run_dir="$BAP_RUN_ROOT/analysis/nesso1_full50"
  if bap_complete_run "$run_dir"; then
    bap_log "Preserving completed Nesso analysis: $run_dir/run.json"
    return
  fi
  ((BAP_EXECUTE)) || bap_die "Nesso analysis is missing: $run_dir"
  bap_refuse_incomplete_run "$run_dir"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$runner" \
    --run-dir "$run_dir" --model nesso1_structure_scoring --cwd "$BAP_REPO_ROOT" \
    --experiment-config "$BAP_EXPERIMENT/experiment.yaml" \
    --benchmark-config "$BAP_BENCHMARK" --model-config "$BAP_REPO_ROOT/configs/models/nesso1.yaml" \
    --manifest "$BAP_MANIFEST" -- \
    env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/analyze_nesso_rnp_distograms.py" \
    --manifest "$BAP_MANIFEST" --ground-truth-root "$BAP_GROUND_TRUTH_ROOT/ground_truth" \
    --predictions-root "$BAP_RUN_ROOT/nesso1/seed42/full50/raw/predictions" \
    --output-dir "$run_dir/results" --experiment-id "$BAP_EXP_ID" \
    --familiarity-source nesso1_sep2021 --symmetry-policy metric_specific
}

run_boltz_analysis() {
  local run_dir="$BAP_RUN_ROOT/analysis/boltz2_full50"
  if bap_complete_run "$run_dir"; then
    bap_log "Preserving completed Boltz analysis: $run_dir/run.json"
    return
  fi
  ((BAP_EXECUTE)) || bap_die "Boltz analysis is missing: $run_dir"
  bap_refuse_incomplete_run "$run_dir"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$runner" \
    --run-dir "$run_dir" --model boltz2_structure_scoring --cwd "$BAP_REPO_ROOT" \
    --experiment-config "$BAP_EXPERIMENT/experiment.yaml" \
    --benchmark-config "$BAP_BENCHMARK" --model-config "$BAP_REPO_ROOT/configs/models/boltz2.yaml" \
    --manifest "$BAP_MANIFEST" -- \
    env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/analyze_boltz2_rnp_structures.py" \
    --manifest "$BAP_MANIFEST" --ground-truth-root "$BAP_GROUND_TRUTH_ROOT/ground_truth" \
    --prediction-root "$BAP_RUN_ROOT/boltz2_msa1024/seed42/full50/raw/boltz_results_all" \
    --output-dir "$run_dir/results"
}

run_paired_analysis() {
  local run_dir="$BAP_RUN_ROOT/analysis/paired_comparison"
  if bap_complete_run "$run_dir"; then
    bap_log "Preserving completed paired analysis: $run_dir/run.json"
    return
  fi
  ((BAP_EXECUTE)) || bap_die "Paired analysis is missing: $run_dir"
  bap_refuse_incomplete_run "$run_dir"
  python3 "$runner" \
    --run-dir "$run_dir" --model paired_structure_analysis --cwd "$BAP_REPO_ROOT" \
    --experiment-config "$BAP_EXPERIMENT/experiment.yaml" \
    --benchmark-config "$BAP_BENCHMARK" --manifest "$BAP_MANIFEST" -- \
    python3 "$BAP_REPO_ROOT/scripts/compare_rnp_paired_models.py" \
    --manifest "$BAP_MANIFEST" --nesso-metrics "$nesso_metrics" \
    --boltz-metrics "$boltz_metrics" --nesso-run "$nesso_run" --boltz-run "$boltz_run" \
    --output-dir "$BAP_RESULT_ROOT" --bootstrap-iterations 2000 --bootstrap-seed 2026082803
}

run_nesso_analysis
run_boltz_analysis
run_paired_analysis
bap_require_file "$BAP_RESULT_ROOT/summary.json"
bap_require_file "$BAP_RESULT_ROOT/paired_metrics.csv"
bap_require_file "$BAP_RESULT_ROOT/familiarity_binned_metrics.csv"
bap_log "Validated the compact paired Runs N' Poses result."
