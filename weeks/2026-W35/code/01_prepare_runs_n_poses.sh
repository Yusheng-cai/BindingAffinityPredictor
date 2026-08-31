#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 01_prepare_runs_n_poses.sh [--check|--execute] [--allow-msa-server]

Run 00_setup.sh first. Stage 00 installs the pinned software and downloads the
minimal checksum-verified Runs N' Poses record needed by this stage.

--check             Validate the frozen manifest, selected references, MSAs,
                    and 50 full/5 gate model inputs without changing them.
--execute           Reproduce the familiarity table and 50-system selection,
                    then recreate missing references and model-native inputs.
--allow-msa-server  Permit submission of protein sequences to the public
                    ColabFold MMseqs2 server if an MSA is missing.
EOF
}

allow_msa=0
for argument in "$@"; do
  case "$argument" in
    --check) BAP_EXECUTE=0 ;;
    --execute) BAP_EXECUTE=1 ;;
    --allow-msa-server) allow_msa=1 ;;
    -h|--help) usage; exit 0 ;;
    *) bap_die "Unknown argument: $argument" ;;
  esac
done

bap_require_file "$BAP_MANIFEST"
bap_require_file "$BAP_BENCHMARK"
bap_require_file "$BAP_EXPERIMENT/experiment.yaml"
bap_require_executable "$BAP_NESSO_ENV/bin/python"
bap_require_executable "$BAP_DATA_PYTHON"
metadata_root="$BAP_REPO_ROOT/data/raw/runs_n_poses_metadata"
bap_require_file "$metadata_root/all_similarity_scores.parquet"
bap_require_file "$metadata_root/annotations.csv"
bap_require_file "$metadata_root/inputs.json"
printf '%s  %s\n' \
  ce771a2439a91210b7e4ebd09a729fb75d95b6c8265f4f38bae95913f4016804 \
  "$metadata_root/all_similarity_scores.parquet" | sha256sum --check -
printf '%s  %s\n' \
  259aa0a8f5ea6008d4036886f1ad5ae255689c223eb839fa1a2e7cb54adae609 \
  "$metadata_root/annotations.csv" | sha256sum --check -
printf '%s  %s\n' \
  85ded25555c1efba0aeb7dc90df2fbf24d9550f390ed3e8ac9b32cf11a0ede4d \
  "$metadata_root/inputs.json" | sha256sum --check -
printf '%s  %s\n' \
  ac98b0d4a0a632d4643beb8d999de59b2c8561c31538ab6565e95b145312e5ed \
  "$BAP_MANIFEST" | sha256sum --check -

if ((BAP_EXECUTE)); then
  reproduction_root="$BAP_REPO_ROOT/data/processed/runs_n_poses/reproduction_w35"
  reproduced_similarity="$reproduction_root/boltz2_2023_similarity.csv"
  reproduced_manifest="$reproduction_root/rnp_boltz2_nesso1_postcutoff50.json"
  mkdir -p "$reproduction_root"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_DATA_PYTHON" "$BAP_REPO_ROOT/scripts/extract_rnp_cutoff_similarity.py" \
    --input "$metadata_root/all_similarity_scores.parquet" \
    --output "$reproduced_similarity" --cutoff 2023-06-01
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/select_rnp_boltz2_nesso1_postcutoff.py" \
    --annotations "$metadata_root/annotations.csv" --inputs "$metadata_root/inputs.json" \
    --similarity "$reproduced_similarity" --output "$reproduced_manifest"
  reproduced_sha="$(sha256sum "$reproduced_manifest" | awk '{print $1}')"
  [[ "$reproduced_sha" == ac98b0d4a0a632d4643beb8d999de59b2c8561c31538ab6565e95b145312e5ed ]] || \
    bap_die "Reproduced selection differs from the frozen manifest: $reproduced_sha"
  bap_log "Reproduced the June-2023 familiarity table and byte-identical 50-system manifest."
fi

if [[ ! -d "$BAP_GROUND_TRUTH_ROOT/ground_truth" ]]; then
  ((BAP_EXECUTE)) || bap_die "Selected ground truth is missing: $BAP_GROUND_TRUTH_ROOT/ground_truth"
  bap_require_file "$BAP_GROUND_TRUTH_ARCHIVE"
  printf '%s  %s\n' \
    1b9f778bc3150f246c0f37b48215588c48369e39e6603c2b15975f0ec1b18d51 \
    "$BAP_GROUND_TRUTH_ARCHIVE" | sha256sum --check -
  env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" \
    "$BAP_REPO_ROOT/scripts/extract_rnp_ground_truth.py" \
    --manifest "$BAP_MANIFEST" --archive "$BAP_GROUND_TRUTH_ARCHIVE" \
    --output-dir "$BAP_GROUND_TRUTH_ROOT"
fi
ground_truth_count="$(find "$BAP_GROUND_TRUTH_ROOT/ground_truth" -mindepth 1 -maxdepth 1 -type d | wc -l)"
[[ "$ground_truth_count" -eq 50 ]] || \
  bap_die "Expected 50 selected ground-truth systems, found $ground_truth_count"

msa_count=0
[[ -d "$BAP_MSA_ROOT" ]] && \
  msa_count="$(find "$BAP_MSA_ROOT" -mindepth 2 -maxdepth 2 -type f -name '*.csv' | wc -l)"
if [[ "$msa_count" -ne 50 ]]; then
  ((BAP_EXECUTE)) || bap_die "Expected 50 MSA CSV files, found $msa_count"
  ((allow_msa)) || bap_die \
    "MSAs are missing; rerun with --execute --allow-msa-server to submit sequences"
  bap_require_executable "$BAP_BOLTZ_ENV/bin/python"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_BOLTZ_ENV/bin/python" "$BAP_REPO_ROOT/scripts/generate_boltz_msas.py" \
    --manifest "$BAP_MANIFEST" --output-dir "$BAP_MSA_ROOT" \
    --server-url https://api.colabfold.com --pairing-strategy greedy --retries 3
fi
[[ "$msa_count" -eq 50 ]] || \
  msa_count="$(find "$BAP_MSA_ROOT" -mindepth 2 -maxdepth 2 -type f -name '*.csv' | wc -l)"
[[ "$msa_count" -eq 50 ]] || bap_die "Expected 50 MSA CSV files after preparation"

if ((BAP_EXECUTE)); then
  gate_ids=(
    7ftm__1__1.A__1.C__1.C
    8aqf__1__1.A__1.B__1.B
    8jmp__1__1.A__1.D__1.D
    8haq__1__1.A__1.C__1.C
    8iqt__1__1.A__1.B__1.B
  )
  gate_args=()
  for sample_id in "${gate_ids[@]}"; do gate_args+=(--sample-id "$sample_id"); done

  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/prepare_rnp_paired_inputs.py" \
    --manifest "$BAP_MANIFEST" --model nesso1 --output-dir "$BAP_EXPERIMENT/inputs/nesso1/all"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/prepare_rnp_paired_inputs.py" \
    --manifest "$BAP_MANIFEST" --model nesso1 --output-dir "$BAP_EXPERIMENT/inputs/nesso1/gate" \
    "${gate_args[@]}"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_BOLTZ_ENV/bin/python" "$BAP_REPO_ROOT/scripts/prepare_rnp_paired_inputs.py" \
    --manifest "$BAP_MANIFEST" --model boltz2 --output-dir "$BAP_EXPERIMENT/inputs/boltz2/all" \
    --msa-root "$BAP_MSA_ROOT"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_BOLTZ_ENV/bin/python" "$BAP_REPO_ROOT/scripts/prepare_rnp_paired_inputs.py" \
    --manifest "$BAP_MANIFEST" --model boltz2 --output-dir "$BAP_EXPERIMENT/inputs/boltz2/gate" \
    --msa-root "$BAP_MSA_ROOT" "${gate_args[@]}"
fi

for model in nesso1 boltz2; do
  bap_require_dir "$BAP_EXPERIMENT/inputs/$model/all"
  bap_require_dir "$BAP_EXPERIMENT/inputs/$model/gate"
  [[ "$(bap_count_files "$BAP_EXPERIMENT/inputs/$model/all" '*.yaml')" -eq 50 ]] || \
    bap_die "Expected 50 $model full inputs"
  [[ "$(bap_count_files "$BAP_EXPERIMENT/inputs/$model/gate" '*.yaml')" -eq 5 ]] || \
    bap_die "Expected five $model gate inputs"
done

env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" \
  "$BAP_REPO_ROOT/scripts/validate_boltz_msas.py" \
  --manifest "$BAP_MANIFEST" --msa-root "$BAP_MSA_ROOT" \
  --output "${TMPDIR:-/tmp}/exp015_msa_inventory.json"

bap_log "Validated 50 references, 50 MSAs, and 50 full/5 gate inputs per model."
