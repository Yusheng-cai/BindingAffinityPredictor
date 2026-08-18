#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 03_prepare_inputs.sh [--check|--execute] [--allow-msa-server]

--check             Validate the existing manifest, MSAs, and model inputs.
--execute           Recreate missing benchmark and model-native inputs.
--allow-msa-server  Explicitly permit submission of protein sequences to the
                    public ColabFold MMseqs2 server. Required with --execute
                    when any Boltz MSA is missing.
EOF
}

BAP_EXECUTE=0
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

raw_root="$BAP_REPO_ROOT/data/raw/fepplus4"
msa_root="$BAP_EXP007_RUN/boltz2/msas"

if [[ ! -f "$BAP_MANIFEST" ]]; then
  ((BAP_EXECUTE)) || bap_die "Manifest is missing: $BAP_MANIFEST"
  git -C "$BAP_OPENFF_SOURCE" cat-file -e "$BAP_OPENFF_REVISION^{commit}" || \
    bap_die "OpenFF source does not contain $BAP_OPENFF_REVISION"
  mkdir -p "$raw_root/rcsb" "$raw_root/protein_pdb"
  declare -A dated_targets=(
    [cdk2]=2019-12-13_cdk2
    [tyk2]=2020-02-07_tyk2
    [jnk1]=2019-09-23_jnk1
    [p38]=2019-12-09_p38
  )
  declare -A pdb_ids=([cdk2]=1H1Q [tyk2]=4GIH [jnk1]=2GMX [p38]=3FLY)
  for target in $BAP_TARGETS; do
    source_dir="data/${dated_targets[$target]}"
    git -C "$BAP_OPENFF_SOURCE" show \
      "$BAP_OPENFF_REVISION:$source_dir/00_data/ligands.yml" \
      > "$raw_root/${target}_ligands_da7c337.yml"
    git -C "$BAP_OPENFF_SOURCE" show \
      "$BAP_OPENFF_REVISION:$source_dir/00_data/target.yml" \
      > "$raw_root/${target}_target_da7c337.yml"
    git -C "$BAP_OPENFF_SOURCE" show \
      "$BAP_OPENFF_REVISION:$source_dir/01_protein/crd/protein.pdb" \
      > "$raw_root/protein_pdb/${target}.pdb"
    curl --fail --location \
      "https://www.rcsb.org/fasta/entry/${pdb_ids[$target]}/display" \
      --output "$raw_root/rcsb/${pdb_ids[$target]}.fasta"
  done
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/build_fepplus4_manifest.py"
fi

bap_require_file "$BAP_MANIFEST"
printf '%s  %s\n' 5f0fb23bb97ea84f3da153e050d5913b97da9bba091fb84cd556b9894fdac1e1 \
  "$BAP_MANIFEST" | sha256sum --check -

if [[ ! -d "$BAP_EXP007/inputs/nesso1" ]]; then
  bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_NESSO_ENV/bin/python" "$BAP_REPO_ROOT/scripts/prepare_fepplus4_inputs.py" \
    --manifest "$BAP_MANIFEST" --model nesso1 \
    --output-dir "$BAP_EXP007/inputs/nesso1"
fi

missing_msa=0
for target in $BAP_TARGETS; do
  [[ -d "$msa_root/$target" ]] || missing_msa=1
done
if ((missing_msa)); then
  ((BAP_EXECUTE)) || bap_die "One or more Boltz MSA directories are missing beneath $msa_root"
  ((allow_msa)) || bap_die "MSAs are missing; rerun with --execute --allow-msa-server"
  env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_BOLTZ_ENV/bin/python" "$BAP_REPO_ROOT/scripts/generate_boltz_msas.py" \
    --manifest "$BAP_MANIFEST" --output-dir "$msa_root" \
    --server-url https://api.colabfold.com --pairing-strategy greedy
fi

if [[ ! -d "$BAP_EXP007/inputs/boltz2" ]]; then
  bap_run env -u PYTHONPATH PYTHONPATH="$BAP_REPO_ROOT/src" \
    "$BAP_BOLTZ_ENV/bin/python" "$BAP_REPO_ROOT/scripts/prepare_fepplus4_inputs.py" \
    --manifest "$BAP_MANIFEST" --model boltz2 \
    --output-dir "$BAP_EXP007/inputs/boltz2" --msa-root "$msa_root"
fi

for target in $BAP_TARGETS; do
  bap_require_dir "$BAP_EXP007/inputs/nesso1/$target"
  bap_require_dir "$BAP_EXP007/inputs/boltz2/$target"
done

nesso_count="$(find "$BAP_EXP007/inputs/nesso1" -name '*.yaml' -type f | wc -l)"
boltz_count="$(find "$BAP_EXP007/inputs/boltz2" -name '*.yaml' -type f | wc -l)"
msa_count="$(find "$msa_root" -maxdepth 2 -name '*.csv' -type f | wc -l)"
[[ "$nesso_count" -eq 87 ]] || bap_die "Expected 87 Nesso inputs, found $nesso_count"
[[ "$boltz_count" -eq 87 ]] || bap_die "Expected 87 Boltz inputs, found $boltz_count"
[[ "$msa_count" -eq 6 ]] || bap_die "Expected six target-chain MSA CSV files, found $msa_count"

bap_log "The manifest, 87 inputs per model, and six target-chain MSAs are ready."
