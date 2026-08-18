#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 02_download_assets.sh [--check|--execute]

--check    Check that all model and released benchmark assets exist.
--execute  Download missing pinned assets, verify checksums, and extract archives.

This stage transfers approximately 11 GB. Large files remain in cache/ and are
ignored by Git; only their identifiers and checksums are committed.
EOF
}

if ! bap_parse_execute_only "${1:---check}"; then
  usage
  exit 0
fi

download_url() {
  local url="$1"
  local destination="$2"
  local expected_sha="$3"
  if [[ ! -f "$destination" ]]; then
    bap_run mkdir -p "$(dirname "$destination")"
    bap_run curl --fail --location --continue-at - --output "$destination" "$url"
  fi
  if ((BAP_EXECUTE)); then
    printf '%s  %s\n' "$expected_sha" "$destination" | sha256sum --check -
  else
    bap_require_file "$destination"
  fi
}

download_url https://model-gateway.boltz.bio/boltz2_conf.ckpt \
  "$BAP_BOLTZ_CACHE/boltz2_conf.ckpt" \
  090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1
download_url https://model-gateway.boltz.bio/boltz2_aff.ckpt \
  "$BAP_BOLTZ_CACHE/boltz2_aff.ckpt" \
  dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e
download_url https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar \
  "$BAP_BOLTZ_CACHE/mols.tar" \
  39e076d96dbec6b4e86982bbda16f3a53a2a60c9bdc17828d88f6f9a0c7d1fd7

if [[ ! -d "$BAP_BOLTZ_CACHE/mols" ]]; then
  bap_run tar -xf "$BAP_BOLTZ_CACHE/mols.tar" -C "$BAP_BOLTZ_CACHE"
fi

if ((BAP_EXECUTE)); then
  mkdir -p "$BAP_BOLTZ_ENV/checkpoints"
  [[ -e "$BAP_BOLTZ_STRUCTURE_CHECKPOINT" ]] || \
    ln -s "$BAP_BOLTZ_CACHE/boltz2_conf.ckpt" "$BAP_BOLTZ_STRUCTURE_CHECKPOINT"
  [[ -e "$BAP_BOLTZ_AFFINITY_CHECKPOINT" ]] || \
    ln -s "$BAP_BOLTZ_CACHE/boltz2_aff.ckpt" "$BAP_BOLTZ_AFFINITY_CHECKPOINT"
fi
bap_require_file "$BAP_BOLTZ_STRUCTURE_CHECKPOINT"
bap_require_file "$BAP_BOLTZ_AFFINITY_CHECKPOINT"

bap_require_executable "$BAP_NESSO_ENV/bin/hf"
if ((BAP_EXECUTE)); then
  HF_HOME="$BAP_NESSO_CACHE/huggingface" "$BAP_NESSO_ENV/bin/hf" download \
    recursionpharma/nesso ccd.pkl v1.0.0/hparams.json v1.0.0/model.safetensors \
    --revision "$BAP_NESSO_MODEL_REVISION" \
    --cache-dir "$BAP_NESSO_CACHE/huggingface"
  HF_HOME="$BAP_NESSO_CACHE/huggingface" "$BAP_NESSO_ENV/bin/hf" download \
    facebook/esm2_t33_650M_UR50D \
    --revision "$BAP_NESSO_ESM_REVISION" \
    --cache-dir "$BAP_NESSO_CACHE/huggingface"
else
  bap_require_dir "$BAP_NESSO_CACHE/huggingface/models--recursionpharma--nesso"
  bap_require_dir "$BAP_NESSO_CACHE/huggingface/models--facebook--esm2_t33_650M_UR50D"
fi

if ((BAP_EXECUTE)); then
  mkdir -p "$BAP_FLASHBIND_CACHE/checkpoints" "$BAP_FLASHBIND_CACHE/datasets"
  "$BAP_NESSO_ENV/bin/hf" download clorf6/FlashBind \
    value_1.ckpt value_2.ckpt \
    --revision "$BAP_FLASHBIND_MODEL_REVISION" \
    --local-dir "$BAP_FLASHBIND_CACHE/checkpoints"
  "$BAP_NESSO_ENV/bin/hf" download clorf6/FlashBind fep4.tar.zst \
    --repo-type dataset \
    --revision "$BAP_FLASHBIND_DATA_REVISION" \
    --local-dir "$BAP_FLASHBIND_CACHE/datasets"
fi

bap_require_file "$BAP_FLASHBIND_CHECKPOINT_1"
bap_require_file "$BAP_FLASHBIND_CHECKPOINT_2"
bap_require_file "$BAP_FLASHBIND_CACHE/datasets/fep4.tar.zst"

if ((BAP_EXECUTE)); then
  printf '%s  %s\n' ca87f84dedd4642da693a508b1119c802caa370e356dad31c47b8a708821e553 \
    "$BAP_FLASHBIND_CHECKPOINT_1" | sha256sum --check -
  printf '%s  %s\n' 59fa2686b0859a632862aec21ea69f2bebb1ab104e48606fe739c929368d488d \
    "$BAP_FLASHBIND_CHECKPOINT_2" | sha256sum --check -
  printf '%s  %s\n' c09047cb8c2e8212816b6bdbe1581f56a696b13773c8ec743da4d12b0c14678d \
    "$BAP_FLASHBIND_CACHE/datasets/fep4.tar.zst" | sha256sum --check -
fi

if [[ ! -d "$BAP_FLASHBIND_DATA" ]]; then
  bap_run mkdir -p "$BAP_FLASHBIND_CACHE/datasets/fep4"
  bap_require_executable "$BAP_FLASHBIND_ENV/bin/unzstd"
  bap_run tar --use-compress-program="$BAP_FLASHBIND_ENV/bin/unzstd" -xf \
    "$BAP_FLASHBIND_CACHE/datasets/fep4.tar.zst" \
    -C "$BAP_FLASHBIND_CACHE/datasets/fep4"
fi
bap_require_file "$BAP_FLASHBIND_DATA/id.json"
bap_require_file "$BAP_FLASHBIND_DATA/ligand_sdf.lmdb/data.mdb"

bap_log "Pinned model and FlashBind released-data assets are ready."
