#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 00_setup.sh [--check|--execute] [options]

--check                         Verify prerequisites, pinned software, model
                                assets, and the minimal Runs N' Poses files.
--execute                       Create or download missing items.
--allow-downloads               Explicitly permit GitHub, PyPI/Conda,
                                Hugging Face, Boltz, and Zenodo downloads.
--accept-af3-output-terms       Confirm that you reviewed and accept the
                                AlphaFold 3 Output Terms governing the Runs
                                N' Poses experimental structure archive.
--component NAME                Limit work to one of: all, data-tools, models,
                                model-assets, or runs-n-poses (default: all).

The Runs N' Poses terms are available before download at:
https://zenodo.org/records/18366081

Existing source checkouts at a different Git revision are never reset or
overwritten. Source is placed under BAP_SOURCE_ROOT, environments under
BAP_PROGRAM_ROOT, caches under BAP_CACHE_ROOT, and downloaded benchmark files
under this repository's ignored data/raw directory.
EOF
}

allow_downloads=0
accept_af3_terms=0
component=all

while (($#)); do
  case "$1" in
    --check) BAP_EXECUTE=0 ;;
    --execute) BAP_EXECUTE=1 ;;
    --allow-downloads) allow_downloads=1 ;;
    --accept-af3-output-terms) accept_af3_terms=1 ;;
    --component)
      shift
      (($#)) || bap_die "--component requires a value"
      component="$1"
      ;;
    --component=*) component="${1#*=}" ;;
    -h|--help) usage; exit 0 ;;
    *) bap_die "Unknown argument: $1" ;;
  esac
  shift
done

case "$component" in
  all|data-tools|models|model-assets|runs-n-poses) ;;
  *) bap_die "Unknown component: $component" ;;
esac

selected() {
  [[ "$component" == all || "$component" == "$1" ]]
}

require_download_permission() {
  ((BAP_EXECUTE)) || bap_die "Missing local item; rerun with --execute"
  ((allow_downloads)) || bap_die \
    "A network download is required; rerun with --execute --allow-downloads"
}

ensure_command() {
  command -v "$1" >/dev/null 2>&1 || bap_die "Required command is unavailable: $1"
}

verify_editable_install() {
  local python="$1"
  local package="$2"
  local expected_source="$3"
  local observed_source
  observed_source="$(env -u PYTHONPATH "$python" -m pip show "$package" | \
    sed -n 's/^Editable project location: //p')"
  [[ "$observed_source" == "$expected_source" ]] || bap_die \
    "$package is installed from '$observed_source'; expected pinned source '$expected_source'"
}

ensure_pinned_repo() {
  local url="$1"
  local revision="$2"
  local destination="$3"

  if [[ -d "$destination/.git" ]]; then
    bap_check_git_revision "$destination" "$revision"
    bap_log "Verified source revision: $destination"
    return
  fi
  [[ ! -e "$destination" ]] || bap_die \
    "Source destination exists but is not a Git checkout: $destination"
  require_download_permission
  mkdir -p "$(dirname "$destination")"
  git init "$destination"
  git -C "$destination" remote add origin "$url"
  git -C "$destination" fetch --depth 1 origin "$revision"
  git -C "$destination" checkout --detach FETCH_HEAD
  bap_check_git_revision "$destination" "$revision"
}

download_verified() {
  local url="$1"
  local destination="$2"
  local expected="$3"
  local partial="$destination.part"

  if [[ -f "$destination" ]]; then
    bap_verify_sha256 "$destination" "$expected"
    return
  fi
  require_download_permission
  mkdir -p "$(dirname "$destination")"
  curl --fail --location --retry 5 --retry-delay 2 \
    --continue-at - --output "$partial" "$url"
  bap_verify_sha256 "$partial" "$expected"
  mv "$partial" "$destination"
}

ensure_data_tools() {
  if [[ ! -x "$BAP_DATA_PYTHON" ]]; then
    require_download_permission
    bap_require_executable "$BAP_CONDA_EXE"
    mkdir -p "$BAP_PROGRAM_ROOT"
    "$BAP_CONDA_EXE" create -y -p "$BAP_DATA_ENV" python=3.9.23 pip
    env -u PYTHONPATH "$BAP_DATA_PYTHON" -m pip install \
      -r "$BAP_REPO_ROOT/environments/week35-data/requirements.txt"
  fi
  bap_require_executable "$BAP_DATA_PYTHON"
  env -u PYTHONPATH "$BAP_DATA_PYTHON" - <<'PY'
import sys
import numpy, pandas, pyarrow
if sys.version_info[:3] != (3, 9, 23):
    raise SystemExit(f"Week 35 data Python mismatch: {sys.version.split()[0]}; expected 3.9.23")
expected = {"numpy": "1.26.4", "pandas": "2.3.3", "pyarrow": "19.0.0"}
observed = {"numpy": numpy.__version__, "pandas": pandas.__version__, "pyarrow": pyarrow.__version__}
if observed != expected:
    raise SystemExit(f"Week 35 data environment mismatch: {observed}; expected {expected}")
print("Week 35 data environment:", observed)
PY
}

ensure_model_software() {
  ensure_command git
  ensure_pinned_repo https://github.com/jwohlwend/boltz.git \
    "$BAP_BOLTZ_REVISION" "$BAP_BOLTZ_SOURCE"
  ensure_pinned_repo https://github.com/recursionpharma/nesso.git \
    "$BAP_NESSO_REVISION" "$BAP_NESSO_SOURCE"

  if [[ ! -x "$BAP_BOLTZ_ENV/bin/python" ]]; then
    require_download_permission
    bap_require_executable "$BAP_CONDA_EXE"
    mkdir -p "$BAP_PROGRAM_ROOT"
    "$BAP_CONDA_EXE" create -y -p "$BAP_BOLTZ_ENV" \
      --file "$BAP_REPO_ROOT/environments/boltz-2.2.1/conda-explicit-linux-64.txt"
    env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/python" -m pip install \
      --extra-index-url https://download.pytorch.org/whl/cu121 \
      -r <(grep -v '^boltz==' "$BAP_REPO_ROOT/environments/boltz-2.2.1/resolved-packages.txt")
    env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/python" -m pip install \
      --no-deps -e "$BAP_BOLTZ_SOURCE"
  fi

  if [[ ! -x "$BAP_NESSO_ENV/bin/python" ]]; then
    require_download_permission
    bap_require_executable "$BAP_CONDA_EXE"
    bap_require_executable "$BAP_BOLTZ_ENV/bin/python"
    "$BAP_CONDA_EXE" create -y -p "$BAP_NESSO_ENV" --clone "$BAP_BOLTZ_ENV"
    env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" -m pip uninstall -y boltz
    env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" -m pip install \
      -r "$BAP_REPO_ROOT/environments/nesso-1.0.0/requirements-week35.txt"
    env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" -m pip install \
      --no-deps -e "$BAP_NESSO_SOURCE"
  fi

  bap_require_executable "$BAP_BOLTZ_ENV/bin/boltz"
  bap_require_executable "$BAP_NESSO_ENV/bin/nesso"
  verify_editable_install "$BAP_BOLTZ_ENV/bin/python" boltz "$BAP_BOLTZ_SOURCE"
  verify_editable_install "$BAP_NESSO_ENV/bin/python" nesso "$BAP_NESSO_SOURCE"
  env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/python" -m pip check
  env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" -m pip check
  env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/python" - <<'PY'
import boltz, torch
if not torch.cuda.is_available():
    raise SystemExit("Boltz environment imports, but PyTorch cannot see a CUDA GPU")
print(f"Boltz environment: torch={torch.__version__}, CUDA={torch.version.cuda}, GPU={torch.cuda.get_device_name(0)}")
PY
  env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" - <<'PY'
import nesso, torch
if not torch.cuda.is_available():
    raise SystemExit("Nesso environment imports, but PyTorch cannot see a CUDA GPU")
print(f"Nesso environment: torch={torch.__version__}, CUDA={torch.version.cuda}, GPU={torch.cuda.get_device_name(0)}")
PY
}

ensure_model_assets() {
  local boltz_conf="$BAP_BOLTZ_CACHE/boltz2_conf.ckpt"
  local boltz_mols="$BAP_BOLTZ_CACHE/mols.tar"
  local nesso_snapshot="$BAP_NESSO_CACHE/huggingface/models--recursionpharma--nesso/snapshots/$BAP_NESSO_MODEL_REVISION"
  local esm_snapshot="$BAP_NESSO_CACHE/huggingface/models--facebook--esm2_t33_650M_UR50D/snapshots/$BAP_NESSO_ESM_REVISION"

  ensure_command curl
  ensure_command tar
  download_verified https://model-gateway.boltz.bio/boltz2_conf.ckpt \
    "$boltz_conf" 090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1
  download_verified https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar \
    "$boltz_mols" 39e076d96dbec6b4e86982bbda16f3a53a2a60c9bdc17828d88f6f9a0c7d1fd7

  if [[ ! -d "$BAP_BOLTZ_CACHE/mols" ]]; then
    ((BAP_EXECUTE)) || bap_die "Extracted Boltz chemical components are missing"
    mkdir -p "$BAP_BOLTZ_CACHE"
    tar -xf "$boltz_mols" -C "$BAP_BOLTZ_CACHE"
  fi
  local mol_count
  mol_count="$(find "$BAP_BOLTZ_CACHE/mols" -type f -name '*.pkl' | wc -l)"
  [[ "$mol_count" -eq 45227 ]] || bap_die \
    "Expected 45,227 Boltz chemical-component files, found $mol_count"

  if [[ ! -f "$nesso_snapshot/v1.0.0/model.safetensors" || \
        ! -f "$nesso_snapshot/v1.0.0/hparams.json" || \
        ! -f "$nesso_snapshot/ccd.pkl" ]]; then
    require_download_permission
    bap_require_executable "$BAP_NESSO_ENV/bin/hf"
    HF_HOME="$BAP_NESSO_CACHE/huggingface" "$BAP_NESSO_ENV/bin/hf" download \
      recursionpharma/nesso ccd.pkl v1.0.0/hparams.json v1.0.0/model.safetensors \
      --revision "$BAP_NESSO_MODEL_REVISION" \
      --cache-dir "$BAP_NESSO_CACHE/huggingface"
  fi
  if [[ ! -f "$esm_snapshot/model.safetensors" ]]; then
    require_download_permission
    bap_require_executable "$BAP_NESSO_ENV/bin/hf"
    HF_HOME="$BAP_NESSO_CACHE/huggingface" "$BAP_NESSO_ENV/bin/hf" download \
      facebook/esm2_t33_650M_UR50D \
      --revision "$BAP_NESSO_ESM_REVISION" \
      --cache-dir "$BAP_NESSO_CACHE/huggingface"
  fi
  bap_verify_sha256 "$nesso_snapshot/v1.0.0/model.safetensors" \
    9928a8a824d147d665e76656804af1cd91c86731516d0b561f9fd1c91ee45622
  bap_verify_sha256 "$nesso_snapshot/v1.0.0/hparams.json" \
    06aa0c44fcd44eaa2c5c2473bfcf0d9a3af892d867e5279859bb22eb0f6e2d72
  bap_verify_sha256 "$nesso_snapshot/ccd.pkl" \
    7ed0ccd3903f19627926a5e41a5c0b5309c127a071cb9498dcae96926799ffd8
  bap_verify_sha256 "$esm_snapshot/model.safetensors" \
    a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0

  if [[ ! -f "$BAP_BOLTZ_CHECKPOINT" ]]; then
    ((BAP_EXECUTE)) || bap_die "Boltz staged checkpoint is missing: $BAP_BOLTZ_CHECKPOINT"
    mkdir -p "$(dirname "$BAP_BOLTZ_CHECKPOINT")"
    ln -s "$boltz_conf" "$BAP_BOLTZ_CHECKPOINT"
  fi
  bap_verify_sha256 "$BAP_BOLTZ_CHECKPOINT" \
    090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1
  bap_log "Verified pinned Boltz-2 and Nesso-1 model assets."
}

ensure_runs_n_poses() {
  local root="$BAP_REPO_ROOT/data/raw/runs_n_poses_metadata"
  local base="https://zenodo.org/api/records/18366081/files"

  ensure_command curl
  download_verified "$base/OUTPUT_TERMS_OF_USE.md/content" \
    "$root/OUTPUT_TERMS_OF_USE.md" \
    7ab1c4ef16a8fe88db3ccef813a5e28dec63b1e2d336876e5ae60b6d6b4686b8
  download_verified "$base/Legally_Binding_Terms_of_Use.txt/content" \
    "$root/Legally_Binding_Terms_of_Use.txt" \
    7de1599ca607a6c8df07a6b5829a9cd6fd8c864b156f95627600b556547b77f4

  if [[ ! -f "$root/annotations.csv" || ! -f "$root/inputs.json" || \
        ! -f "$root/all_similarity_scores.parquet" || \
        ! -f "$root/ground_truth.tar.gz" ]]; then
    ((accept_af3_terms)) || bap_die \
      "Runs N' Poses data are missing. Review the Zenodo terms, then rerun with --accept-af3-output-terms."
  fi
  download_verified "$base/annotations.csv/content" "$root/annotations.csv" \
    259aa0a8f5ea6008d4036886f1ad5ae255689c223eb839fa1a2e7cb54adae609
  download_verified "$base/inputs.json/content" "$root/inputs.json" \
    85ded25555c1efba0aeb7dc90df2fbf24d9550f390ed3e8ac9b32cf11a0ede4d
  download_verified "$base/all_similarity_scores.parquet/content" \
    "$root/all_similarity_scores.parquet" \
    ce771a2439a91210b7e4ebd09a729fb75d95b6c8265f4f38bae95913f4016804
  download_verified "$base/ground_truth.tar.gz/content" "$root/ground_truth.tar.gz" \
    1b9f778bc3150f246c0f37b48215588c48369e39e6603c2b15975f0ec1b18d51
  bap_log "Verified the pinned minimal Runs N' Poses record (Zenodo 18366081)."
}

ensure_command sha256sum
if selected data-tools; then ensure_data_tools; fi
if selected models; then ensure_model_software; fi
if selected model-assets; then ensure_model_assets; fi
if selected runs-n-poses; then ensure_runs_n_poses; fi

bap_log "Stage 00 is ready for component: $component"
