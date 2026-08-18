#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 01_setup_models.sh [--check|--execute]

--check    Verify the already installed source revisions and CLI entry points.
--execute  Clone missing source trees and create missing isolated environments.

The local CUDA-12.1 environments intentionally use --no_kernels during
inference. This script does not download model checkpoints.
EOF
}

if ! bap_parse_execute_only "${1:---check}"; then
  usage
  exit 0
fi

clone_pinned() {
  local url="$1"
  local revision="$2"
  local destination="$3"
  if [[ -d "$destination/.git" ]]; then
    bap_check_git_revision "$destination" "$revision"
    bap_log "Pinned source already present: $destination"
    return
  fi
  bap_run git clone "$url" "$destination"
  bap_run git -C "$destination" checkout --detach "$revision"
}

clone_pinned https://github.com/jwohlwend/boltz.git \
  "$BAP_BOLTZ_REVISION" "$BAP_BOLTZ_SOURCE"
clone_pinned https://github.com/recursionpharma/nesso.git \
  "$BAP_NESSO_REVISION" "$BAP_NESSO_SOURCE"
clone_pinned https://github.com/AIDD-Lab/FlashBind.git \
  "$BAP_FLASHBIND_REVISION" "$BAP_FLASHBIND_SOURCE"

if [[ ! -d "$BAP_OPENFF_SOURCE/.git" ]]; then
  bap_run git clone https://github.com/openforcefield/protein-ligand-benchmark.git \
    "$BAP_OPENFF_SOURCE"
fi
if ((BAP_EXECUTE)); then
  if ! git -C "$BAP_OPENFF_SOURCE" cat-file -e "$BAP_OPENFF_REVISION^{commit}"; then
    git -C "$BAP_OPENFF_SOURCE" fetch origin "$BAP_OPENFF_REVISION"
  fi
elif ! git -C "$BAP_OPENFF_SOURCE" cat-file -e "$BAP_OPENFF_REVISION^{commit}"; then
  bap_die "$BAP_OPENFF_SOURCE does not contain benchmark revision $BAP_OPENFF_REVISION"
fi

if [[ ! -x "$BAP_BOLTZ_ENV/bin/python" ]]; then
  bap_require_executable "$BAP_CONDA_EXE"
  bap_run "$BAP_CONDA_EXE" create -y -p "$BAP_BOLTZ_ENV" \
    --file "$BAP_REPO_ROOT/environments/boltz-2.2.1/conda-explicit-linux-64.txt"
  if ((BAP_EXECUTE)); then
    env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/python" -m pip install \
      --extra-index-url https://download.pytorch.org/whl/cu121 \
      -r <(grep -v '^boltz==' "$BAP_REPO_ROOT/environments/boltz-2.2.1/resolved-packages.txt")
    env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/python" -m pip install \
      --no-deps -e "$BAP_BOLTZ_SOURCE"
  fi
fi

if [[ ! -x "$BAP_NESSO_ENV/bin/python" ]]; then
  bap_require_executable "$BAP_CONDA_EXE"
  bap_require_executable "$BAP_BOLTZ_ENV/bin/python"
  bap_run "$BAP_CONDA_EXE" create -y -p "$BAP_NESSO_ENV" --clone "$BAP_BOLTZ_ENV"
  if ((BAP_EXECUTE)); then
    env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" -m pip uninstall -y boltz
    env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" -m pip install \
      numpy==2.1.3 lightning==2.6.5 transformers==5.15.0 \
      safetensors==0.8.0 huggingface_hub==1.27.0
    env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" -m pip install \
      --no-deps -e "$BAP_NESSO_SOURCE"
  fi
fi

if [[ ! -x "$BAP_FLASHBIND_ENV/bin/python" ]]; then
  bap_require_executable "$BAP_CONDA_EXE"
  bap_require_executable "$BAP_NESSO_ENV/bin/python"
  bap_run "$BAP_CONDA_EXE" create -y -p "$BAP_FLASHBIND_ENV" --clone "$BAP_NESSO_ENV"
  if ((BAP_EXECUTE)); then
    env -u PYTHONPATH "$BAP_FLASHBIND_ENV/bin/python" -m pip install \
      lmdb==1.6.2 timeout-decorator==0.5.0 torch-geometric==2.6.1
  fi
fi

if ((BAP_EXECUTE == 0)); then
  bap_check_git_revision "$BAP_BOLTZ_SOURCE" "$BAP_BOLTZ_REVISION"
  bap_check_git_revision "$BAP_NESSO_SOURCE" "$BAP_NESSO_REVISION"
  bap_check_git_revision "$BAP_FLASHBIND_SOURCE" "$BAP_FLASHBIND_REVISION"
  git -C "$BAP_OPENFF_SOURCE" cat-file -e "$BAP_OPENFF_REVISION^{commit}"
  bap_require_executable "$BAP_BOLTZ_ENV/bin/boltz"
  bap_require_executable "$BAP_NESSO_ENV/bin/nesso"
  bap_require_executable "$BAP_FLASHBIND_ENV/bin/python"
  env -u PYTHONPATH "$BAP_BOLTZ_ENV/bin/python" -c \
    'import boltz, torch; print(f"Boltz environment: torch={torch.__version__}, cuda={torch.version.cuda}")'
  env -u PYTHONPATH "$BAP_NESSO_ENV/bin/python" -c \
    'import nesso, torch; print(f"Nesso environment: torch={torch.__version__}, cuda={torch.version.cuda}")'
  env -u PYTHONPATH "$BAP_FLASHBIND_ENV/bin/python" -c \
    'import importlib.metadata as m, torch; print("FlashBind environment:", "torch=" + torch.__version__, "pyg=" + m.version("torch-geometric"), "lmdb=" + m.version("lmdb"))'
fi

bap_log "Model source and environment stage is ready."
