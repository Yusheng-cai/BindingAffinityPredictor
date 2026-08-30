#!/usr/bin/env bash
# Small shared Bash helpers for the 2026-W35 workflow.

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

BAP_EXECUTE=0

bap_log() {
  printf '[2026-W35] %s\n' "$*"
}

bap_die() {
  printf '[2026-W35] ERROR: %s\n' "$*" >&2
  exit 1
}

bap_require_file() {
  [[ -f "$1" ]] || bap_die "Required file is missing: $1"
}

bap_require_dir() {
  [[ -d "$1" ]] || bap_die "Required directory is missing: $1"
}

bap_require_executable() {
  [[ -x "$1" ]] || bap_die "Required executable is missing: $1"
}

bap_complete_run() {
  local run_record="$1/run.json"
  [[ -f "$run_record" ]] && grep -q '"status": "complete"' "$run_record"
}

bap_refuse_incomplete_run() {
  local run_dir="$1"
  [[ ! -e "$run_dir" ]] || bap_die \
    "Run directory exists without a complete run record; inspect it before retrying: $run_dir"
}

bap_count_files() {
  find "$1" -maxdepth 1 -type f -name "$2" | wc -l
}
