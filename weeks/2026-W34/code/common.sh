#!/usr/bin/env bash
# Small shared Bash helpers for the 2026-W34 workflow.

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

BAP_EXECUTE=0

bap_log() {
  printf '[2026-W34] %s\n' "$*"
}

bap_die() {
  printf '[2026-W34] ERROR: %s\n' "$*" >&2
  exit 1
}

bap_quote_command() {
  printf '  '
  printf '%q ' "$@"
  printf '\n'
}

bap_run() {
  bap_quote_command "$@"
  if ((BAP_EXECUTE)); then
    "$@"
  fi
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

bap_parse_execute_only() {
  case "${1:---check}" in
    --check)
      BAP_EXECUTE=0
      ;;
    --execute)
      BAP_EXECUTE=1
      ;;
    -h|--help)
      return 2
      ;;
    *)
      bap_die "Expected --check or --execute, received: $1"
      ;;
  esac
}

bap_git_revision() {
  git -C "$1" rev-parse HEAD 2>/dev/null || true
}

bap_check_git_revision() {
  local path="$1"
  local expected="$2"
  local observed
  bap_require_dir "$path/.git"
  observed="$(bap_git_revision "$path")"
  [[ "$observed" == "$expected" ]] || bap_die \
    "$path is at $observed; expected pinned revision $expected"
}

bap_skip_complete_run() {
  local run_dir="$1"
  local run_record="$run_dir/run.json"
  if [[ -f "$run_record" ]] && grep -q '"status": "complete"' "$run_record"; then
    bap_log "Complete run already exists; preserving it: $run_record"
    return 0
  fi
  if [[ -e "$run_dir" ]]; then
    bap_die "Run directory already exists but is not recorded complete: $run_dir"
  fi
  return 1
}
