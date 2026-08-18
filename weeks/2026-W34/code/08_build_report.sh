#!/usr/bin/env bash
set -euo pipefail

code_dir="$(cd "$(dirname "$0")" && pwd)"
exec "$code_dir/build_report.sh" "$@"
