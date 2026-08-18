#!/usr/bin/env bash
set -euo pipefail

code_dir="$(cd "$(dirname "$0")" && pwd)"
week_dir="$(cd "$code_dir/.." && pwd)"

python3 "$code_dir/build_assets.py"
mkdir -p "$week_dir/report/build"
(
  cd "$week_dir/report"
  latexmk \
    -pdf \
    -interaction=nonstopmode \
    -halt-on-error \
    -outdir=build \
    main.tex
)
install -m 0644 \
  "$week_dir/report/build/main.pdf" \
  "$week_dir/report/report.pdf"

sha256sum "$week_dir/report/report.pdf"
