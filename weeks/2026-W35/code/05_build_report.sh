#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../../.." && pwd)"
report_dir="$repo_root/weeks/2026-W35/report"

python3 "$script_dir/build_assets.py"
install -m 0644 "$repo_root/literature/references.bib" "$report_dir/references.bib"
install -d "$report_dir/build"

(
  cd "$report_dir"
  latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
)

install -m 0644 "$report_dir/build/main.pdf" "$report_dir/report.pdf"
pdfinfo "$report_dir/report.pdf" | awk '/^Pages:/ {print "pages:", $2}'
sha256sum "$report_dir/report.pdf"
