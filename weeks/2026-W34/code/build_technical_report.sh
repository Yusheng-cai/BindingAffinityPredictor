#!/usr/bin/env bash
set -euo pipefail

code_dir="$(cd "$(dirname "$0")" && pwd)"
week_dir="$(cd "$code_dir/.." && pwd)"
report_dir="$week_dir/technical-report"
repo_root="$(cd "$week_dir/../.." && pwd)"

install -m 0644 "$repo_root/literature/references.bib" \
  "$report_dir/references.bib"
mkdir -p "$report_dir/build"
(
  cd "$report_dir"
  latexmk \
    -pdf \
    -interaction=nonstopmode \
    -halt-on-error \
    -outdir=build \
    main.tex
)
install -m 0644 "$report_dir/build/main.pdf" "$report_dir/report.pdf"

pdfinfo "$report_dir/report.pdf" | awk '/^Pages:/ {print "pages:", $2}'
sha256sum "$report_dir/report.pdf"
