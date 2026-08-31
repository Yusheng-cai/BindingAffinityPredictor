#!/usr/bin/env bash
set -euo pipefail

cat >&2 <<'EOF'
This is the weekly bootstrap template. Before publishing this week, replace it
with a read-only --check mode and an explicit --execute mode that installs all
pinned software and acquires every external artifact required by the workflow.

Requirements: docs/decisions/0004-weekly-fresh-machine-bootstrap.md
EOF
exit 2
