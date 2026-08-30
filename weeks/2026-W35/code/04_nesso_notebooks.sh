#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 04_nesso_notebooks.sh [--check|--execute]

--check    Validate the two saved-output notebooks and run the focused mini-Nesso tests.
--execute  Rebuild/execute the paired result notebook, then validate both notebooks.

The Nesso architecture notebook is intentionally not overwritten: its saved
interactive checkpoint trace requires the local ignored exp012 tensors. Its
durable builder and tested teaching operations are validated here.
EOF
}

case "${1:---check}" in
  --check) BAP_EXECUTE=0 ;;
  --execute) BAP_EXECUTE=1 ;;
  -h|--help) usage; exit 0 ;;
  *) bap_die "Expected --check or --execute" ;;
esac

paired_notebook="$BAP_WEEK_DIR/notebooks/boltz2_nesso1_rnp_postcutoff50.ipynb"
architecture_notebook="$BAP_WEEK_DIR/notebooks/nesso_architecture_visualized.ipynb"
bap_require_file "$BAP_REPO_ROOT/scripts/build_rnp_paired_comparison_notebook.py"
bap_require_file "$BAP_REPO_ROOT/scripts/build_nesso_architecture_notebook.py"
bap_require_file "$BAP_REPO_ROOT/src/affinity_benchmark/educational/mini_nesso.py"
bap_require_file "$BAP_REPO_ROOT/tests/test_mini_nesso.py"

python3 - \
  "$BAP_REPO_ROOT/scripts/build_rnp_paired_comparison_notebook.py" \
  "$BAP_REPO_ROOT/scripts/build_nesso_architecture_notebook.py" <<'PY'
from pathlib import Path
import sys

for name in sys.argv[1:]:
    source = Path(name).read_text(encoding="utf-8")
    compile(source, name, "exec")
    print(f"{name}: syntax OK")
PY

if ((BAP_EXECUTE)); then
  python3 "$BAP_REPO_ROOT/scripts/build_rnp_paired_comparison_notebook.py"
fi

PYTHONPATH="$BAP_REPO_ROOT/src" python3 -m pytest -q "$BAP_REPO_ROOT/tests/test_mini_nesso.py"
python3 - "$paired_notebook" "$architecture_notebook" <<'PY'
import json
import sys

for name in sys.argv[1:]:
    with open(name, encoding="utf-8") as handle:
        notebook = json.load(handle)
    outputs = sum(bool(cell.get("outputs")) for cell in notebook.get("cells", []))
    if outputs == 0:
        raise SystemExit(f"{name}: no saved outputs")
    print(f"{name}: {len(notebook['cells'])} cells, {outputs} cells with saved outputs")
PY

bap_log "Validated the paired result notebook and Nesso architecture notebook."
