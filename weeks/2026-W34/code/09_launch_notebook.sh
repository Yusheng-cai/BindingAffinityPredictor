#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/common.sh"

usage() {
  cat <<'EOF'
Usage: 09_launch_notebook.sh [--check|--execute]

--check    Validate the notebook and the saved artifacts it reads.
--execute  Start the local Jupyter Notebook server with the canonical
           FEP+4 analysis notebook selected.

The notebook analyzes saved predictions. It never launches model inference.
EOF
}

if ! bap_parse_execute_only "${1:---check}"; then
  usage
  exit 0
fi

notebook="$BAP_REPO_ROOT/weeks/2026-W34/notebooks/fepplus4_nesso1_analysis.ipynb"
notebook_python="${BAP_NOTEBOOK_PYTHON:-python3}"
bap_require_file "$notebook"

if ((BAP_EXECUTE)); then
  cd "$BAP_REPO_ROOT"
  exec env -u PYTHONPATH "$notebook_python" -m notebook "$notebook"
fi

env -u PYTHONPATH "$notebook_python" - "$notebook" "$BAP_REPO_ROOT" <<'PY'
import importlib.util
import json
import sys
from pathlib import Path

notebook_path = Path(sys.argv[1])
root = Path(sys.argv[2])
document = json.loads(notebook_path.read_text(encoding="utf-8"))
assert document.get("nbformat") == 4
assert len(document.get("cells", [])) == 23

for package in ("notebook", "matplotlib", "numpy", "pandas"):
    if importlib.util.find_spec(package) is None:
        raise ModuleNotFoundError(f"Notebook dependency is missing: {package}")

required = (
    root / "data/manifests/fepplus4_87.json",
    root / "runs/exp007_fepplus4_boltz2_nesso1/nesso1/seed42/predictions.csv",
    root / "reports/exp007_fepplus4_boltz2_nesso1/nesso1_metrics.json",
)
missing = [str(path) for path in required if not path.is_file()]
if missing:
    raise FileNotFoundError(f"Notebook inputs are missing: {missing}")

print(f"Notebook check passed: {notebook_path} ({len(document['cells'])} cells)")
PY
