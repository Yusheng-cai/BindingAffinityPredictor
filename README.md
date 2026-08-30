# Binding Affinity Predictor

Reproducible experiments for understanding and comparing protein-ligand affinity models.

The project begins with Boltz-2 and is designed to add Nesso-1, FlashBind/FlashAffinity, classical scoring methods, and simple physicochemical baselines without changing the experiment interface.

## Initial scientific questions

1. Can Boltz-2 recover a known protein-ligand pose from sequence and SMILES?
2. How stable are its structures and affinity predictions across random seeds?
3. Can it rank ligands measured in one consistent assay?
4. Does it separate known binders from measured weak or inactive compounds?
5. How much information does it add beyond molecular weight and cLogP?

## Repository map

| Path | Purpose | Git policy |
| --- | --- | --- |
| `configs/` | Reusable model and benchmark settings | Tracked |
| `data/manifests/` | Small, versioned descriptions of inputs and labels | Tracked |
| `data/raw/`, `data/processed/` | Downloaded and derived datasets | Ignored |
| `docs/` | Scope, workflow, model notes, and decisions | Tracked |
| `environments/` | Environment specifications and lock files | Tracked |
| `experiments/` | Frozen questions, protocols, and experiment configurations | Tracked |
| `src/affinity_benchmark/` | Reusable adapters, data handling, and metrics | Tracked |
| `scripts/` | Thin command-line entry points and operational helpers | Tracked |
| `runs/` | Raw predictions, logs, and run metadata | Ignored except documentation |
| `reports/` | Compact reviewed summaries, tables, and figures | Tracked |
| `weeks/` | Self-contained weekly code, notebooks, compact results, and reports | Tracked |
| `notebooks/` | Cross-project and educational notebooks | Tracked |
| `tests/` | Unit and small integration tests | Tracked |

External model repositories and installed programs are not vendored here. Their exact source revisions, checkpoints, and environment details will be recorded in model configurations.

## Current status

The 2026-W34 workflow contains a completed released-model reproduction on the
87-compound FEP+ four-kinase benchmark for Nesso-1, Boltz-2 with an explicitly
documented 1,024-sequence MSA subsample, and FlashBind using the authors'
released FABind+ poses and representations. A 16-complex crystal subset also
compares Boltz-2 and FlashBind ligand-pose RMSD.

Start at `weeks/2026-W34/code/README.md` for the complete pinned setup,
download, input-preparation, inference, analysis, and report workflow. Large
external artifacts remain in ignored cache and run directories; Git tracks the
scientific protocols, source revisions, checksums, canonical code, compact
results, and weekly report.

For a presentation-first view, start at
[`weeks/2026-W34/README.md`](weeks/2026-W34/README.md).

The latest structural benchmark, `exp015`, compares Nesso-1 and Boltz-2 on a
frozen set of 50 diverse Runs N' Poses complexes released after June 2023.
Nesso-1 completed 50/50 systems and Boltz-2 completed 49/50 on the local 10 GiB
GPU. See the [paired result](reports/exp015_boltz2_nesso1_rnp_postcutoff50/README.md)
and its [executed notebook](weeks/2026-W35/notebooks/boltz2_nesso1_rnp_postcutoff50.ipynb).

## Local checks (no downloads)

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
PYTHONPATH=src python3 scripts/validate_manifest.py \
  data/manifests/exp001_tyk2_4gih.json
```

See `experiments/exp001_boltz2_4gih_smoke/README.md` for the scientific design,
`docs/model-notes/boltz2.md` for model interpretation, and
`docs/literature.md` for the annotated literature map.
