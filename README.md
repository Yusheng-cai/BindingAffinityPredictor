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
| `notebooks/` | Exploratory analysis that calls reusable project code | Tracked |
| `tests/` | Unit and small integration tests | Tracked |

External model repositories and installed programs are not vendored here. Their exact source revisions, checkpoints, and environment details will be recorded in model configurations.

## Current status

The repository scaffold and the protocol for `exp001_boltz2_4gih_smoke` are in
place. The local, dependency-free tests validate the frozen one-complex manifest
and the affinity-unit conversions. The pinned Boltz source, isolated Python
environment, chemical-component data, and both Boltz-2 checkpoints are installed
and verified. No MSA, experimental coordinate file, or prediction has been
generated yet.

The first model run is intentionally one TYK2 complex and one seed. It is a
pipeline and interpretation smoke test, not an affinity benchmark. If that run
is technically successful, the next stages are repeated seeds and then the
complete single-assay TYK2 congeneric series.

## Local checks (no downloads)

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
PYTHONPATH=src python3 scripts/validate_manifest.py \
  data/manifests/exp001_tyk2_4gih.json
```

See `experiments/exp001_boltz2_4gih_smoke/README.md` for the scientific design,
`docs/model-notes/boltz2.md` for model interpretation, and
`docs/literature.md` for the annotated literature map.
