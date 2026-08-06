# Project Agent Instructions

These instructions apply to the entire repository.

## Before working

1. Read `README.md` and `docs/agentic-workflow.md`.
2. Read the relevant experiment definition and model/benchmark configuration.
3. Check `git status` and preserve unrelated user changes.
4. Ask before downloading datasets or checkpoints, installing software, changing an experiment's scientific protocol, or performing broad/destructive edits.

## Sources of truth

- Experiment intent and protocol: `experiments/<experiment-id>/`.
- Reusable model settings and pinned provenance: `configs/models/`.
- Benchmark definitions and metric policy: `configs/benchmarks/`.
- Input identity and experimental labels: `data/manifests/`.
- Reusable implementation: `src/affinity_benchmark/`.
- Reviewed conclusions: `reports/`.

Do not make a notebook or an ad hoc shell command the only record of an experiment.

## Reproducibility requirements

Every substantive run must record:

- experiment ID and Git revision;
- model source URL and revision;
- checkpoint identifier and checksum when available;
- exact command and resolved configuration;
- input-manifest version or checksum;
- random seed and sampling/recycling settings;
- hardware, software environment, runtime, peak memory, and failure status;
- raw model outputs before score normalization.

Keep binder probability and continuous affinity outputs separate. Record score units and direction explicitly; do not silently mix `pX`, `log10(X/M)`, and `log10(X/uM)`.

## Data and generated artifacts

- Never commit checkpoints, model caches, MSAs, downloaded raw datasets, processed bulk datasets, or raw prediction directories.
- Commit only small manifests, configurations, tests, reviewed reports, and compact figures/tables.
- Preserve experimental qualifiers such as `<` and `>` rather than silently treating them as exact measurements.
- Use the exact experimental protein construct and stereochemical ligand representation when available.

## Implementation

- Give each model a thin adapter under `src/affinity_benchmark/adapters/`.
- Keep model-specific dependencies in separate external environments.
- Put shared normalization, schemas, metrics, and provenance capture in reusable code rather than individual experiment scripts.
- Add a focused test when changing parsing, score conversion, manifest validation, or metrics.
- Record durable scientific or architectural choices under `docs/decisions/`.

Downloaded or cloned source belongs under `/home/yusheng/source`; locally installed programs belong under `/home/yusheng/programs`, unless the user approves another location.

