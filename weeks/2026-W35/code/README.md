# Week 2026-W35 reproduction workflow

This folder contains only the code needed for the three computational
components presented in Week 35:

1. the paired 50-system Runs N' Poses calculation for Nesso-1 and Boltz-2;
2. the timing measurements collected during those same model runs; and
3. the paired-results and Nesso-1 architecture notebooks.

The numbered Bash files are thin orchestration scripts. They call the tested,
canonical Python implementations under `scripts/` and
`src/affinity_benchmark/`; the scientific code is not duplicated inside the
weekly folder. Raw model outputs, MSAs, reference structures, checkpoints, and
environments remain in ignored directories and are not committed.

The familiarity-table step uses the existing Anaconda Python because it has a
Parquet reader (`pyarrow`); the diverse-selection step uses the isolated Nesso
environment because it has RDKit. `BAP_DATA_PYTHON`, `BAP_PROGRAM_ROOT`, and
`BAP_CACHE_ROOT` can be overridden for another machine without changing the
frozen scientific settings.

## Reproduce in order

Run these commands from the repository root. Every stage defaults to a
read-only check. Expensive work is performed only when `--execute` is supplied,
and completed run directories are preserved rather than overwritten.

```bash
# 1. Check source-data checksums, the frozen cohort, 50 references/MSAs, and inputs.
bash weeks/2026-W35/code/01_prepare_runs_n_poses.sh --check

# Recalculate June-2023 familiarity, reproduce the byte-identical 50-system
# selection, and recreate local inputs. The additional flag explicitly permits
# submission of protein sequences to the public ColabFold MSA server if needed.
bash weeks/2026-W35/code/01_prepare_runs_n_poses.sh \
  --execute --allow-msa-server

# 2. Check the recorded full Nesso-1 and Boltz-2 runs and print raw wall times.
bash weeks/2026-W35/code/02_run_models_and_time.sh --check

# On a new machine, run the five-system gate and then the full 50-system jobs.
bash weeks/2026-W35/code/02_run_models_and_time.sh --execute

# 3. Check or reproduce structure scoring and the paired/familiarity analysis.
bash weeks/2026-W35/code/03_analyze_runs_n_poses.sh --check
bash weeks/2026-W35/code/03_analyze_runs_n_poses.sh --execute

# 4. Validate the two notebooks and the miniature Nesso tensor operations.
bash weeks/2026-W35/code/04_nesso_notebooks.sh --check

# Re-execute the paired result notebook from compact reviewed results.
bash weeks/2026-W35/code/04_nesso_notebooks.sh --execute

# 5. Regenerate the two report figures and compile the PDF.
bash weeks/2026-W35/code/05_build_report.sh
```

## What the timing means

`02_run_models_and_time.sh` does not run a separate synthetic timing test. The
same 50-system inference commands used for the scientific benchmark are wrapped
by `scripts/run_with_provenance.py`. Each `run.json` records the exact command,
wall time, peak host/GPU memory, hardware, software, Git revision, settings
checksums, and failure status. The report divides each batched wall time by 50
attempted systems:

- Nesso-1: 326.07 s total, or 6.52 s per attempted system;
- Boltz-2 MSA-1024: 1181.73 s total, or 23.63 s per attempted system.

The Boltz-2 time excludes the earlier remote ColabFold MSA search, whereas
Nesso-1's local ESM-2 feature calculation is included. These are descriptive
single runs on one RTX 3080, not repeated latency distributions.

## Durable source files

The experiment definition and frozen data identity live outside the weekly
folder:

- `experiments/exp015_boltz2_nesso1_rnp_postcutoff50/` — protocol and exact
  gate/full model commands;
- `configs/benchmarks/rnp_boltz2_nesso1_postcutoff50.yaml` — cohort and metric
  policy;
- `data/manifests/rnp_boltz2_nesso1_postcutoff50.json` — 50 selected systems;
- `scripts/run_with_provenance.py` — timing and provenance capture;
- `scripts/extract_rnp_cutoff_similarity.py` and
  `scripts/select_rnp_boltz2_nesso1_postcutoff.py` — familiarity calculation
  and deterministic balanced/diverse cohort selection;
- `scripts/analyze_nesso_rnp_distograms.py` — Nesso distance/contact scoring;
- `scripts/analyze_boltz2_rnp_structures.py` — Boltz structure scoring;
- `scripts/compare_rnp_paired_models.py` — paired statistics and familiarity;
- `scripts/build_rnp_paired_comparison_notebook.py` — executed results notebook;
- `scripts/build_nesso_architecture_notebook.py` and
  `src/affinity_benchmark/educational/mini_nesso.py` — architecture companion;
- `reports/exp015_boltz2_nesso1_rnp_postcutoff50/` — reviewed compact result.

The saved Nesso architecture notebook is kept intact because its interactive
checkpoint trace reads ignored raw tensors from `exp012`. Rebuilding its source
cells would remove those saved outputs; the builder and the numerical teaching
operations are therefore tested without overwriting the presentation copy.
