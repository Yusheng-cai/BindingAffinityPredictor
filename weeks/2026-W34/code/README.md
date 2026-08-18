# Complete 2026-W34 computational workflow

This directory is the single starting point for the calculations reported in
the 2026-W34 paper. It covers source acquisition, isolated Python environments,
model weights, benchmark data, input preparation, model inference, common
output tables, affinity statistics, crystal-pose RMSD, an interactive notebook,
figures, and PDF build.

The numbered files are Bash orchestration scripts. They call the tested Python
implementations under `scripts/` and `src/affinity_benchmark/`; those canonical
implementations are not duplicated here. Model repositories, environments,
checkpoints, MSAs, and raw predictions are intentionally excluded from Git.

## What is installed

Each model has four separate pieces:

1. a pinned source checkout under `/home/yusheng/source`;
2. an isolated Python environment under `/home/yusheng/programs`;
3. trained weights and model assets beneath the ignored project `cache/`;
4. a Python command or script that loads the source, dependencies, and weights.

Nesso and Boltz expose small Python command-line launchers named `nesso` and
`boltz`. FlashBind is invoked through its upstream `scripts/predict.py` file.
PyTorch supplies the compiled CUDA kernels used by all three.

## Safe execution model

Setup, download, input generation, analysis, and inference are separate. The
scripts do not download or run expensive calculations unless `--execute` is
given. Existing complete `run.json` records are preserved and skipped.

MSA generation additionally requires `--allow-msa-server`, because it submits
the four benchmark protein sequences to the public ColabFold MMseqs2 service.
Crystal-reference acquisition requires `--allow-rcsb-download`.

## Reproduce in order

From the repository root:

```bash
# Inspect the existing installation and cached assets.
bash weeks/2026-W34/code/01_setup_models.sh --check
bash weeks/2026-W34/code/02_download_assets.sh --check
bash weeks/2026-W34/code/03_prepare_inputs.sh --check

# On a new machine, explicitly perform the external setup steps.
bash weeks/2026-W34/code/01_setup_models.sh --execute
bash weeks/2026-W34/code/02_download_assets.sh --execute
bash weeks/2026-W34/code/03_prepare_inputs.sh \
  --execute --allow-msa-server

# Preview the exact three-model commands, then run them on a suitable GPU.
bash weeks/2026-W34/code/04_run_inference.sh --model=all --check
bash weeks/2026-W34/code/04_run_inference.sh --model=all --execute

# Parse native outputs and calculate the predefined statistics.
bash weeks/2026-W34/code/05_collect_outputs.sh --execute
bash weeks/2026-W34/code/06_analyze_affinity.sh --execute
bash weeks/2026-W34/code/07_analyze_poses.sh \
  --execute --allow-rcsb-download

# Regenerate compact figures and compile the weekly PDF.
bash weeks/2026-W34/code/08_build_report.sh

# Compile the separate technical methods report.
bash weeks/2026-W34/code/build_technical_report.sh

# Validate or launch the interactive saved-output walkthrough.
bash weeks/2026-W34/code/09_launch_notebook.sh --check
bash weeks/2026-W34/code/09_launch_notebook.sh --execute
```

A single target can be run first:

```bash
bash weeks/2026-W34/code/04_run_inference.sh \
  --model=nesso1 --target=tyk2 --execute
```

## Inputs and outputs

The shared affinity input is `data/manifests/fepplus4_87.json`: 87 neutral
ligands across CDK2, TYK2, JNK1, and p38. Boltz-2 and Nesso-1 receive matched
protein sequences and stereochemical SMILES. Boltz additionally receives the
reused target-chain MSAs. The FlashBind reproduction receives the authors'
released receptor structures, FABind+ ligand poses, and precomputed features.

Raw model files remain unchanged beneath `runs/`. `05_collect_outputs.sh`
creates one common CSV per model while preserving native fields. Binder
probability is never substituted for continuous affinity. Score units and
direction are defined in `configs/models/` and the experiment YAML files.

## Interactive notebook

The canonical notebook is
[`weeks/2026-W34/notebooks/fepplus4_nesso1_analysis.ipynb`](../notebooks/fepplus4_nesso1_analysis.ipynb).
It is an executable, presentation-friendly walkthrough of the 87 saved Nesso-1
predictions, target-level metrics, centered MAE, aggregate statistics, native
outputs, and recorded inference command. It reads canonical files and calls the
tested metric implementation; it does not rerun Nesso or hide a separate metric
implementation inside the notebook.

The current notebook focuses on Nesso-1 because it was created for the first
completed analysis. Boltz-2 and FlashBind remain fully covered by the numbered
pipeline and report; a future notebook revision can add their side-by-side
tables without changing the saved scientific results.

## Provenance and portability

`workflow_manifest.yaml` pins every source and model revision, checkpoint
checksum, benchmark revision, seed, and canonical Python entry point. Local
defaults are in `config.sh` and can be overridden without editing tracked code,
for example:

```bash
BAP_SOURCE_ROOT=/work/source \
BAP_PROGRAM_ROOT=/work/programs \
BAP_CACHE_ROOT=/scratch/binding-cache \
bash weeks/2026-W34/code/01_setup_models.sh --check
```

Environment reconstruction is Linux/CUDA-specific. The frozen scientific
protocol is portable, but another GPU driver may require a compatible PyTorch
build. Such a change must be recorded as an environment deviation; it must not
silently change seeds, MSA sampling, recycling, diffusion, or metric policy.

## GitHub and weekly updates

Commit this directory together with the canonical experiment definitions,
configs, manifests, tested Python code, compact results, report source, and
PDF. Do not commit `cache/`, `runs/`, external model repositories, environments,
or checkpoints. Every new `weeks/YYYY-Www/` folder should preserve the
week's report and workflow manifest at the project Git revision used to build
it. See `weeks/README.md` for the weekly publishing checklist.
