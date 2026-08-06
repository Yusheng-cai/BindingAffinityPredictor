# exp001: Boltz-2 TYK2–ligand 46 smoke test

## Status

**Protocol frozen; local checks pass; the model environment and checksummed
checkpoints are ready. No MSA, experimental coordinate, or prediction has been
generated.**

## Scientific question

Given only the exact deposited TYK2 construct sequence, the neutral ligand-46
SMILES, and an MSA—not the holo coordinates or binding-site residues—can the
pinned Boltz-2 setup complete inference and produce an interpretable complex,
confidence record, continuous affinity output, and binder probability on local
hardware?

Pose comparison to PDB 4GIH is a secondary readout. Affinity comparison to the
reported `Ki = 4.8 nM` is descriptive only.

## Why this system

- 4GIH is a 2.00 Å X-ray structure with one protein chain and a small, neutral,
  achiral ligand, minimizing avoidable input ambiguity.
- Ligand 46 belongs to an openly curated TYK2 congeneric benchmark and has a
  reported Ki with an error estimate.
- The complex is small enough to be a plausible first test on the available
  10 GiB RTX 3080, although memory feasibility remains to be measured.

Its weakness is equally important: it is a 2013 public structure and may be
familiar to Boltz-2. Therefore this experiment is a software/behavior check, not
a novelty-controlled benchmark.

## Withheld-information boundary

The prediction input contains sequence, SMILES, and an MSA. It contains no 4GIH
coordinates, structural template, pocket list, contact constraint, or ligand
starting pose. The crystal structure is acquired only for post-prediction
scoring. This prevents direct pose leakage through our inference input, but does
not remove possible training-data familiarity.

## Staged procedure

1. Validate the frozen manifest and logarithmic unit conversions locally.
2. Install the pinned Boltz source outside this repository and resolve/checksum
   the automatically selected model checkpoints. **Completed 2026-08-06.**
3. After separate approval for external sequence submission, generate the MSA
   with the configured ColabFold server. Cache its exact response and checksum.
4. Run seed 42 with the exact settings in `configs/models/boltz2.yaml`.
5. Confirm expected output files, finite affinity values, structure parseability,
   runtime, and peak GPU memory before any larger run.
6. Download 4GIH only after prediction and calculate preregistered structure and
   pose diagnostics.

No stage silently falls back to single-sequence mode; the official documentation
warns that it reduces accuracy.

The environment uses PyTorch CUDA 12.1 and the documented `--no_kernels`
fallback. This disables optional cuEquivariance acceleration kernels, not GPU
inference itself. The host's AmberTools `PYTHONPATH` must be unset for every
Boltz command; the exact invocation is recorded in the environment notes.

## Readouts

### Engineering readouts

- process completion and explicit failure reason;
- exact command, Git/source/checkpoint revisions, seed, environment, and hardware;
- wall time and peak GPU memory;
- presence and parseability of mmCIF, confidence JSON, and affinity JSON;
- finite native continuous-affinity and binder-probability fields.

### Structural readouts

- protein Cα RMSD over common experimentally resolved residues after optimal
  superposition;
- symmetry-corrected ligand heavy-atom RMSD after alignment on the protein
  pocket, with 2 Å or less recorded as the conventional pose-recovery heuristic;
- native confidence fields, ligand contact recovery, and PoseBusters checks.

A low ligand RMSD with chemically invalid geometry is not a successful pose.

### Affinity readouts

Record, without relabeling:

- experimental: `Ki = 0.0048 uM`, `log10(Ki/uM) = -2.31876`, and
  `pKi = 8.31876`;
- model continuous head: approximate `log10(IC50/uM)`, lower is stronger;
- model binary head: binder probability, higher is more binder-like.

No Pearson/Spearman correlation, RMSE, calibration statistic, or ranking claim
is valid at `n = 1`.

## Success criteria

The smoke test succeeds if the pinned run completes on local hardware, all
required provenance and native outputs are captured, and outputs parse without
schema or finite-value errors. Pose metrics are reported as scientific results,
not required for engineering success. Affinity agreement is not a pass/fail
criterion.

## Predeclared next step

If exp001 is technically successful, run the same protocol with two additional
seeds to estimate stochastic spread. Only then expand to all 16 consistently
annotated TYK2 ligands for within-series ranking. That follow-up must include a
ligand-only baseline and must retain the leakage warning.
