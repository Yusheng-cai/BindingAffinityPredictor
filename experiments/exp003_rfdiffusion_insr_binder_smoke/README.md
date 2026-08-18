# exp003: RFdiffusion insulin-receptor binder smoke test

## Status

**Engineering smoke test completed on 2026-08-11.** One RFdiffusion backbone,
both native trajectories, and four ProteinMPNN sequences were generated and
audited. Independent complex prediction and experimental validation have not
been performed.

## Scientific question

Can the official insulin-receptor PPI example be reproduced locally as an
auditable one-backbone experiment, including RFdiffusion's native reverse
trajectory and a subsequent ProteinMPNN sequence-design step?

This is a pipeline and representation test. It is not evidence that the output
binds insulin receptor.

## Inputs and conditioning

RFdiffusion receives the bundled 150-residue insulin-receptor target structure.
The target is present during generation and is therefore not a withheld test
structure. The compact contig instruction is:

```text
[A1-150/0 70-100]
```

This retains target residues A1 through A150, introduces a non-covalent chain
break, and requests a new 70- to 100-residue chain. Target residues A59, A83,
and A91 guide the interface location. They are geometrical RFdiffusion hotspot
labels, not experimental energetic measurements.

## Frozen procedure

1. Generate exactly one backbone using the complex checkpoint and deterministic
   design number/seed 42.
2. Use the official example's zero additional translation and frame noise.
3. Preserve both native 50-frame trajectories:
   - `Xt-1`: the sampled reverse state after each update;
   - `pX0`: the network's clean-structure estimate at each update.
4. Confirm the final chains, residue counts, trajectory frame counts, checkpoint,
   hardware, runtime, and native confidence array.
5. Keep the target sequence fixed and sample four binder sequences with
   ProteinMPNN `v_48_020`, seed 42, temperature 0.1, and no backbone noise.
6. For the teaching animation only, display the sample with the lowest
   designed-chain mean negative log probability. This is a sequence-model score,
   not a binding-affinity estimate.

## Why preserve two trajectories?

At a reverse step, `pX0` is RFdiffusion's current estimate of the fully denoised
structure. `Xt-1` is the state passed into the next step. They answer different
questions and must not be silently mixed. Neither series is molecular dynamics
or physical time.

## Success and interpretation

Engineering success requires a GPU run with a parseable final complex, one
generated chain in the requested length range, two 50-model trajectory files,
and four ProteinMPNN sequences for only the generated chain.

Even perfect engineering success establishes neither folding nor binding. A
separate structure-prediction filter is a later computational test, and direct
binding requires an experiment such as BLI or SPR.

## Observed result

- RFdiffusion exited successfully on an NVIDIA GeForce RTX 3080. Model-reported
  generation time was 66.01 s; total command wall time including model loading
  and schedule setup was 79.32 s.
- The selected length draw was 90 residues. The final PDB contains target chain
  A with 150 residues and generated chain B with 90 glycine placeholders.
- Both `Xt-1` and `pX0` trajectory PDBs contain 50 coordinate blocks. The files
  are written final-to-initial for PyMOL; the tutorial reverses them for display
  from noise level 50 to final step 1.
- RFdiffusion's pLDDT-like mean output for the binder rose from 0.0936 at
  `t=50` to 0.9919 at `t=1`. This is an internal structure-confidence output,
  not a calibrated probability of folding or binding.
- ProteinMPNN produced four 90-residue sequences in 4.46 s. Designed-chain
  scores were 0.7786, 0.8071, 0.7828 and 0.8508. The frozen selection rule picks
  sample 1, the lowest mean negative log-probability under ProteinMPNN.

The engineering success criteria were met. No conclusion can yet be made about
folding, expression, solubility, specificity, affinity, or experimental binding.

The tutorial maps ProteinMPNN sample 1 identities onto the generated binder's
N–Cα–C–O coordinates and colors those atoms by element. This is not an all-atom
model: binder side-chain conformations have not been packed or predicted.

## Local artifacts

Raw outputs are intentionally under the git-ignored directory:

```text
runs/exp003_rfdiffusion_insr_binder_smoke/
├── rfdiffusion/seed42/raw/
│   ├── design_ppi_42.pdb
│   ├── design_ppi_42.trb
│   └── traj/
└── proteinmpnn/seed42/raw/
    ├── probs/
    ├── scores/
    └── seqs/design_ppi_42.fa
```

The compact browser artifact is regenerated from these files with:

```bash
env -u PYTHONPATH external/envs/rfdiffusion/bin/python \
  scripts/export_rfdiffusion_tutorial_data.py \
  --run-root runs/exp003_rfdiffusion_insr_binder_smoke \
  --output docs/rfdiffusion_binder_tutorial/real-run-data.js
```

Artifact hashes and measured resource use are recorded in `experiment.yaml`.
