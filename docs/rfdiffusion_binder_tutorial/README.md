# RFdiffusion binder tutorial

Dependency-free interactive explanation of the RFdiffusion → ProteinMPNN
protein-binder design pipeline, using the official insulin-receptor example.

The page now contains two complementary viewers:

- a real, rotatable 50-step Cα trajectory exported from experiment `exp003`,
  with both RFdiffusion trajectory conventions and the four actual ProteinMPNN
  candidates;
- a target-aligned overlay of ten real RFdiffusion backbones from experiment
  `exp004`, with per-design geometry and all 45 pairwise comparisons;
- a deliberately schematic animation that isolates the representations at each
  pipeline stage.

Both distinguish generative diffusion timesteps from molecular dynamics and
computational filtering from experimental evidence.

The companion page `proteinmpnn_guidance.html` explains why a chemistry-aware
interface score belongs primarily beside ProteinMPNN, compares reranking with
logit-guided and coupled sequence optimization, and defines a conservative
first validation experiment. It explicitly distinguishes an interface-quality
score from a binding free energy.

## Open locally

From this directory:

```bash
python3 -m http.server 8000
```

Then open <http://localhost:8000>.

The page also works when `index.html` is opened directly because it has no
external runtime dependencies.

## Regenerate the real-run asset

From the repository root, after running experiment `exp003`:

```bash
env -u PYTHONPATH external/envs/rfdiffusion/bin/python \
  scripts/export_rfdiffusion_tutorial_data.py \
  --run-root runs/exp003_rfdiffusion_insr_binder_smoke \
  --output docs/rfdiffusion_binder_tutorial/real-run-data.js
```

The compact JavaScript asset contains rounded Cα trajectories, final N–Cα–C–O
coordinates, ProteinMPNN residue identities and scores. It does not contain
invented side-chain coordinates. Raw PDB, TRB and NPZ outputs remain under the
git-ignored `runs/` tree.

Regenerate the ten-seed diversity asset and report tables with:

```bash
env -u PYTHONPATH external/envs/rfdiffusion/bin/python \
  scripts/analyze_rfdiffusion_diversity.py \
  --seed42-root runs/exp003_rfdiffusion_insr_binder_smoke/rfdiffusion/seed42/raw \
  --ensemble-root runs/exp004_rfdiffusion_insr_structural_diversity/rfdiffusion/raw \
  --report-dir reports/exp004_rfdiffusion_insr_structural_diversity \
  --web-output docs/rfdiffusion_binder_tutorial/diversity-data.js
```
