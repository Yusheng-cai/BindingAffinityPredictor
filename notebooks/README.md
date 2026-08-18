# Notebooks

Notebooks are for exploration and visualization. They should read canonical predictions and call reusable functions from `src/affinity_benchmark/`; they should not be the sole implementation of parsing, normalization, or metrics.

## Nesso-1 FEP+4 benchmark walkthrough

The executable FEP+4 analysis-first walkthrough is stored with its weekly
report at
`weeks/2026-W34/notebooks/fepplus4_nesso1_analysis.ipynb`. It loads the frozen
87-compound manifest and completed Nesso-1 outputs, explains the meaning of
target, ligand count, IC50/Ki endpoint, construct, and reference PDB, and
derives Pearson, Spearman, Kendall, ordinary MAE, centered MAE, and the
compound-weighted aggregate. It also reproduces the report plots and displays
recorded run provenance.

The notebook never launches model inference automatically. Start it from the
repository root with:

```bash
env -u PYTHONPATH python3 -m notebook weeks/2026-W34/notebooks/fepplus4_nesso1_analysis.ipynb
```

## Educational architecture notebook

`mini_rfdiffusion_proteinmpnn.ipynb` is a synthetic, executable teaching
exercise. It illustrates C-alpha diffusion, equivariant coordinate updates,
residue-graph message passing, random-order sequence decoding, and transparent
interface-score reranking. The structural views use Plotly and support drag to
rotate, scroll to zoom, hover inspection, and animated trajectory sliders.
Reusable components live in
`src/affinity_benchmark/educational/mini_binder.py` and have focused tests.

The models are intentionally small approximations. They are not the production
RFdiffusion or ProteinMPNN architectures, are not trained on the PDB, and do not
predict experimental affinity or binding free energy.

Launch it from the repository root with the PyTorch-enabled Python interpreter:

```bash
python3 -m notebook notebooks/mini_rfdiffusion_proteinmpnn.ipynb
```

Using `python3 -m notebook` is intentional on this workstation: the standalone
`jupyter` executable resolves an older Python kernel that does not contain
PyTorch.

## Real RFdiffusion + ProteinMPNN walkthrough

`real_rfdiffusion_proteinmpnn_insr.ipynb` reads the actual, audited outputs of
`exp003_rfdiffusion_insr_binder_smoke`. It explains the target PDB, contig and
hotspot inputs; displays the real pX0 and Xt-1 RFdiffusion trajectories; audits
the final backbone and internal confidence; and inspects the four real
ProteinMPNN sequences, probability archive and scores. It also reads the
matched `exp005` hotspot-conditioning ablation, with rotatable guided versus
unguided seed pairs, an overlay of both ten-backbone ensembles, and the paired
minimum-distance result. The experimental-reference section aligns the 4OGA
L1–insulin–αCT site-1 complex to the target and compares each generated binder
with the experimentally occupied site. Apparent overlaps are explicitly treated
as steric/competitive geometry rather than affinity evidence.

The expensive commands are printed but are not launched automatically. The
notebook reuses the local raw artifacts under the Git-ignored `runs/` tree and
uses Plotly for rotatable structures and trajectory playback. Launch it from
the repository root with:

```bash
python3 -m notebook notebooks/real_rfdiffusion_proteinmpnn_insr.ipynb
```

This notebook is a genuine model-output walkthrough, but the experiment is
still only an engineering smoke test: no independent structure prediction or
experimental binding validation has been completed.
