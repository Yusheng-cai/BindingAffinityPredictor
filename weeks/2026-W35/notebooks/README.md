# Notebook links

The only notebooks supporting the Week 35 report are:

- `weeks/2026-W35/notebooks/nesso_architecture_visualized.ipynb`
- `weeks/2026-W35/notebooks/boltz2_nesso1_rnp_postcutoff50.ipynb`

The paired Boltz-2/Nesso-1 notebook is the presentation-ready result summary
for this report. It reads compact tracked outputs and embeds both the direct
familiarity-versus-accuracy comparison and paired model comparison without
rerunning inference. The older Nesso-only notebooks remain as experiment
provenance but are not used in the current report.

The architecture companion starts with a labeled miniature calculation and
then loads one real saved Nesso-1 checkpoint trace. It visualizes the separate
Nesso input-token and ESM-2 pathways, initial pair tensor, triangle updates,
pair attention, 64-bin distograms, pocket selection, affinity pooling, and the
two continuous affinity outputs. Its 13 Plotly figures are saved in the
executed notebook and remain interactive when opened; Nesso inference does not
need to be rerun. The miniature has random teaching weights and is explicitly
not a scientific predictor.

The durable source for the companion is split across:

- `scripts/build_nesso_architecture_notebook.py`: notebook builder;
- `src/affinity_benchmark/educational/mini_nesso.py`: transparent tensor
  operations; and
- `tests/test_mini_nesso.py`: focused numerical checks.

The Boltz-2/Nesso-1 notebook is an executed companion to the paired 50-system
post-cutoff Runs N' Poses experiment. It shows completion and runtime, the
49-system paired distance/contact/pocket comparison, model-specific familiarity
correlations, and Boltz-only protein and ligand RMSD results. Its durable result
record is `reports/exp015_boltz2_nesso1_rnp_postcutoff50/`.

Open it with the working Anaconda Jupyter installation on this machine:

```bash
/home/yusheng/anaconda3/bin/jupyter-notebook \
  weeks/2026-W35/notebooks/
```

Older Nesso-only Runs N' Poses explorations and the RFdiffusion architecture
notebook are staged under `weeks/2026-W36/notebooks/`, not mixed with this
week's final evidence.
