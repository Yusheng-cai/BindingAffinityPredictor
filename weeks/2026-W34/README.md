# Week of August 17, 2026: three affinity models on FEP+4

This week asks three practical questions:

1. What do Boltz-2, Nesso-1, and FlashBind take as input and predict?
2. Can their released checkpoints reproduce the reported within-target
   affinity-ranking performance on the same 87-compound, four-kinase
   benchmark?
3. For the 16 compounds with exact cocrystal references, how accurately do
   Boltz-2 and the released FABind+ ligand-position predictions, generated with
   experimental PDB protein structures, recover the experimental ligand pose?

## Headline results

The compound-count-weighted within-target Pearson correlations were 0.807 for
Nesso-1, 0.627 for Boltz-2 with a documented 1,024-row MSA subsample, and 0.533
for FlashBind using the authors' released FABind+ poses and representations.
The corresponding paper values were 0.800, 0.660, and 0.530.

On the 16 exact structure-matched pairs, the median symmetry-corrected ligand
RMSD after protein-pocket alignment was 0.898 Å for Boltz-2 and 0.634 Å for
the released FABind+ poses. These are not equivalent input conditions:
Boltz-2 starts from sequence, ligand chemistry, and an MSA, whereas FlashBind
uses an experimental PDB protein structure and a FABind+-predicted ligand
position and orientation (pose).

## Navigate this week

- [`report/report.pdf`](report/report.pdf): presentation-ready scientific report
- [`report/main.tex`](report/main.tex): editable LaTeX source
- [`technical-report/report.pdf`](technical-report/report.pdf): detailed
  reproducibility and computational-procedure report
- [`technical-report/main.tex`](technical-report/main.tex): editable technical
  report source
- [`code/README.md`](code/README.md): complete setup, inference, analysis, and
  report-building workflow
- [`notebooks/fepplus4_nesso1_analysis.ipynb`](notebooks/fepplus4_nesso1_analysis.ipynb):
  interactive analysis walkthrough
- [`results/benchmark_summary.csv`](results/benchmark_summary.csv): compact
  affinity-metric snapshot
- [`results/pose_summary.csv`](results/pose_summary.csv): compact pose-recovery
  snapshot
- [`manifest.yaml`](manifest.yaml): experiment, model, run, and report
  provenance

Reusable scientific implementations remain in the repository-level `src/`,
`scripts/`, `configs/`, and `experiments/` directories. Raw predictions,
checkpoints, MSAs, downloaded datasets, and external model installations are
intentionally excluded from Git.
