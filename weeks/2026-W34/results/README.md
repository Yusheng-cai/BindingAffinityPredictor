# Compact reviewed results

These CSV files are small presentation snapshots of the reviewed values in
the weekly report. They do not replace the canonical experiment reports or raw
model outputs.

- `benchmark_summary.csv` contains aggregate affinity metrics for the common
  87-compound benchmark.
- `pose_summary.csv` contains overall pose-recovery statistics on the 16 exact
  cocrystal matches.
- `pose_rmsd_by_target.csv` contains the target-level medians plotted in the
  weekly report.

Canonical metric files remain under `reports/exp007_*`,
`reports/exp009_*`, and `reports/exp010_*`. Raw native model outputs and run
records remain under the Git-ignored `runs/` tree.
