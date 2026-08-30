# exp014: Combined 150-system Nesso Runs N' Poses analysis

This analysis pools the frozen 100-system `exp012` manifest and the additional
50-system `exp013` manifest. No Nesso inference was repeated. Each component
manifest was frozen before its respective inference run, and their exact sample
IDs do not overlap.

The combined analysis contains 150 scored systems with no failures. Because 14
Runs N' Poses structural clusters occur in both components, all pooled
uncertainty intervals resample complete clusters rather than individual rows.

See the [combined report](../../reports/exp014_nesso1_rnp_combined150/README.md),
[combined manifest](../../data/manifests/rnp_nesso1_combined150.json), and
[benchmark configuration](../../configs/benchmarks/rnp_nesso1_combined150.yaml).
