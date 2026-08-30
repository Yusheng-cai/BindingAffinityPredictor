# exp013: Nesso-1 Runs N' Poses confirmation cohort

## Objective

Run Nesso-1 on 50 additional Runs N' Poses systems and determine whether the
familiarity trend found in `exp012` reproduces in records that were not part of
the original 100-system pilot.

## Frozen selection

The selection excludes the exact 100 sample IDs in
[`rnp_nesso1_pilot100.json`](../../data/manifests/rnp_nesso1_pilot100.json),
then applies the same Nesso-compatible eligibility filters to the parent Runs
N' Poses collection. It does not require the new systems to belong to clusters
that were absent from the original cohort. Within the new cohort, clusters,
PDB entries, and canonical ligand SMILES are unique.

Across the original and additional cohorts, 14 structural clusters and two
ligand SMILES recur. No exact sample ID or PDB entry recurs. The pooled
150-record uncertainty analysis therefore resamples complete Runs N' Poses
clusters.

The 50 systems are balanced across the paper's familiarity intervals: six per
interval, plus one additional system in each of the 40--50 and 50--60 bins to
examine the previously observed change around 50. The manifest was frozen
before any inference or reference-coordinate inspection:
[`rnp_nesso1_confirmation50.json`](../../data/manifests/rnp_nesso1_confirmation50.json).

## Protocol

Nesso receives only the exact deposited protein sequence and stereochemical
ligand SMILES. It receives no experimental structure, pocket, MSA, or template.
Inference and coordinate-free scoring are unchanged from `exp012`: model
revision v1.0.0, seed 42, five recycling steps, and the same distogram, contact,
pocket, and entropy metrics.

Results will be reported for the original 100, the new 50, and the combined
150. A 40--50 versus 50--60 comparison is exploratory and will not be presented
as evidence of a universal threshold without replication.

## Result

All 50 Nesso runs and all 50 coordinate comparisons completed without failure.
The additional cohort reproduces the continuous familiarity relationship:
distance MAE decreases with familiarity (Spearman rho -0.630), while contact AP
(rho +0.529) and pocket F1 (rho +0.625) increase. Nesso's native interface
entropy remains associated with distance MAE (rho +0.701).

The sharp change around score 50 seen in the first 100 does not reproduce for
all metrics. Only pocket F1 improves consistently from the 40--50 to the 50--60
bin in both cohorts. The stronger conclusion is therefore a continuous
familiarity dependence, not a universal threshold at 50. See the
[reviewed result](../../reports/exp013_nesso1_rnp_confirmation50/README.md).
