# exp012: Nesso-1 on a diverse Runs N' Poses pilot

## Objective

Test whether Nesso-1's coordinate-free structural representation generalizes
as protein-pocket and ligand-pose familiarity decreases. Nesso-1 does not write
a predicted all-atom complex; the primary object evaluated here is its
protein--ligand distance distribution (`pdistogram`). Pocket recovery and
structural uncertainty are evaluated without first reconstructing Cartesian
coordinates.

## Frozen pilot

The 100 systems were selected from the official Runs N' Poses annotations
before any Nesso inference. Selection uses the paper's eight
`sucos_shape_pocket_qcov` intervals. The four lowest-familiarity intervals
contain 13 systems each and the remaining intervals contain 12 each. Every
selected system has a unique Runs N' Poses cluster, PDB entry, and canonical
ligand SMILES.

Eligibility was deliberately simple for the first implementation: exactly one
protein chain and one proper ligand, 40--900 standard amino acids, 6--50 ligand
heavy atoms, a valid sequence/SMILES input record, and a non-null September
2021 familiarity score. Of 4,235 annotation rows, 650 passed all filters.

The canonical selection is
[`data/manifests/rnp_nesso1_pilot100.json`](../../data/manifests/rnp_nesso1_pilot100.json).
The selection algorithm is
[`scripts/select_rnp_nesso_pilot.py`](../../scripts/select_rnp_nesso_pilot.py).

## Evidence boundary

The official `annotations.csv`, `inputs.json`, and `ground_truth.tar.gz` were
downloaded after user approval. Only the 100 systems in the frozen manifest
were extracted from the ground-truth archive. Released prediction archives,
MSAs, and the full all-pairs similarity archive were not downloaded. The
per-sample `ground_truth_coordinates_downloaded: false` values in the manifest
record the state when selection was frozen; the immutable selection file was
not rewritten after the approved download.

Runs N' Poses has no systematic affinity labels, so this experiment cannot
measure affinity Pearson correlation or affinity calibration.

## One-system gate

The lowest-familiarity system, `7ftm__1__1.A__1.C__1.C` (similarity score
0.0), passed the inference and mapping gate with no model or scoring failure.
The 362-residue input maps exactly to the deposited construct; 356 Nesso
protein distance atoms are resolved. The 17-heavy-atom ligand has two exact
graph-isomorphic mappings to the SDF, and an independently calculated
any-heavy-atom 6 A pocket contains 13 residues, matching the Runs N' Poses
annotation.

This hard example is a candidate model failure rather than a pipeline failure:
its 6 A physical-pocket F1 is 0.0 and its token-contact AUROC is 0.451. The full
frozen panel was therefore launched unchanged to test whether performance
improves with structural familiarity.

The native Nesso affinity pocket mask uses a 15 A cutoff. It is retained as a
native output but is not mislabeled as a 6 A pocket prediction. For the 6 A
benchmark, the predicted residue pocket is recomputed from the minimum expected
protein-token--ligand distance at a 6 A cutoff.

## Full result

All 100 model runs and all 100 final coordinate comparisons completed. The
central finding is a strong dependence on pre-cutoff structural familiarity:
distogram NLL decreases with familiarity (Spearman rho -0.704; 95% cluster-
bootstrap CI [-0.796, -0.577]), while physical-pocket F1 increases (rho +0.647;
CI [+0.503, +0.760]). Nesso's native protein--ligand entropy also tracks
expected-distance MAE (rho +0.824; CI [+0.707, +0.903]).

See the [reviewed result](../../reports/exp012_nesso1_rnp_distogram_generalization/README.md),
[compact summary](../../reports/exp012_nesso1_rnp_distogram_generalization/results/summary.json),
and [per-system table](../../reports/exp012_nesso1_rnp_distogram_generalization/results/per_system_metrics.csv).
