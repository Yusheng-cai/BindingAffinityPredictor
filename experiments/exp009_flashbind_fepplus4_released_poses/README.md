# exp009: FlashBind on FEP+4 with released FABind+ poses

## Scientific objective

Test whether the two released FlashBind value checkpoints reproduce the
reported within-target affinity ranking on exactly the same 87 neutral
compounds used for the Boltz-2 and Nesso-1 comparison.

This first experiment deliberately uses the authors' precomputed FABind+
ligand poses, protein and ligand representations, and fixed receptor PDBs. It
therefore tests the released FlashBind scoring stage and released preprocessing
as a unit. It does not independently test whether a fresh FABind+ run would
recover the same pockets or poses.

## Identity and label validation

All 87 released records match the canonical manifest by target and InChIKey:
16 CDK2, 16 TYK2, 21 JNK1, and 34 p38 compounds. The released labels equal the
canonical `log10_value_uM` values to floating-point precision, and all four
receptor PDBs are byte-identical to the local benchmark copies. Metrics are
nevertheless calculated from the canonical manifest, not the released label
file.

## Execution order

1. Run one CDK2 record through both checkpoints and validate all native fields.
2. Run all 87 records without changing the model or cropping settings.
3. Preserve each checkpoint's JSON and the arithmetic-mean ensemble JSON.
4. Convert the ensemble output to the canonical prediction table.
5. Calculate per-target Pearson, Spearman, Kendall, pairwise MAE, ordinary MAE,
   and centered MAE, then aggregate target metrics by compound count.
6. Bootstrap compounds within targets with the same fixed policy used for
   Boltz-2 and Nesso-1.

The preregistered paper reference is weighted Pearson R = 0.53 and weighted
Kendall tau = 0.38. These are comparison values, not acceptance thresholds.
