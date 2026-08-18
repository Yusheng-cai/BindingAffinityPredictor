# FlashBind released-pose FEP+4 result

The released two-checkpoint FlashBind value ensemble completed all 87 records
using the authors' precomputed FABind+ poses and representations.

## Primary result

- Compound-weighted mean of per-target Pearson correlations: **0.5333**
- Compound-bootstrap 95% interval: **0.3821–0.6693**
- Paper reference: **0.53**
- Complete predictions: **87/87**

Per-target Pearson correlations were 0.612 for CDK2, 0.134 for TYK2, 0.546
for JNK1, and 0.677 for p38. The complete metric record, including Spearman,
Kendall, and centered-error results, is in
`flashbind_metrics.json`.

This is a released-pose scoring reproduction, not independent end-to-end
validation. The receptor structures, FABind+ poses, ESM3 representations, and
TorchDrug features came from the authors' released archive. A subsequent
experiment would need to rerun FABind+ from receptor PDBs and ligand SMILES to
test docking reproducibility separately.
