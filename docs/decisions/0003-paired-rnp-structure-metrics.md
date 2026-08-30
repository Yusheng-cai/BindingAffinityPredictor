# Decision 0003: Paired structural metrics for Nesso-1 and Boltz-2

Date: 2026-08-28

Status: accepted for `exp015_boltz2_nesso1_rnp_postcutoff50`

## Context

Nesso-1 returns protein-ligand distance distributions, whereas Boltz-2 returns
an explicit coordinate model. RMSD is therefore not a valid common metric.
Both outputs can nevertheless be reduced to protein-token--ligand-heavy-atom
distances and residue-pocket labels. Ligand graph automorphisms create multiple
chemically equivalent atom correspondences, so a numbering-dependent mapping
would make the paired comparison arbitrary.

The earlier Nesso-only `exp012`--`exp014` protocol chose one symmetry mapping by
minimum distogram NLL and used it for every metric. That historical protocol is
retained as the default `legacy_nll` mode of the Nesso scorer so its reviewed
results remain exactly reproducible.

## Decision

- The common protein token is C-beta, or C-alpha for glycine; ligand tokens are
  heavy atoms.
- Interface distance MAE is evaluated for resolved experimental token pairs at
  distance at most 15 A.
- Token-contact precision, recall, and F1 use a 6 A token distance.
- Physical pocket precision, recall, and F1 use residues with any heavy atom at
  most 6 A from any ligand heavy atom.
- Ligand correspondence is restricted to exact graph-isomorphic mappings.
- For `exp015`, interface MAE is minimized and contact F1 is maximized over
  those mappings independently. Ties use the lexicographically first mapping.
  The same rule is applied to both models.
- Nesso uses expected distances from its 64-bin distribution. Boltz uses
  distances in its single sampled coordinate model. This is the closest common
  representation, not an assertion that the native outputs are identical.
- Boltz-only protein and ligand RMSD remain additional metrics. Ligand RMSD is
  symmetry-corrected after a fixed global protein C-alpha alignment; ligand
  atoms do not participate in the fit.
- The primary paired comparison includes only systems completed by both models.
  Model completion/failure rates are reported on all 50 attempted systems.

## Consequences

The paired metrics are invariant to arbitrary numbering of chemically
equivalent ligand atoms and do not invent Nesso coordinates. Pocket F1 favors
models capable of predicting correct side-chain geometry, so it is an
output-capability comparison as well as a pocket-localization comparison. The
Nesso scorer requires `--symmetry-policy metric_specific` to reproduce exp015;
omitting it deliberately preserves the earlier Nesso-only protocol.
