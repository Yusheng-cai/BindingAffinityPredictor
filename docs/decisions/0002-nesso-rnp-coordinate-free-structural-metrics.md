# Decision 0002: Coordinate-free Nesso-1 structural metrics on Runs N' Poses

Date: 2026-08-20

Status: accepted for `exp012_nesso1_rnp_distogram_generalization`

## Context

Nesso-1 predicts a categorical inter-token distance distribution but does not
write a Cartesian protein--ligand complex. A coordinate RMSD would therefore
require a separate reconstruction algorithm and would mix reconstruction error
with the representation being tested. Runs N' Poses supplies experimental
complex coordinates and a structural-familiarity score, but not systematic
experimental affinity labels.

## Decision

- Nesso receives only the deposited protein sequence and ligand SMILES.
  Experimental coordinates, pockets, and similarity neighbors are used only
  after inference.
- The protein distogram token is C-beta, except C-alpha for glycine. Every
  ligand heavy atom is one token, matching the released Nesso tokenizer.
- The primary metric is mean categorical negative log likelihood for resolved
  protein-token--ligand-atom pairs whose experimental distance is at most
  15 A. The released 64-bin scheme has 63 linearly spaced boundaries from
  2 to 22 A and uses the two open end bins.
- With refinement enabled, the released writer saves a full-system distogram:
  locally refined pairs are merged back into the initial global prediction.
  All resolved true-interface pairs are scored. A true pocket that Nesso fails
  to retain in its 256-token refinement crop is not silently excluded.
- Token-contact average precision and AUROC use a 6 A experimental
  C-beta/C-alpha--ligand-heavy-atom definition.
- Physical pocket precision, recall, and F1 use the Runs N' Poses definition:
  a protein residue is positive when any resolved protein heavy atom is within
  6 A of any ligand heavy atom. The predicted residue pocket follows Nesso's
  own minimum-expected-token-distance construction, evaluated at 6 A.
- Nesso's saved native affinity pocket mask uses 15 A. It is retained as a
  native output and reported separately; it is not relabeled as a 6 A pocket.
- Ligand atom correspondence is restricted to exact RDKit graph isomorphisms.
  When chemical symmetry gives multiple mappings, the mapping with minimum
  interface distogram NLL is used and the number and chosen mapping are saved.
- Uncertainty intervals resample Runs N' Poses clusters. The frozen pilot has
  one system per cluster, so this is operationally a system-level bootstrap.
  The analysis uses 2,000 replicates and seed 20260820.
- Native affinity values are preserved, but affinity accuracy is not claimed.
  Coordinate RMSD is not calculated in this experiment.
- As an explicitly exploratory calibration check, entropy--error partial
  Spearman correlations are also calculated after rank-residualizing (i)
  familiarity and (ii) familiarity, protein length, and ligand heavy-atom
  count. These were added after the primary analysis and are not a new primary
  endpoint.

## Consequences

The evaluation directly tests Nesso's released structural representation and
its usefulness for pocket localization without inventing coordinates. Pocket
F1 and token-contact metrics answer related but non-identical questions and
must remain separately labeled. Symmetry minimization makes ligand scoring
invariant to arbitrary numbering of chemically equivalent atoms, at the cost
of reporting the best exact symmetry assignment.
