# Boltz-2: what is being tested

## Model in one sentence

Boltz-2 is a full-atom biomolecular cofolding model with a separate learned
affinity module: from molecular identities (for this experiment, one protein
sequence and one ligand SMILES), it samples a three-dimensional complex and
reports structure-confidence quantities plus two distinct affinity outputs.

The primary sources are the [Boltz-2 technical
report](https://doi.org/10.1101/2025.06.14.659707), the [official source
repository](https://github.com/jwohlwend/boltz), and the [official prediction
documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md).

## What “full-atom cofolding” means here

“Cofolding” means that the protein conformation, ligand conformation, and their
relative placement are predicted jointly, rather than docking a ligand into one
rigid, already supplied receptor structure. “Full atom” means that the generated
complex is represented at atomic resolution, including ligand atoms and protein
side-chain atoms, not merely a residue-level backbone or binding score.

This does **not** mean that Boltz-2 explicitly integrates Newtonian dynamics,
solvent trajectories, or a thermodynamic partition function. Its coordinates
are samples from a learned generative model. Its affinity scalar is a learned
prediction, not a rigorous alchemical binding free energy.

## High-level computation

1. The input featurization represents polymer sequence/MSA information and the
   ligand chemical graph.
2. A trunk derived from the Boltz-1/AlphaFold-3 family builds single- and
   pairwise representations through PairFormer-style processing.
3. A diffusion model samples atomic coordinates for the complete complex.
4. Confidence heads assess the predicted structure.
5. The affinity module consumes learned pair representations together with
   predicted structural information and produces continuous-affinity and binary
   binder outputs.

The architecture should be read as coupled structure and property prediction,
not as a physical energy decomposition. AlphaFold 3 is the key architectural
background ([Abramson et al., 2024](https://doi.org/10.1038/s41586-024-07487-w));
Boltz-1 documents the open predecessor ([Wohlwend et al.,
2024](https://doi.org/10.1101/2024.11.19.624167)).

## The two affinity outputs are different tasks

| Native field | Intended interpretation | Direction |
| --- | --- | --- |
| `affinity_pred_value` | Approximate `log10(IC50 / uM)` continuous value for comparing binder potency, especially within related compounds | Lower is stronger |
| `affinity_probability_binary` | Probability-like binder-versus-decoy classification output | Higher is more binder-like |

They are trained with substantially different supervision and must never be
merged into one score. In exp001, the experimental value is a (K_i), not an
(IC_{50}). We retain both labels and show them side by side; numerical
difference is descriptive, not a calibrated error estimate.

For the reference (K_i=0.0048\,\mu\mathrm{M}):

\[
\log_{10}(K_i/\mu\mathrm{M})=-2.31876,
\qquad
pK_i=-\log_{10}(K_i/\mathrm{M})=8.31876.
\]

These are the same concentration on different logarithmic scales:
(pK_i=6-\log_{10}(K_i/\mu\mathrm{M})). They are not interchangeable without
the sign and six-unit offset.

## What exp001 can and cannot establish

It can test whether the pinned software accepts the intended inputs, produces a
complex and both affinity heads, records provenance, recovers a known pose, and
fits local hardware. A ligand heavy-atom RMSD at or below 2 Å is recorded as the
conventional pose-recovery heuristic, but physical validity and contact recovery
must also be checked.

It cannot establish affinity accuracy, ranking, calibration, selectivity, or
generalization. Correlation is undefined for one observation. Moreover, 4GIH is
an old public structure and is not a novelty-controlled test for a modern model.

## Frozen exp001 choices

- Exact deposited 4GIH entity sequence plus neutral, achiral 0X5 SMILES.
- No holo template, receptor coordinates, pocket residues, or contact restraints.
- MSA generated through the configured public server only after explicit approval.
- One structure sample and seed 42 for the first resource check.
- Upstream default structure settings, with `use_potentials: false`.
- Five affinity diffusion samples, with molecular-weight correction disabled.
- Crystal 4GIH coordinates withheld until post-prediction scoring.

Changing potentials or the molecular-weight correction defines a different
method and therefore requires a separate protocol variant.
