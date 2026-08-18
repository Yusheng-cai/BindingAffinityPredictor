# exp004 result: RFdiffusion structural diversity

## Outcome

Ten deterministic seed samples from the same target-conditioned RFdiffusion
protocol produced materially different binder backbone shapes and target-relative
poses. All ten runs passed the frozen engineering checks; none has been shown to
fold or bind.

## Ensemble summary

| Quantity | Observed |
|---|---:|
| Successful / attempted backbones | 10 / 10 |
| Binder length | 71–98 residues; mean 82.2 |
| Cα radius of gyration | 12.01–15.82 Å |
| Target contact residues within 10 Å Cα distance | 20–35 |
| Final binder pLDDT-like output | 0.9909–0.9929 |
| Intrinsic trace RMSD across 45 pairs | median 13.36 Å; range 4.68–16.45 Å |
| Target-aligned pose Chamfer distance | median 5.08 Å; range 2.66–10.40 Å |
| Binder center displacement | median 6.37 Å; range 0.54–14.20 Å |
| Contact-set Jaccard similarity | median 0.50; range 0.36–0.78 |
| Designs with interchain backbone pairs below 2 Å | 0 |

The closest pair by intrinsic trace shape was seeds 42/48 at 4.68 Å; the most
different was 42/43 at 16.45 Å. The closest target-relative pose pair was 49/50
at 2.66 Å Chamfer distance; the most different was 43/50 at 10.40 Å. Different
metrics answer different geometrical questions, so these pair labels are not a
single ranking of “best” designs.

## Interpretation

The main scientific lesson is that one RFdiffusion output is one stochastic
sample. Identical conditioning does not determine a unique binder fold, length,
or interface footprint. The near-constant final RFdiffusion confidence values
also show why confidence alone cannot summarize ensemble diversity.

These results do **not** imply ten independent binding predictions. RFdiffusion
generated poly-glycine backbone proposals. Only seed 42 has undergone the
ProteinMPNN handoff, and none has undergone an independent sequence-based
complex prediction or an experiment.

## Metric cautions

- Intrinsic trace RMSD resamples unequal-length Cα curves to 64 arc-length
  positions before Kabsch alignment. It is a transparent descriptive metric,
  not TM-score and not a claim of fold homology.
- Pose Chamfer is a symmetric nearest-neighbor Cα distance after aligning the
  target. It tolerates unequal lengths but has no direct energetic meaning.
- A target contact is defined with a coarse 10 Å Cα cutoff. Contact-set overlap
  does not establish atomic complementarity or favorable thermodynamics.
- The clash screen uses only backbone N/Cα/C/O atoms because the RFdiffusion
  binder outputs are glycine placeholders.

See [per_design.csv](per_design.csv), [pairwise.csv](pairwise.csv), and
[summary.json](summary.json) for the complete compact outputs.
