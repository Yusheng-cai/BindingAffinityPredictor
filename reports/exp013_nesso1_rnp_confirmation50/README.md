# Nesso-1 on 50 additional Runs N' Poses systems

## What was tested

We selected 50 additional protein--ligand systems from the Runs N' Poses
collection, excluding the exact 100 records used previously. The new records
were balanced across the eight structural-familiarity bins, with seven rather
than six systems in both the 40--50 and 50--60 bins. Within the new 50, every
cluster, PDB entry, and ligand SMILES is unique. Across the two cohorts, 14
structural clusters and two ligand SMILES recur, but no exact sample or PDB
entry does. Combined confidence intervals therefore resample whole structural
clusters.

Nesso-1 received only the protein sequence and ligand SMILES. Experimental
coordinates were withheld until scoring. Model revision, seed, recycling,
cropping, and all metric definitions were unchanged from the original
100-system experiment.

## Main result

The broad familiarity trend reproduced in the non-overlapping set of records.

| Relationship | Original 100 rho (95% CI) | New 50 rho (95% CI) | Combined 150 rho (95% CI) |
| --- | ---: | ---: | ---: |
| Familiarity vs distogram NLL | -0.704 [-0.796, -0.577] | -0.513 [-0.736, -0.260] | -0.651 [-0.749, -0.529] |
| Familiarity vs distance MAE | -0.709 [-0.814, -0.584] | -0.630 [-0.806, -0.388] | -0.686 [-0.781, -0.571] |
| Familiarity vs contact AP | +0.682 [+0.543, +0.787] | +0.529 [+0.258, +0.738] | +0.649 [+0.517, +0.752] |
| Familiarity vs pocket F1 | +0.647 [+0.503, +0.760] | +0.625 [+0.410, +0.775] | +0.638 [+0.508, +0.734] |
| Nesso entropy vs distance MAE | +0.824 [+0.710, +0.901] | +0.701 [+0.486, +0.846] | +0.785 [+0.685, +0.859] |

Here, negative correlations are favorable for NLL and MAE because lower is
better. Positive correlations are favorable for contact AP and pocket F1
because higher is better. Every confirmation-cohort interval excludes zero.

## Was there a transition at familiarity 50?

Not consistently. In the original 100, the median 50--60 bin improved over the
40--50 bin for all four highlighted structural metrics. In the new 50:

| Metric | New-50 median at 40--50 | New-50 median at 50--60 | Replicated sharp improvement? |
| --- | ---: | ---: | --- |
| Distogram NLL | 2.237 | 3.203 | No |
| Distance MAE (A) | 0.876 | 0.907 | No; essentially unchanged |
| Contact AP | 0.905 | 0.779 | No |
| Pocket F1 | 0.414 | 0.760 | Yes |

Thus, the scientifically supported statement is that Nesso's structural
representation performs better on more familiar systems **across the full
similarity range**. The current data do not support treating 50 as a universal
changepoint. Pocket recovery may have a particularly strong change in this
region, but that is a narrower hypothesis for a future, more densely sampled
test.

## Reproducibility

- Inference: 50/50 successful; scoring: 50/50 successful.
- Nesso version/model revision: v1.0.0; inference seed: 42.
- Wall time: 401.9 s, including preprocessing and output writing.
- Peak GPU allocation observed system-wide: 9,294 MiB on the RTX 3080.
- Frozen manifest SHA-256: `14b03457edfa2d93e90b24334729ac86d8a01d05e270b6bb7e6280400634bc26`.
- Compact results: [`results/summary.json`](results/summary.json),
  [`results/comparison.json`](results/comparison.json), and
  [`results/combined_per_system_metrics.csv`](results/combined_per_system_metrics.csv).

Raw predictions remain under the ignored `runs/` tree and are not committed.
