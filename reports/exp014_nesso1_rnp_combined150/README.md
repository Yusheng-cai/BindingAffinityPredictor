# Nesso-1 structural familiarity across 150 Runs N' Poses systems

## Result

We pooled all 150 scored systems into one analysis. Nesso-1 structural accuracy
depends strongly on similarity to pre-September-2021 protein-pocket and
ligand-pose examples.

| Relationship across 150 systems | Spearman rho | Cluster-bootstrap 95% CI |
| --- | ---: | ---: |
| Familiarity vs distogram NLL | -0.651 | [-0.749, -0.529] |
| Familiarity vs distance MAE | -0.686 | [-0.781, -0.571] |
| Familiarity vs contact AP | +0.649 | [+0.517, +0.752] |
| Familiarity vs contact AUROC | +0.625 | [+0.491, +0.736] |
| Familiarity vs pocket F1 | +0.638 | [+0.508, +0.734] |
| Nesso entropy vs distance MAE | +0.785 | [+0.685, +0.859] |

Negative correlations are favorable for NLL and MAE because these are errors.
Positive correlations are favorable for contact and pocket metrics. Every
interval excludes zero.

## Performance across familiarity bins

| Familiarity | N | Distogram NLL | Distance MAE (A) | Contact AP | Pocket F1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0--20 | 19 | 5.886 | 5.833 | 0.048 | 0.000 |
| 20--30 | 19 | 4.126 | 3.067 | 0.229 | 0.235 |
| 30--40 | 19 | 3.447 | 2.208 | 0.500 | 0.350 |
| 40--50 | 20 | 3.235 | 1.600 | 0.532 | 0.417 |
| 50--60 | 19 | 2.874 | 0.805 | 0.877 | 0.720 |
| 60--70 | 18 | 1.852 | 0.492 | 0.932 | 0.700 |
| 70--80 | 18 | 1.924 | 0.505 | 0.944 | 0.686 |
| 80--100 | 18 | 1.877 | 0.471 | 0.959 | 0.683 |

These are bin medians. The principal result is the continuous rank correlation,
not a claim that every successive bin must improve.

## Interpretation

For highly unfamiliar systems, Nesso often assigns little probability to the
experimental protein--ligand distance pattern and frequently fails to recover
the 6 A pocket. For familiar systems, median expected-distance error falls
below 1 A and contact ranking becomes strong.

Nesso's native protein--ligand entropy is also informative: higher entropy is
associated with greater distance error. Familiarity and entropy therefore
provide related but distinct warnings about structural reliability.

This remains a structural-representation test, not an affinity benchmark.
Runs N' Poses does not provide a systematic experimental affinity label for
these systems.

## Provenance

- 150/150 Nesso predictions and 150/150 coordinate comparisons succeeded.
- The model received protein sequence and ligand SMILES only; experimental
  coordinates entered only during scoring.
- Model revision v1.0.0, inference seed 42, five recycling steps.
- The 150 systems contain 150 unique PDB entries, 136 structural clusters, and
  148 unique ligand SMILES.
- Confidence intervals use 2,000 cluster-level bootstrap replicates.
- Compact outputs: [`results/summary.json`](results/summary.json) and
  [`results/per_system_metrics.csv`](results/per_system_metrics.csv).
