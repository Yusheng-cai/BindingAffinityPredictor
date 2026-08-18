# exp003 result: one RFdiffusion backbone and ProteinMPNN handoff

## Outcome

The frozen engineering smoke test passed. RFdiffusion produced one 90-residue
binder backbone against the bundled 150-residue insulin-receptor target, and
ProteinMPNN produced four sequences for the generated chain while keeping the
target fixed.

This result establishes that the pinned local pipeline executes and produces
the expected representations. It does not establish that the candidate folds
or binds.

## Audited outputs

| Check | Observed | Result |
|---|---:|---|
| GPU execution | NVIDIA GeForce RTX 3080 | pass |
| Target residues | 150 | pass |
| Binder residues | 90 (requested 70–100) | pass |
| `Xt-1` coordinate blocks | 50 | pass |
| `pX0` coordinate blocks | 50 | pass |
| ProteinMPNN candidates | 4 | pass |
| Designed chain | B only | pass |

RFdiffusion reported 66.01 s for the design; command wall time including model
and schedule initialization was 79.32 s. ProteinMPNN required 4.46 s.

## ProteinMPNN candidates

| Sample | Designed-chain score ↓ | Global score | Frozen selection |
|---:|---:|---:|---|
| 1 | 0.7786 | 1.2917 | selected |
| 2 | 0.8071 | 1.3053 | — |
| 3 | 0.7828 | 1.2927 | — |
| 4 | 0.8508 | 1.3264 | — |

These are negative-log-probability sequence-model scores. They are not energies,
affinities, confidence in experimental success, or values comparable across
different ProteinMPNN protocols without additional controls.

## Evidence boundary and next decision

The next scientifically distinct step would be an independent sequence-based
complex prediction and prespecified interface analysis. It should only be run
after choosing the prediction protocol—including whether a novel designed
sequence may be submitted to an external MSA service. Experimental evidence
would still be required after that computational filter.
