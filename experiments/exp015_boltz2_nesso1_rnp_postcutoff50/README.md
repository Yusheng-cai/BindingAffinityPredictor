# exp015: Paired Boltz-2/Nesso-1 Runs N' Poses benchmark

This experiment compares Boltz-2 v2.2.1 and Nesso-1 v1.0.0 on exactly the same
50 protein-ligand systems. Every complex was released after 1 June 2023. The
set contains ten systems in each Boltz-2 structural-familiarity band (0–20,
20–40, 40–60, 60–80, and 80–100) and has no repeated PDB entry, ligand,
protein sequence, or Runs N' Poses structural cluster.

## What each model receives

Both models receive only the deposited protein construct sequence and the
stereochemical ligand SMILES. Boltz-2 also receives an MSA generated from the
protein sequence through the ColabFold server. Neither model receives the
experimental protein coordinates, ligand coordinates, pocket residues, or a
template structure.

This is a **structure benchmark**, not an affinity benchmark: Runs N' Poses
does not supply a uniform experimental affinity endpoint for these 50 systems.
Nesso-1 produces protein-ligand distance distributions. Boltz-2 produces an
explicit three-dimensional complex. Their shared comparison therefore uses
interface-distance and pocket/contact metrics. Coordinate RMSD is reported as
an additional Boltz-2 metric, not forced onto Nesso-1.

## Why the familiarity score is model-specific

Runs N' Poses measures how similar each test complex is to structures available
before a chosen date. For Boltz-2 we use 1 June 2023; for Nesso-1 we use 30
September 2021. The same physical complexes are scored against two different
historical structure pools, so the familiarity values need not match. Each
model's error is correlated with the score based on its own cutoff.

## Frozen five-system gate

Before the full run, one median-protein-length example from each Boltz-2
familiarity band is processed end to end:

| Band | Sample | Protein length | Ligand heavy atoms |
| --- | --- | ---: | ---: |
| 0–20 | `7ftm__1__1.A__1.C__1.C` | 362 | 17 |
| 20–40 | `8aqf__1__1.A__1.B__1.B` | 323 | 34 |
| 40–60 | `8jmp__1__1.A__1.D__1.D` | 260 | 28 |
| 60–80 | `8haq__1__1.A__1.C__1.C` | 277 | 31 |
| 80–100 | `8iqt__1__1.A__1.B__1.B` | 306 | 30 |

The run expands to all 50 only if both models produce finite, parseable output
for all five systems and the predicted/reference atom mapping is validated.

The source of truth is the [experiment configuration](experiment.yaml), the
[benchmark configuration](../../configs/benchmarks/rnp_boltz2_nesso1_postcutoff50.yaml),
and the [frozen manifest](../../data/manifests/rnp_boltz2_nesso1_postcutoff50.json).

## Outcome

Nesso-1 completed 50/50 systems. Boltz-2 completed 49/50; the 818-residue 8UK6
system exceeded the 10 GiB GPU limit under the frozen MSA-1024 protocol. The
primary paired comparison therefore contains the 49 systems completed by both
models. The failure remains part of the benchmark failure accounting and was
not silently rescued with a different protocol.

See the [reviewed result](../../reports/exp015_boltz2_nesso1_rnp_postcutoff50/README.md).
