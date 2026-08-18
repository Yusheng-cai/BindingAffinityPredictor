# exp007: Boltz-2 versus Nesso-1 on FEP+4

## Scientific objective

Reproduce the shared 87-neutral-compound, four-kinase affinity benchmark from
the Boltz-2 and Nesso-1 reports using pinned released models and identical
sequence/SMILES inputs. The primary readout is the compound-count-weighted mean
of per-target Pearson correlations, matching the papers' aggregation policy.

This is a test of within-series ranking on public, largely chemically
in-distribution data. It is not a binding-free-energy calculation and must not
be generalized to prospective medicinal chemistry without an independent
out-of-distribution benchmark.

## Dataset identity

We use the older OpenFF Git commit
`da7c3372256446222e424368be38ef3d2b55a67b` because it contains the full
87-example benchmark: 16 CDK2, 16 TYK2, 21 JNK1, and 34 p38 neutral ligands.

Protein inputs use all deposited protein or peptide entities in the curated
complex: CDK2 plus cyclin A2, JNK1 plus the JIP1 peptide, and the single-chain
TYK2 and p38 constructs. Experimental coordinates, ligand poses, pocket
residues, contact restraints, and affinity labels are withheld from inference.

## Execution order

1. Validate the Nesso-1 installation with a one-input GPU smoke test.
2. Generate one Boltz MSA and one Nesso ESM embedding per unique protein chain.
3. Run all 16 TYK2 ligands through both models and validate the output contract.
4. Without changing the frozen protocol, run the remaining 71 compounds.
5. Preserve raw outputs and calculate per-target and weighted metrics with
   within-target compound bootstrap intervals.
6. Report differences from the published values and investigate failures and
   target-level heterogeneity.

Boltz-2 has not released its exact evaluation pipeline or processed affinity
benchmark inputs. The result is therefore a released-model reproduction with
reconstructed preprocessing, not a bitwise reproduction of the paper.

## Interactive analysis

The completed Nesso-1 calculation can be explored without rerunning inference
in `weeks/2026-W34/notebooks/fepplus4_nesso1_analysis.ipynb`. The notebook reads the canonical
manifest and local prediction table, calls the reusable tested metric
implementation, and reproduces the reported target-level and aggregate values.
