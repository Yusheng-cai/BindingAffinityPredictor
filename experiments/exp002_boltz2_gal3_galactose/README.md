# exp002: Boltz-2 galectin-3–galactose observability test

## Status

**Seed-42 prediction and withheld-reference comparison completed on
2026-08-07; results are pending review before any expansion.**

## Scientific question

Given only the exact 139-residue human galectin-3 carbohydrate-recognition
domain sequence, stereochemically explicit beta-D-galactopyranose, and an MSA,
can pinned Boltz-2 generate an interpretable complex and recover the binding pose
in PDB 9D63?

This experiment also asks an engineering question: which meaningful intermediate
artifacts can be preserved without modifying Boltz's neural-network code?

## Why 9D63

- Wild-type human galectin-3 CRD and beta-D-galactose match the requested
  biological system.
- The 1.15 A X-ray structure provides a high-resolution pose reference.
- It was released on 2026-03-18, after Boltz-2 was developed, substantially
  reducing the risk that the deposited coordinates were used directly in model
  training.
- At 139 protein residues and 12 ligand atoms including oxygens, it is a small
  first test for the available 10 GiB GPU.

The experiment has no matched quantitative affinity measurement. The native
continuous and binary affinity outputs will be recorded, but they cannot be
called accurate or inaccurate for this system.

## Withheld-information boundary

Inference receives sequence, ligand SMILES, and the server-generated MSA. It
receives no experimental coordinates, template, pocket residue, ligand pose, or
contact constraint. The 9D63 coordinate file is downloaded only after the
prediction has completed, for scoring.

## What “detailed outputs” means

The unmodified CLI preserves the raw MSA, parsed records, processed molecular
objects, model-ready structure/MSA arrays, full PAE and PDE arrays, per-atom
pLDDT, single (`s`) and pair (`z`) trunk embeddings, predicted mmCIF,
confidence summary, and native affinity summary. Logs, resolved parameters,
software versions, timings, and GPU-memory observations are retained beside
them.

These outputs expose each scientific pipeline stage, but not every activation
from every neural-network layer or all 200 transient diffusion coordinate
states. Obtaining those would require a separate source-instrumentation variant
and would change performance and storage requirements.

## Stages

1. Validate the input sequence and ligand stereochemistry.
2. Submit only the protein sequence to the ColabFold MMseqs2 server and preserve
   the returned MSA.
3. Parse and tokenize the protein and ligand; preserve Boltz's processed arrays.
4. Featurize the tokens, atoms, bonds, MSA, and masks.
5. Run three trunk recycling iterations for structure prediction and write the
   final single and pair embeddings. At this pinned revision the separate
   affinity path uses five recycling iterations internally.
6. Run one 200-step structure diffusion sample.
7. Run the confidence head and write pLDDT, PAE, PDE, and scalar summaries.
8. Run the affinity pathway with five 200-step affinity diffusion samples and
   write both native outputs.
9. Only then download 9D63 coordinates and calculate pose diagnostics.

The run stops after seed 42 for review.

## Resolved run details

The structure pass uses the requested three recycles. Inspection of the pinned
CLI source and emitted affinity checkpoint hyperparameters shows that the
separate affinity pass fixes its recycling count at five. The MSA module used
all 6,387 parsed sequences; MSA subsampling resolved to false for this run.

The compact reviewed result is in
[`reports/exp002_boltz2_gal3_galactose/report.md`](../../reports/exp002_boltz2_gal3_galactose/report.md).
Raw predictions, intermediate arrays, logs, and the reference coordinate file
remain ignored by Git under `runs/` and `data/raw/`.
