# Input Manifests

A canonical affinity manifest will normally include:

- sample, target, ligand, and assay identifiers;
- exact protein construct sequence;
- isomeric SMILES and ligand provenance;
- measurement type, value, units, and qualifier;
- source dataset/version and source record;
- optional experimental structure and receptor-template identifiers;
- split, inclusion status, and exclusion reason.

The initial JSON schema is implemented in
`src/affinity_benchmark/data/manifest.py`. It intentionally checks only facts
that can be validated without chemistry or structure dependencies. RDKit-based
standardization and coordinate-level checks will be added before the first
multi-ligand experiment.
