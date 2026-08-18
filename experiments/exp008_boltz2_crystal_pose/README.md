# exp008: Boltz-2 pose recovery on the FEP+4 crystal subset

## Scientific objective

Measure how closely the existing seed-42 Boltz-2 predictions reproduce the
experimental ligand poses for the 16 FEP+4 protein--ligand identities with an
exact deposited cocrystal structure.

This is a retrospective structural reproduction test. The reference coordinates
were withheld from Boltz-2 inference, but the PDB entries predate the model and
may overlap its training data. The result therefore does not establish
generalization to unseen protein--ligand chemistry.

## Frozen comparison

The kinase chain is mapped between prediction and reference by sequence. Two
protein-only Kabsch fits are calculated:

1. all sequence-mapped kinase C-alpha atoms;
2. sequence-mapped C-alpha atoms belonging to reference residues with any heavy
   atom within 5 Angstrom of the crystallographic ligand.

Each protein-derived transform is applied unchanged to the predicted ligand.
Ligand heavy atoms are matched by chirality-aware graph isomorphism, and the
minimum RMSD across symmetry-equivalent atom mappings is reported. Ligand atoms
never contribute to the protein fit.

The primary metric is pocket-aligned ligand heavy-atom RMSD. Global-aligned
ligand RMSD, protein fit RMSDs, ligand internal-conformation RMSD, and 4-Angstrom
contact recovery are secondary diagnostics.

## Inputs and outputs

- Crystal-pair manifest: `data/manifests/fepplus4_crystal16.json`
- Metric policy: `configs/benchmarks/fepplus4_crystal16_pose.yaml`
- Prediction source: `runs/exp007_fepplus4_boltz2_nesso1/boltz2_msa1024/seed42/`
- Raw RCSB downloads: `data/raw/exp008_boltz2_crystal_pose/` (ignored)
- Raw analysis output: `runs/exp008_boltz2_crystal_pose/` (ignored)
- Reviewed summary: `reports/exp008_boltz2_crystal_pose/`
