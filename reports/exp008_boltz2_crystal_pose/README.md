# Boltz-2 pose recovery on 16 FEP+4 cocrystals

## Result

The existing seed-42 Boltz-2 MSA-1024 predictions were compared with all 16
exact FEP+4 protein--ligand identities for which an experimental cocrystal was
identified. All 16 comparisons completed.

The median symmetry-corrected ligand heavy-atom RMSD after local pocket
alignment was **0.898 Angstrom** (mean 0.977 Angstrom; range 0.181--2.790
Angstrom). **14/16 poses (87.5%) were below the descriptive 2-Angstrom
threshold.** Median 4-Angstrom native-contact precision and recall were 0.882
and 0.904, respectively.

The two poses above 2 Angstrom were CDK2--N76/PDB 1OIU (2.152 Angstrom) and
p38--533/PDB 3FMH (2.790 Angstrom). Their ligand-only internal-conformation
RMSDs were 1.974 and 2.223 Angstrom, indicating that ligand conformation makes
a substantial contribution to these errors.

## What was aligned

The prediction was not superposed using the ligand. For each pair:

1. the predicted and experimental kinase chains were matched by sequence;
2. a Kabsch transform was calculated from protein C-alpha atoms only;
3. that fixed transform was applied to the predicted ligand;
4. ligand heavy atoms were matched by chirality-aware graph isomorphism;
5. RMSD was minimized only over chemically symmetry-equivalent atom mappings,
   without any ligand-only rigid-body fit.

The primary result uses C-alpha atoms from crystallographic pocket residues,
defined as residues with any heavy atom within 5 Angstrom of the experimental
ligand. A secondary global result uses every sequence-matched kinase C-alpha
atom. The median global-aligned ligand RMSD was 0.697 Angstrom.

## Per-complex results

| Sample | PDB | Pocket-aligned ligand RMSD (A) | Global-aligned ligand RMSD (A) | Internal ligand RMSD (A) | 4-A contact Jaccard |
|---|---:|---:|---:|---:|---:|
| cdk2_1h1q | 1H1Q | 0.722 | 0.739 | 0.617 | 0.895 |
| cdk2_1h1r | 1H1R | 1.502 | 1.465 | 1.379 | 0.737 |
| cdk2_1h1s | 1H1S | 0.313 | 0.398 | 0.224 | 0.947 |
| cdk2_1oi9 | 1OI9 | 0.854 | 0.774 | 0.623 | 0.722 |
| cdk2_1oiu | 1OIU | **2.152** | 2.144 | 1.974 | 0.842 |
| cdk2_1oiy | 1OIY | 0.297 | 0.326 | 0.236 | 0.944 |
| tyk2_ejm_46 | 4GIH | 0.181 | 0.285 | 0.133 | 1.000 |
| tyk2_jmc_23 | 4GJ2 | 0.296 | 0.327 | 0.240 | 0.944 |
| jnk1_17124-1 | 2GMX | 0.502 | 0.480 | 0.445 | 0.875 |
| p38_p38a_3fln | 3FLN | 1.137 | 0.654 | 0.162 | 0.750 |
| p38_p38a_3flq | 3FLQ | 1.108 | 1.220 | 0.967 | 0.684 |
| p38_p38a_3flw | 3FLW | 1.111 | 0.483 | 0.319 | 0.682 |
| p38_p38a_3fly | 3FLY | 1.175 | 0.255 | 0.184 | 0.778 |
| p38_p38a_3flz | 3FLZ | 0.547 | 0.767 | 0.339 | 0.875 |
| p38_p38a_3fmh | 3FMH | **2.790** | 2.423 | 2.223 | 0.684 |
| p38_p38a_3fmk | 3FMK | 0.941 | 0.882 | 0.702 | 0.727 |

## Interpretation limits

These are old, publicly deposited structures and may overlap Boltz-2 training
data. The calculation is therefore a strong technical pose-reproduction sanity
check, but not a prospective or time-split generalization benchmark. In
addition, only compounds with deposited structures enter this subset, creating
selection bias relative to the complete 87-compound affinity benchmark.

The calculation used one stochastic Boltz structure sample and one seed per
pair. It does not measure sampling variability or top-k pose success.

## Reproduction

The protocol is frozen in
`configs/benchmarks/fepplus4_crystal16_pose.yaml`; exact PDB/CCD identities are
in `data/manifests/fepplus4_crystal16.json`. Raw RCSB files and their checksums
are local under `data/raw/exp008_boltz2_crystal_pose/`. Raw per-atom mappings and
metrics are under `runs/exp008_boltz2_crystal_pose/seed42/`.

```bash
PYTHONPATH=src /home/yusheng/programs/boltz-2.2.1/bin/python \
  scripts/fetch_fepplus4_crystal_references.py

PYTHONPATH=src /home/yusheng/programs/boltz-2.2.1/bin/python \
  scripts/analyze_boltz2_crystal_pose.py
```
