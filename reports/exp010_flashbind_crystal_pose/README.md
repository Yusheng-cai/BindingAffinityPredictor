# FlashBind/FABind+ pose recovery on 16 FEP+4 cocrystals

## Result

I compared the authors' 16 released FABind+ ligand poses with the same exact
experimental cocrystals and the same frozen pose metric used for Boltz-2. All
16 comparisons completed.

The FlashBind/FABind+ median symmetry-corrected ligand heavy-atom RMSD after
local protein-pocket alignment was **0.634 Angstrom** (mean 0.651 Angstrom;
range 0.077--1.595 Angstrom). **16/16 poses were below the descriptive
2-Angstrom threshold.** On the same pairs, Boltz-2 MSA-1024 had a median of
0.898 Angstrom and 14/16 poses below 2 Angstrom.

This is not a controlled claim that FlashBind is more accurate than Boltz-2.
The released FlashBind workflow starts from one experimental receptor
structure per target and supplies authors' archived FABind+ poses, whereas the
Boltz-2 structures were cofolded here from sequence, ligand chemistry, and an
MSA. FlashBind had the lower RMSD on 9 pairs and Boltz-2 on 7; the difference
is driven partly by the two Boltz-2 poses above 2 Angstrom rather than a
uniform pairwise advantage.

## Protein and ligand decomposition

The ligand RMSD is only one component of the Boltz-2 cofolding result. Across
the 16 pairs, Boltz-2 had a median fitted global protein C-alpha RMSD of
**1.231 Angstrom**, a median locally fitted pocket C-alpha RMSD of **0.400
Angstrom**, and a median pocket-aligned ligand heavy-atom RMSD of **0.898
Angstrom**. The corresponding FlashBind receptor values were 0.272, 0.208,
and 0.634 Angstrom, respectively.

Those protein numbers have different meanings for the two workflows. Boltz-2
generated its receptor, so its protein RMSDs are prediction errors. FlashBind
was supplied a fixed experimental receptor; its protein RMSDs only measure the
difference between that supplied receptor and each ligand-specific cocrystal.
No joint protein--ligand fit was performed, because the much larger protein
atom set would dominate such a number and ligand coordinates would then help
define the transformation used to evaluate the ligand.

Boltz-2's locally fitted pocket C-alpha median was target-dependent: 0.348
Angstrom for CDK2, 0.243 for TYK2, 0.460 for JNK1, and **2.305 for p38**. The
p38 ligand result therefore coincides with a substantial predicted-pocket
error and should not be interpreted as ligand placement alone.

## What was calculated

For each prediction--crystal pair, I:

1. matched the predicted and experimental kinase chains by sequence;
2. defined the experimental pocket as protein residues with any heavy atom
   within 5 Angstrom of the crystallographic ligand;
3. calculated a Kabsch transform from matched pocket C-alpha atoms only;
4. applied that fixed protein-derived transform to the predicted ligand;
5. enumerated complete bond-order- and chirality-compatible ligand graph
   mappings and reported the minimum heavy-atom RMSD over symmetry-equivalent
   mappings.

The ligand was never refitted to the crystal ligand. Consequently, the primary
RMSD measures the error of the ligand pose relative to the aligned protein
pocket, not merely whether the two isolated ligand conformations can be
superposed.

## Paired results

| Pair | PDB | FlashBind/FABind+ RMSD (A) | Boltz-2 RMSD (A) | Lower RMSD |
|---|---:|---:|---:|---|
| cdk2_1h1q | 1H1Q | 0.566 | 0.722 | FlashBind |
| cdk2_1h1r | 1H1R | 1.595 | 1.502 | Boltz-2 |
| cdk2_1h1s | 1H1S | 0.485 | 0.313 | Boltz-2 |
| cdk2_1oi9 | 1OI9 | 0.185 | 0.854 | FlashBind |
| cdk2_1oiu | 1OIU | 0.398 | 2.152 | FlashBind |
| cdk2_1oiy | 1OIY | 0.303 | 0.297 | Boltz-2 |
| tyk2_ejm_46 | 4GIH | 0.223 | 0.181 | Boltz-2 |
| tyk2_jmc_23 | 4GJ2 | 0.077 | 0.296 | FlashBind |
| jnk1_17124-1 | 2GMX | 0.609 | 0.502 | Boltz-2 |
| p38_p38a_3fln | 3FLN | 0.660 | 1.137 | FlashBind |
| p38_p38a_3flq | 3FLQ | 1.170 | 1.108 | Boltz-2 |
| p38_p38a_3flw | 3FLW | 0.704 | 1.111 | FlashBind |
| p38_p38a_3fly | 3FLY | 0.768 | 1.175 | FlashBind |
| p38_p38a_3flz | 3FLZ | 0.726 | 0.547 | Boltz-2 |
| p38_p38a_3fmh | 3FMH | 1.062 | 2.790 | FlashBind |
| p38_p38a_3fmk | 3FMK | 0.880 | 0.941 | FlashBind |

## Input audit and limitations

All 16 released SDF records matched the corresponding crystallographic ligand
in heavy-atom count and admitted a complete bond-order graph mapping. The
FlashBind documentation describes these SDFs as FABind+ docking outputs, but
the released archive does not include the original docking commands, random
seeds, or per-pose generation logs. This analysis therefore evaluates the
archived poses exactly as supplied; it does not independently reproduce or
audit FABind+ pose generation.

The released p38 receptor PDB has a missing chain identifier after residue
168. This causes structure parsers to split one continuously numbered kinase
construct into named and blank-chain fragments. I joined those two fragments
by deposited residue number in memory, without changing any coordinate, and
recorded the repair for every p38 comparison.

These 16 public complexes are a selected retrospective subset and may overlap
training data. The FlashBind and Boltz-2 inputs also differ materially:
FlashBind uses fixed experimental receptor coordinates, while Boltz-2 predicts
the entire complex. The result is useful for checking the supplied poses, but
it is not a prospective or receptor-input-matched docking benchmark.

## Reproduction

The protocol is frozen in `configs/benchmarks/fepplus4_crystal16_pose.yaml`.
Raw atom mappings and metrics are under
`runs/exp010_flashbind_crystal_pose/released_poses/`. Compact paired results
are in `paired_pose_rmsd.csv` and `paired_pose_comparison.json` beside this
file.

```bash
env PYTHONPATH=src \
  /home/yusheng/programs/flashbind-f161268/bin/python \
  scripts/analyze_flashbind_crystal_pose.py

/home/yusheng/programs/flashbind-f161268/bin/python \
  scripts/compare_pose_models.py \
  --flashbind runs/exp010_flashbind_crystal_pose/released_poses/pose_metrics.csv \
  --boltz2 reports/exp008_boltz2_crystal_pose/pose_metrics.csv \
  --output-dir reports/exp010_flashbind_crystal_pose
```
