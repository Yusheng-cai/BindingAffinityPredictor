# Experiment 007: FEP+ four-kinase reproduction

## Current result

Both released checkpoints completed all 87 compounds at seed 42 under the
feasible `boltz2_msa1024` protocol. The compound-count-weighted mean of the
four within-target Pearson correlations is **0.8066** for Nesso-1 (95%
bootstrap interval **0.7308--0.8731**) and **0.6265** for Boltz-2 (interval
**0.4845--0.7443**). These are close to the paper references of 0.80 and 0.66,
respectively.

The paired Nesso-minus-Boltz Pearson difference is **+0.1801**. Resampling the
same ligands for both models gives a 95% interval of **+0.0580--+0.3189**.
This result is target-heterogeneous: Boltz is stronger on CDK2 and JNK1,
whereas Nesso is much stronger on p38.

| Target | N | Nesso Pearson r | Boltz Pearson r | Nesso centered MAE | Boltz centered MAE |
|---|---:|---:|---:|---:|---:|
| CDK2 | 16 | 0.7480 | 0.8308 | 0.7234 | 0.5403 |
| TYK2 | 16 | 0.9167 | 0.8074 | 0.7492 | 0.6631 |
| JNK1 | 21 | 0.6663 | 0.7581 | 0.4500 | 0.4418 |
| p38 | 34 | 0.8692 | 0.3639 | 0.4247 | 0.7814 |
| Weighted | 87 | **0.8066** | **0.6265** | **0.5454** | **0.6333** |

Canonical model and paired metrics are in
[boltz2_msa1024_vs_nesso1_metrics.json](boltz2_msa1024_vs_nesso1_metrics.json).
Raw model outputs and canonical prediction tables remain under ignored
`runs/exp007_fepplus4_boltz2_nesso1/` paths.

## Boltz-2 feasibility history

The unchanged 8,192-sequence-MSA protocol was not feasible on the local 10 GiB
RTX 3080. All 16 TYK2 structure batches raised CUDA out-of-memory warnings.
The subsequent affinity pass requires `pre_affinity_*.npz` intermediates from
the structure pass and therefore terminated without scores. This is an
engineering/hardware failure, not an affinity-accuracy result.

The failed-run provenance is retained locally at
`runs/exp007_fepplus4_boltz2_nesso1/boltz2/seed42/tyk2/run.json`. A proposed
separately named `boltz2_msa1024` variant explicitly enabled Boltz's built-in
MSA subsampling to 1,024 sequences while retaining the parsed 8,192-row cap.
After a one-ligand and then 16-ligand TYK2 gate, the identical variant was run
on CDK2, JNK1, and p38. It completed 87/87 structures and 87/87 affinity
predictions with no model failures. The four target-run records total 5,921.8
s (98.7 min) on one RTX 3080, excluding smoke tests and one-time MSA creation;
the largest observed total GPU allocation was 9,988 MiB. The original
unsampled failure remains preserved rather than being overwritten.

## Scope of the result

This is a released-model reproduction with reconstructed preprocessing, not a
bitwise reproduction of either paper. Boltz's 1,024-row stochastic trunk MSA
is also a hardware-feasibility variant, not necessarily the authors' exact
evaluation setting. The labels are heterogeneous biochemical potencies (IC50
for CDK2/JNK1/p38 and Ki for TYK2), not equilibrium binding free energies. The
reported correlations quantify within-series ranking on four public
congeneric kinase sets; they do not demonstrate universal or
out-of-distribution affinity prediction.
