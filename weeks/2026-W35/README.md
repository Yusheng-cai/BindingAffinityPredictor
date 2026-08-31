# Week 2026-W35: Nesso-1 and Boltz-2 on Runs N' Poses

This week I studied four related questions: model speed, structural
generalization on Runs N' Poses, the benchmark sets commonly used by the FEP
community, and the internal architecture of Nesso-1.

## 1. I measured Nesso-1 and Boltz-2 inference time

I timed the same 50-system Runs N' Poses workload on one NVIDIA RTX 3080. The
timing came from the complete recorded model runs rather than a separate toy
benchmark.

| Model | Completed | Total wall time | Seconds per attempted system |
| --- | ---: | ---: | ---: |
| Nesso-1 v1.0.0 | 50/50 | 326.1 s | 6.52 s |
| Boltz-2 v2.2.1, MSA-1024 | 49/50 | 1181.7 s | 23.63 s |

In this specific local workflow, Nesso-1 was **3.62 times faster**. The Boltz-2
number excludes the earlier remote ColabFold MSA search, while Nesso-1's local
ESM-2 feature calculation is included. This is therefore a descriptive
single-run comparison, not a complete end-to-end latency benchmark or a set of
repeated timing trials.

## 2. I compared both models on Runs N' Poses

I selected 50 diverse protein-ligand complexes released after 1 June 2023,
with ten systems in each Boltz-2 familiarity interval from 0--20 through
80--100. Both models received the same deposited protein sequence and ligand
SMILES; Boltz-2 additionally received a sequence-derived MSA. Neither model
received the experimental structure or pocket during inference.

Nesso-1 completed all 50 systems. Boltz-2 completed 49; the 818-residue 8UK6
system exceeded the 10 GiB GPU limit. On the 49 systems completed by both
models:

- Boltz-2 recovered the experimental pocket more reliably: median 6 Å residue-pocket
  F1 was **0.889**, compared with **0.514** for Nesso-1, and Boltz-2 was better
  on 44/49 systems.
- Median 6 Å token-pair contact F1 was **0.780** for Boltz-2 and **0.462** for Nesso-1.
- Interface-distance MAE was mixed. Although the cohort medians were 0.881 Å
  for Boltz-2 and 1.387 Å for Nesso-1, Nesso-1 was better on 27 systems and
  Boltz-2 on 22, so there was no consistent paired advantage.
- Greater pre-cutoff structural familiarity was strongly associated with
  better accuracy for both models. For 6 Å residue-pocket F1, Spearman rho was **+0.677**
  for Nesso-1 and **+0.690** for Boltz-2. Each model used its own PDB
  structural-training cutoff: September 2021 for Nesso-1 and June 2023 for
  Boltz-2.

Only Boltz-2 generates explicit coordinates. Its median protein C-alpha RMSD
was 0.894 Å, and its median ligand heavy-atom RMSD after global protein
alignment was 3.355 Å. These are structural results; the selected Runs N'
Poses systems do not provide one uniform experimental affinity endpoint.

See the complete [paired result](../../reports/exp015_boltz2_nesso1_rnp_postcutoff50/README.md).

## 3. I reviewed benchmarks used by the FEP community

I compared the purpose, size, data vintage, and public availability of the
original FEP+/JACS benchmark, the four-series FEP+4 subset used in Week 34, D3R
Grand Challenges, the Hahn/OpenFF benchmark, the larger Ross collection, the
OpenFE public study, and the CASP16 affinity task.

The main conclusion is that most FEP benchmarks test a focused lead-optimization
question: given a prepared binding site, a plausible ligand pose, and a series
of related compounds, can a method recover their **relative** affinity
differences? They generally do not test blind pocket discovery or de novo pose
generation. The Week 34 FEP+4 benchmark is therefore four within-target kinase
ranking tests, not 87 independent tests of generalization to new proteins.
Larger Hahn/OpenFF and Ross collections provide useful follow-up diversity, but
their prepared structures and congeneric ligand-series design remain part of
the interpretation.

## 4. I studied the Nesso-1 architecture

I followed the released Nesso-1 forward pass from its input to its affinity
outputs:

```text
protein sequence + ligand SMILES
    -> residue and ligand-atom tokens
    -> ESM-2 protein context
    -> token-pair representation
    -> distance distributions
    -> automatically selected pocket
    -> continuous affinity and binder-probability outputs
```

The key architectural points are:

- ESM-2 provides a 1,280-dimensional contextual vector for each protein
  residue. Nesso learns a projection to its 384-dimensional token features.
- The central object is a pair tensor `z` with shape `N x N x 128`, describing
  every ordered token pair. It is not simply a dot-product similarity matrix.
- Forty-eight Pairformer blocks update `z` using triangle multiplication,
  triangle attention, and pair-transition layers.
- A 64-bin distogram predicts probability distributions over token-token
  distances. These predicted distances are used to identify a protein region
  within 22 Å, subject to a 256-token crop, for recycling and refinement.
- Two independent eight-block affinity modules operate on the predicted
  interface. Their two continuous estimates are averaged into the reported
  `affinity_pred_value`; binder probability remains a separate output.

Unlike Boltz-2, Nesso-1 stops at a probabilistic residue-center/ligand-heavy-atom
distance map rather than generating a complete full-atom Cartesian complex.
The executed architecture notebook combines a transparent miniature
calculation with a real saved Nesso-1 checkpoint trace.

## Files for this week

- [`report/main.tex`](report/main.tex): editable LaTeX source
- [`report/report.pdf`](report/report.pdf): compiled draft report
- [`report/figures/nesso_boltz_familiarity.pdf`](report/figures/nesso_boltz_familiarity.pdf): median bar charts with bootstrap intervals using each model's historical cutoff
- [`report/figures/nesso_boltz_runtime.pdf`](report/figures/nesso_boltz_runtime.pdf): observed Nesso-1 versus Boltz-2 MSA-1024 seconds per attempted system on the same RTX 3080
- [`report/figures/nesso_boltz_runtime.csv`](report/figures/nesso_boltz_runtime.csv): compact timing values from the paired `exp015` run
- [`code/README.md`](code/README.md): complete check/reproduction workflow for
  fresh-machine setup, the paired Runs N' Poses calculation, its recorded
  timing, and the two Nesso notebooks
- [`code/00_setup.sh`](code/00_setup.sh): pinned, checksum-verified acquisition
  of the software, model assets, and minimal Runs N' Poses data required here
- [`manifest.yaml`](manifest.yaml): draft scope and provenance
- [`../../data/manifests/rnp_boltz2_nesso1_postcutoff50.json`](../../data/manifests/rnp_boltz2_nesso1_postcutoff50.json): frozen paired analysis set
- [`notebooks/nesso_architecture_visualized.ipynb`](notebooks/nesso_architecture_visualized.ipynb): executed, interactive companion that follows token features through the pair tensor, Pairformer-style updates, distograms, cropping, and affinity outputs, ending with a real saved Nesso-1 tensor trace
- [`notebooks/boltz2_nesso1_rnp_postcutoff50.ipynb`](notebooks/boltz2_nesso1_rnp_postcutoff50.ipynb): executed paired comparison of 50 post-cutoff Runs N' Poses systems

Only those two notebooks belong to Week 35. Earlier Nesso-only familiarity
explorations and the RFdiffusion architecture notebook were moved intact to
`weeks/2026-W36/notebooks/` for possible follow-up work.
