# exp006: ProteinMPNN fixed-backbone sequence ensemble

## Status

**Generation completed and audited on 2026-08-16.**

This experiment expands the four-sequence engineering smoke test in exp003
into a prespecified ensemble of 100 ProteinMPNN samples. It does not generate a
new backbone and does not evaluate binding affinity.

## Scientific question

Given one fixed RFdiffusion target–binder backbone, how diverse are the
sequences ProteinMPNN considers compatible with it, and which positions are
constrained by the receptor interface?

## Fixed input

The sole structural input is the audited RFdiffusion seed-42 complex:

    runs/exp003_rfdiffusion_insr_binder_smoke/
    └── rfdiffusion/seed42/raw/design_ppi_42.pdb

Chain A is the 150-residue insulin-receptor target and remains fixed. Chain B is
the 90-residue RFdiffusion binder backbone and is the only designed chain. Both
chains remain visible to ProteinMPNN.

This backbone was selected because it was frozen in the earlier smoke-test
protocol, not because it had a favorable downstream score.

## Frozen generation protocol

- ProteinMPNN checkpoint: v_48_020
- designed chain: B
- fixed context chain: A
- sampling temperature: 0.1
- backbone noise: 0.0 Å
- batches: 10
- samples per batch: 10
- batch seeds: 6000 through 6009 inclusive
- native probability and score archives: retained
- execution device for this run: CPU, because an unrelated GROMACS simulation
  is actively using the GPU
- Python environment: /home/yusheng/anaconda3/bin/python
- PyTorch version: 2.5.1
- CPU safeguards: GPU hidden, one OpenMP/BLAS thread, low process and I/O
  priority

The earlier RFdiffusion environment contains PyTorch 1.9.1 with 3.5 GB of
CUDA-linked libraries on the busy HDD. A first CPU attempt spent 326.6 seconds
paging libraries without producing a sequence and was interrupted. A
disposable compatibility test showed that the existing NVMe-based PyTorch
2.5.1 environment loads the same pinned ProteinMPNN checkpoint and produces the
expected FASTA, probability, and score outputs in 3.7 seconds. The production
run therefore uses PyTorch 2.5.1 and records the complete environment. Its
stochastic samples are not expected to be bitwise comparable to the earlier
PyTorch 1.9.1 smoke-test samples.

The four exp003 samples are not part of this 100-sequence ensemble.

## Native score interpretation

ProteinMPNN's designed-chain score is a mean negative log-probability. Lower is
more probable under the sequence model for the supplied structural context. It
is not a binding free energy, affinity, folding probability, or probability of
experimental success.

## Success criteria

Engineering success requires:

1. all ten batch commands exit with status zero;
2. each batch contains exactly ten sampled sequences in addition to the native
   poly-glycine reference record;
3. every sampled binder sequence has length 90;
4. all native probability and score archives are retained;
5. the aggregate table contains exactly 100 samples;
6. commands, hashes, seeds, timing, environment, hardware snapshot, and failure
   state are recorded.

Sequence uniqueness is measured but is not a success criterion; duplicate
samples are legitimate observations from the model distribution.

## Stopping rule

Stop after the 100 samples are generated and audited. Do not select a “best”
sequence or begin all-atom prediction until the ensemble analysis and selection
policy are reviewed.

## Completed sequence analysis

Pairwise sequence identity, average-linkage clustering, and per-position
amino-acid frequencies and entropy have been calculated. Cluster counts are
reported over several dendrogram cuts because there is no unique biological
identity threshold. Entropy is mapped to the real chain-B Cα coordinates in a
rotatable HTML view.

Basic sequence liabilities and separate summaries for provisional interface
and non-interface positions remain pending. A receptor-omission control is
scientifically useful but is not part of this generation run.

## Observed generation result

All ten CPU batches completed successfully. The run generated 100 sampled
90-residue binder sequences, all of which were unique in this finite sample.
The aggregate designed-chain ProteinMPNN score ranged from 0.7159 to 0.9584
with a median of 0.8110. Pairwise sequence identity ranged from 37.8% to 90.0%
with a median of 66.7%.

Empirical per-position entropy across the sampled sequences ranged from 0 to
1.273 nats, with a median of 0.648 nats. Seventeen of 90 positions were
invariant in the 100 samples. These are finite-sample descriptions of the
ProteinMPNN distribution at temperature 0.1, not measurements of biological
conservation.

The ensemble is compositionally narrow: Ala, Glu, and Lys account for 72.1% of
all designed residues. This observation motivates explicit composition,
structure-recovery, and interface analyses before candidate selection. It does
not establish that the sequences fold or bind.

The production generation required 162.899 seconds using CPU-only PyTorch
2.5.1. Peak batch resident memory was 658,016 KiB. The GPU was hidden from
ProteinMPNN.
