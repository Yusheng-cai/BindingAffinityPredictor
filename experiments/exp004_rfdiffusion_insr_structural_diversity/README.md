# exp004: RFdiffusion structural diversity across ten seeds

## Status

**Completed on 2026-08-11.** Seed 42 was reused unchanged from `exp003`; seeds
43–51 were generated successfully, giving ten total designs with seeds 42–51.

## Question

When target coordinates, contig, hotspots, checkpoint and denoiser settings are
identical, how much do independently seeded RFdiffusion samples differ in their
length, intrinsic backbone shape, and pose against the target?

This is an exploratory sampling experiment. Ten structures cannot establish
convergence or a probability distribution over all possible binders.

## What is held fixed and what varies?

Held fixed:

- insulin-receptor target residues A1–A150;
- hotspot labels A59, A83 and A91;
- requested binder length range 70–100 residues;
- complex checkpoint, 50 diffusion steps and zero added denoiser noise;
- software environment and RTX 3080 hardware.

Varied:

- deterministic design number/random seed, 42 through 51;
- the length sampled within 70–100 residues;
- initial binder noise and therefore the generated backbone and interface pose.

## Why more than one diversity metric?

A binder can have approximately the same fold but occupy a different location or
orientation on the target. Conversely, two binders can contact similar target
residues while having different folds. The analysis therefore keeps two ideas
separate:

1. **Target-relative pose diversity.** Align all target Cα atoms, then compare
   binder centers, target-contact sets and symmetric nearest-neighbor distances.
2. **Intrinsic shape diversity.** Resample each ordered binder Cα trace to 64
   arc-length positions and optimally superpose the resampled curves before
   calculating RMSD.

The resampled-trace RMSD is a transparent descriptive metric for unequal-length
chains. It is not TM-score and should not be interpreted using TM-score
thresholds. Cα contact counts are likewise coarse geometric descriptors rather
than atomistic interface energies.

No new ProteinMPNN or independent structure-prediction runs are part of this
experiment. That avoids conflating RFdiffusion backbone diversity with later
sequence sampling and validation stages.

## Result in brief

The ensemble is structurally diverse under both prespecified views:

- sampled lengths span 71–98 residues;
- intrinsic resampled-trace RMSD has median 13.36 Å and range 4.68–16.45 Å;
- target-aligned pose Chamfer distance has median 5.08 Å and range 2.66–10.40 Å;
- binder center displacement has median 6.37 Å and reaches 14.20 Å;
- contact-set Jaccard similarity has median 0.50 and range 0.36–0.78.

Final RFdiffusion binder confidence occupies a very narrow 0.9909–0.9929 range.
Thus, this internal output does not discriminate the strong geometric variation
in this small panel. It is not a binding confidence.

All ten backbones contact the hotspot neighborhood at the coarse Cα level. No
design contains an interchain backbone-atom pair below 2 Å, although side chains
have not yet been designed for seeds 43–51. These checks do not establish an
energetically favorable or experimentally realizable interface.

Exact per-design and all 45 pairwise values are in the tracked report tables.
The website displays the target-aligned ensemble and an interactive heat map.
