# RFdiffusion hotspot-conditioning ablation

## Result in one sentence

With all configured variables except the hotspot tensor held fixed, the
hotspot-guided binder backbone came closer to at least one of A59, A83 or A91
than its no-hotspot counterpart in 9 of 10 matched seeds.

## Design

Both arms used the insulin-receptor example target (A1–A150), the
`Complex_base_ckpt.pt` checkpoint, contig `[A1-150/0 70-100]`, design
indices/seeds 42–51, 50 reverse-diffusion steps, deterministic execution, and
zero additional Cα/frame noise. The guided arm supplied
`ppi.hotspot_res=[A59,A83,A91]`; the unguided arm omitted that feature. All ten
paired seeds were retained.

This comparison required an explicit checkpoint override in the unguided arm.
Without it, RFdiffusion automatically selects the monomer `Base_ckpt.pt` when
hotspots are absent. Using that different checkpoint would confound hotspot
conditioning with model weights. The first command stopped before producing a
design; the successful run explicitly pinned the complex checkpoint.

## Prespecified geometric results

| Quantity | Guided | Unguided |
| --- | ---: | ---: |
| Minimum Cα distance to any nominal hotspot, median | 6.38 Å | 9.72 Å |
| Minimum Cα distance to any nominal hotspot, mean | 6.24 Å | 10.92 Å |
| Binder center to hotspot centroid, median | 14.93 Å | 17.41 Å |
| Designs with an interchain backbone pair below 2 Å | 0/10 | 0/10 |

Matched-arm comparisons had a median contact-set Jaccard similarity of 0.245,
a median target-aligned binder-pose symmetric Chamfer distance of 9.46 Å, and a
median binder-center displacement of 16.32 Å. Thus, removing the hotspot
feature can materially change both contacted target residues and the overall
generated pose. Individual values and SHA-256 hashes for every PDB and `.trb`
file are recorded in `paired_designs.csv`.

The center-to-centroid metric was not directionally consistent across matched
pairs, whereas the local minimum-hotspot distance favored the guided arm in
9/10 pairs. This distinction matters: hotspot conditioning asks the model to
place an interface near nominated residues; it does not require the geometric
center of a variable-length binder to sit on the hotspot centroid.

## Scientific interpretation

This result supports a mechanistic statement about conditioning: the hotspot
feature steers generated backbone placement toward the nominated target
region. It does **not** show that either arm binds, that the guided structures
have higher affinity, or that the unguided complex model performs validated
pocket discovery. RFdiffusion itself warns that this complex checkpoint was
trained with hotspot inputs and is being evaluated without them in the
unguided arm.

The structures remain poly-glycine backbones. ProteinMPNN sequence design,
all-atom reconstruction, independent complex prediction, physical interface
checks, and experimental validation are all downstream requirements.

## Experimental site-1 reference

We additionally aligned the experimental 4OGA insulin-receptor site-1 complex
to the RFdiffusion target. Chain E residues 6–155 of 4OGA correspond directly
to target A1–A150 and align at 0.423 Å Cα RMSD. This confirms that the tutorial
target is a close Rosetta model of the experimentally observed L1 geometry.

After applying the same transform to insulin and the receptor αCT peptide, the
guided ensemble covered a median 58% of modeled insulin Cα positions and 80%
of modeled αCT Cα positions within 8 Å. The unguided medians were 45% and 40%,
respectively. Several minimum generated-binder-to-experimental-component Cα
distances are below 2 Å. Those values signify steric overlap if the structures
were present simultaneously—not favorable atomic contacts. The generated
poses are therefore most naturally viewed as occupying or mimicking portions
of the composite insulin/αCT site, potentially competitively.

This remains a coarse Cα comparison. It does not establish binding or whether
a binder would inhibit, activate, or have no functional effect on the intact
glycosylated receptor. Per-design results are recorded in
`experimental_reference_4oga.csv`; the downloaded coordinates are identified
by the tracked manifest `data/manifests/4oga_insulin_receptor_site1.json`.

## Reproduction

The exact generation command, checkpoint checksum, source revision, hardware,
runtime, memory, failure record, and stopping rule are in
`experiments/exp005_rfdiffusion_insr_hotspot_ablation/experiment.yaml`. Recreate
the compact tables from the Git-ignored raw outputs with:

```bash
python3 scripts/analyze_rfdiffusion_hotspot_ablation.py \
  --guided-seed42-root runs/exp003_rfdiffusion_insr_binder_smoke/rfdiffusion/seed42/raw \
  --guided-ensemble-root runs/exp004_rfdiffusion_insr_structural_diversity/rfdiffusion/raw \
  --unguided-root runs/exp005_rfdiffusion_insr_hotspot_ablation/rfdiffusion/no_hotspots/raw \
  --report-dir reports/exp005_rfdiffusion_insr_hotspot_ablation
```

Recreate the experimental-reference comparison with:

```bash
python3 scripts/analyze_insr_experimental_reference.py \
  --experimental-pdb data/raw/experimental_references/4OGA.pdb \
  --guided-seed42-root runs/exp003_rfdiffusion_insr_binder_smoke/rfdiffusion/seed42/raw \
  --guided-ensemble-root runs/exp004_rfdiffusion_insr_structural_diversity/rfdiffusion/raw \
  --unguided-root runs/exp005_rfdiffusion_insr_hotspot_ablation/rfdiffusion/no_hotspots/raw \
  --report-dir reports/exp005_rfdiffusion_insr_hotspot_ablation
```
