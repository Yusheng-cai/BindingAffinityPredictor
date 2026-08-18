# exp002 result: Boltz-2 galectin-3–galactose

## Outcome

Pinned Boltz-2 completed both structure and affinity inference for the human
galectin-3 CRD with beta-D-galactose on one RTX 3080. After the prediction was
complete, comparison with the withheld 1.15 A PDB 9D63 structure showed a
0.475 A ligand heavy-atom RMSD after global protein alignment. This is a strong
single-case pose recovery, not evidence of affinity accuracy or generalization.

## Exact system and withheld information

- Protein: exact 139-residue wild-type human galectin-3 CRD construct from
  9D63 polymer entity 1.
- Ligand: beta-D-galactopyranose, using the stereochemical SMILES from RCSB CCD
  `GAL`.
- Inference inputs: sequence, SMILES, and server-generated MSA only.
- Excluded from inference: 9D63 coordinates, templates, pocket residues,
  contacts, constraints, experimental ligand pose, and affinity labels.
- Reference: [PDB 9D63](https://www.rcsb.org/structure/9D63), released on
  2026-03-18 and downloaded only after model completion.

The post-Boltz release of 9D63 lowers the risk of direct coordinate leakage.
It does not make this a strict protein- or ligand-novelty test: galectin-3,
related structures, and galactose chemistry were already public.

## What happened at each stage

| Stage | Recorded result | Meaning |
| --- | --- | --- |
| Input parsing | 2 chains, 139 protein residues, 1 ligand residue, 1,127 real atoms, 12 ligand bonds | Boltz constructed the intended protein plus one galactose; no template or constraint was present. |
| MSA request | One 139-residue protein submitted to the public ColabFold MMseqs2 endpoint | Neither ligand nor reference coordinates were sent. |
| MSA parsing | 6,772 raw CSV sequence rows; 6,387 processed sequences; 887,793 aligned residue symbols; 9,667 deletion records | Boltz removed or consolidated some raw rows during parsing. All 6,387 processed sequences were below the 8,192 cap. |
| Tokenization | 151 tokens = 139 residue tokens + 12 ligand-atom tokens | The ligand is atom-tokenized while each standard protein residue contributes one main token. |
| Atom features | 1,127 atoms padded to 1,152; 66 RDKit bound entries; 5 chiral-center constraints | Padding is computational; the extra 25 atom slots are masked. Ligand stereochemistry is explicitly represented. |
| Structure trunk | 3 recycles; final `s` shape `1 x 151 x 384`; final `z` shape `1 x 151 x 151 x 128` | `s` is the final single-token embedding; `z` is the final pair embedding. These raw values are not probabilities or energies. |
| Structure diffusion | 1 sample, seed 42, 200 steps, step scale 1.5 | The stock CLI writes the final coordinates, not all 200 transient coordinate states. |
| Confidence | Full pLDDT, PAE, PDE arrays and scalar confidence JSON | Confidence is model self-assessment, not experimental validation. |
| Affinity path | 5 recycles hard-coded by the pinned CLI; 5 diffusion samples, 200 steps | This is a separate checkpoint/pathway that reloads the top predicted structure. |
| Withheld scoring | Protein and ligand graph alignment against 9D63 after inference | This tests pose recovery without giving Boltz the holo coordinates. |

The model-ready unbatched feature dictionary, including every tensor name,
shape, and data type, is stored in the ignored local
`runs/.../seed_42/intermediate_audit.json`. Full arrays remain in the raw Boltz
directory.

## Structure and pose results

| Readout | Result |
| --- | ---: |
| Protein C-alpha RMSD, all 139 residues | 0.506 A |
| Protein pocket C-alpha RMSD, 7 reference-contact residues | 0.063 A |
| Ligand heavy-atom RMSD after global protein fit | 0.475 A |
| Ligand heavy-atom RMSD after pocket protein fit | 0.538 A |
| Ligand internal RMSD after ligand-only fit | 0.052 A |
| 4 A contact precision | 1.000 (6/6) |
| 4 A contact recall | 0.857 (6/7) |
| 4 A contact Jaccard index | 0.857 |

The unique chirality-aware graph mapping matched all 12 ligand heavy atoms.
Predicted contacts were His47, Asn49, Arg51, Asn63, Trp70, and Glu73 in input
sequence numbering; all six occur in the experimental contact set. The only
experimental contact not crossing the predicted 4 A threshold was Val61:
3.960 A experimentally versus 4.129 A in the prediction. This illustrates why
contact recovery depends on the cutoff.

All five coordinate-derived ligand stereocenters matched the reference CIP
assignments (`R,R,R,S,R` in the predicted-molecule atom order). Bond-length
differences relative to the crystal ligand had 0.022 A RMSE and 0.044 A maximum
absolute deviation. PoseBusters was not installed, so these are focused graph,
chirality, bond-length, and contact checks rather than the complete PoseBusters
suite.

## Confidence outputs

| Native quantity | Value |
| --- | ---: |
| confidence score | 0.9863 |
| pTM | 0.9890 |
| ipTM | 0.9867 |
| ligand ipTM | 0.9867 |
| complex pLDDT | 0.9862 |
| mean ligand-token pLDDT | 0.9823 |
| complex PDE | 0.2601 A |

The PAE array is directional. Its mean protein-to-ligand block was 1.01 A,
whereas the ligand-to-protein block was 9.79 A, so collapsing the full PAE
matrix to one symmetric number would hide model behavior. The independently
computed pose RMSD is the relevant external check here.

## Affinity outputs

| Native field | Ensemble | Component 1 | Component 2 |
| --- | ---: | ---: | ---: |
| `affinity_pred_value` | 2.1040 | 2.2871 | 1.9209 |
| `affinity_probability_binary` | 0.5745 | 0.6734 | 0.4757 |

On Boltz's documented approximate `log10(IC50/uM)` interpretation, 2.104 would
correspond to roughly 127 uM. This is only a model output: 9D63 supplies no
matched quantitative assay, and free galactose must not borrow a lactose or
lactoside affinity from another paper. The component binary heads straddle 0.5,
which also makes the ensemble classification borderline despite the excellent
pose.

Boltz recorded `affinity.mw = 173.1` for this SMILES, while the RCSB/RDKit
formula `C6H12O6` gives 180.156 Da. At this revision the parser adds explicit
hydrogens, removes them with sanitization disabled, and then calls RDKit
`MolWt`, leaving an inconsistent hydrogen state. Molecular-weight correction
was disabled, so the erroneous stored value did not modify this prediction.
The correction option should not be benchmarked until this behavior is isolated
in a dedicated test or fixed upstream.

## Runtime and provenance

- Boltz package 2.2.1; source commit
  `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc`.
- Structure checkpoint SHA-256:
  `090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`.
- Affinity checkpoint SHA-256:
  `dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e`.
- Python 3.11.15; PyTorch 2.5.1+cu121; CUDA inference with bfloat16 mixed
  precision and `--no_kernels`.
- GPU: NVIDIA RTX 3080, 10,240 MiB; driver 535.230.02.
- Wall time: 2 min 21.21 s; logged structure inference: 14 s; logged affinity
  inference: 20 s; MSA server response: approximately 2.5 s.
- Maximum host resident set: 8,053,016 KiB; no swap; exit status 0.
- Peak GPU memory was not captured, a provenance gap to correct before the next
  run.

The repository revision at execution was
`122da1b0c3d87fb616aa1f56a4f19ae588c4ae71`; the working tree was dirty because
the exp002 protocol and prior documentation edits were not yet committed.
Checksums for the exact input and local artifacts are retained in the ignored
run record.

## Interpretation and next decision

This case demonstrates that the pipeline runs locally and that one sampled pose
recovers this galectin-3–galactose crystal geometry extremely closely. It also
demonstrates that high structural confidence and excellent pose recovery need
not yield a decisive binary affinity output.

The next scientifically clean step is to review this run, then choose between:

1. repeat two additional seeds to measure pose and affinity stochasticity; or
2. define a small galectin-3 ligand series with measurements from one assay,
   which is necessary to test affinity ranking.

No correlation, calibration, affinity-accuracy, selectivity, or broad
generalization claim is supported by this single unlabeled complex.
