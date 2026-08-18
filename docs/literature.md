# Literature map

This is an annotated starting set, grouped by the decision each source informs.
Model claims are kept separate from independent benchmark guidance and target
data.

## Model and architectural background

1. **Passaro et al. (2025), “Boltz-2: Towards Accurate and Efficient Binding
   Affinity Prediction.”** The primary model report: architecture, training
   objectives, affinity heads, and the authors' benchmarks.
   [doi:10.1101/2025.06.14.659707](https://doi.org/10.1101/2025.06.14.659707)
2. **Wohlwend et al. (2024), “Boltz-1: Democratizing Biomolecular Interaction
   Modeling.”** Open predecessor and useful background for the structure trunk.
   [doi:10.1101/2024.11.19.624167](https://doi.org/10.1101/2024.11.19.624167)
3. **Abramson et al. (2024), “Accurate structure prediction of biomolecular
   interactions with AlphaFold 3.”** Architectural context for full-atom
   diffusion-based biomolecular cofolding.
   [doi:10.1038/s41586-024-07487-w](https://doi.org/10.1038/s41586-024-07487-w)
4. **Official Boltz prediction documentation.** Operational source of truth for
   YAML inputs, MSA behavior, inference flags, output fields, and affinity-score
   semantics. [GitHub](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)

## Protein surfaces, site prediction, and flexible structural modeling

5. **Gainza et al. (2020), “Deciphering interaction fingerprints from protein
   molecular surfaces using geometric deep learning.”** Introduces MaSIF and
   its three native tasks: seven-class cofactor prediction for already defined
   binding pockets,
   protein–protein interface-site prediction, and surface-fingerprint matching.
   It is a surface representation framework, not a small-molecule affinity or
   free-energy predictor.
   [doi:10.1038/s41592-019-0666-6](https://doi.org/10.1038/s41592-019-0666-6)
6. **Sverrisson et al. (2021), “Fast End-to-End Learning on Protein Surfaces.”**
   Introduces dMaSIF, which generates a differentiable surface point cloud from
   atomic coordinates and chemical types on the fly, avoiding MaSIF's
   precomputed mesh-feature bottleneck.
   [CVPR paper](https://openaccess.thecvf.com/content/CVPR2021/html/Sverrisson_Fast_End-to-End_Learning_on_Protein_Surfaces_CVPR_2021_paper.html)
7. **Gazizov et al. (2026), “AF2BIND: predicting small-molecule binding sites
   using the pair representation of AlphaFold2.”** A protein-only site model
   that predicts ligand-contact residues without the true ligand, an MSA, or a
   template. It is a useful modern comparator for the site-discovery layer.
   [doi:10.1038/s41592-026-03011-2](https://doi.org/10.1038/s41592-026-03011-2)
8. **Lu et al. (2024), “DynamicBind: predicting ligand-specific protein-ligand
   complex structure with a deep equivariant generative model.”** Flexible
   docking from an apo-like receptor structure and a ligand, with iterative
   ligand and protein-coordinate updates and a separate PDBbind-trained affinity
   module.
   [doi:10.1038/s41467-024-45461-2](https://doi.org/10.1038/s41467-024-45461-2)

## The exp001 biological and experimental reference

9. **Liang et al. (2013), “Lead Identification of Novel and Selective TYK2
   Inhibitors.”** Primary source for the ligand-46 (K_i) used here (Table 4,
   entry 46).
   [doi:10.1016/j.ejmech.2013.03.070](https://doi.org/10.1016/j.ejmech.2013.03.070)
10. **Liang et al. (2013), “Lead Optimization of a 4-Aminopyridine Benzamide
   Scaffold To Identify Potent, Selective, and Orally Bioavailable TYK2
   Inhibitors.”** Companion medicinal-chemistry series and assay context.
   [doi:10.1021/jm400266t](https://doi.org/10.1021/jm400266t)
11. **PDB 4GIH.** The 2.00 Å TYK2 JH1–inhibitor 46 coordinate reference used only
   after prediction. [RCSB](https://www.rcsb.org/structure/4GIH) and
   [doi:10.2210/pdb4GIH/pdb](https://doi.org/10.2210/pdb4GIH/pdb)
12. **OpenFF Protein-Ligand Benchmark.** Curated target, ligand, and measurement
   metadata from which the one-case manifest is frozen.
   [Repository](https://github.com/openforcefield/protein-ligand-benchmark)

## The exp002 galectin-3–galactose reference

13. **PDB 9D63.** Wild-type human galectin-3 carbohydrate-recognition domain
   with beta-D-galactose; a 1.15 A X-ray structure deposited in 2024 and first
   released on 2026-03-18. The coordinates are used only after prediction.
   [RCSB](https://www.rcsb.org/structure/9D63) and
   [doi:10.2210/pdb9D63/pdb](https://doi.org/10.2210/pdb9D63/pdb)
14. **Denavit et al. (2023), “Selectively Modified Lactose and
    N-Acetyllactosamine Analogs at Three Key Positions to Afford Effective
    Galectin-3 Ligands.”** Experimental and structural context for galectin-3
    glycan recognition; reports an ITC (K_d) of 91 uM for methyl
    beta-D-lactoside and stronger modified ligands. This is context, not a
    quantitative label for free galactose in exp002.
    [PubMed](https://pubmed.ncbi.nlm.nih.gov/36835132/)
15. **Diehl et al. (2009), “Conformational entropy changes upon lactose binding
    to the carbohydrate recognition domain of galectin-3.”** NMR and molecular
    simulation evidence that binding thermodynamics include conformational
    entropy even when gross structural changes are small; useful caution against
    treating pose recovery as affinity recovery.
    [doi:10.1007/s10858-009-9356-5](https://doi.org/10.1007/s10858-009-9356-5)

## Benchmark design and statistics

16. **Hahn et al. (2022), “Best Practices for Constructing, Preparing, and
   Evaluating Protein-Ligand Binding Affinity Benchmarks.”** The main guide for
   assay consistency, chemical curation, uncertainty, metric choice, and
   reproducible reporting.
   [doi:10.33011/livecoms.4.1.1497](https://doi.org/10.33011/livecoms.4.1.1497)
17. **Schindler et al. (2020), “Large-Scale Assessment of Binding Free Energy
    Calculations in Active Drug Discovery Projects.”** Origin and context for a
    widely reused multi-target FEP benchmark that includes TYK2.
    [doi:10.1021/acs.jcim.0c00900](https://doi.org/10.1021/acs.jcim.0c00900)
18. **Ross et al. (2023), “The maximal and current accuracy of rigorous
    protein-ligand binding free energy calculations.”** Important context on
    experimental reproducibility and the ceiling it places on apparent method
    accuracy.
    [doi:10.1038/s42004-023-01019-9](https://doi.org/10.1038/s42004-023-01019-9)
19. **Buttenschoen et al. (2024), “PoseBusters: AI-based docking methods fail to
    generate physically valid poses or generalise to novel sequences.”** Shows
    why ligand RMSD alone is insufficient and motivates chemical/physical pose
    checks.
    [doi:10.1039/D3SC04185A](https://doi.org/10.1039/D3SC04185A)

## Leakage and generalization cautions

20. **Li et al. (2023), “Leak Proof PDBBind.”** Introduces splits designed to
    reduce protein- and ligand-similarity leakage; useful when we graduate from
    smoke testing to claims about generalization.
    [arXiv:2308.09639](https://arxiv.org/abs/2308.09639)
21. **Joeres et al. (2025), “Data splitting to avoid information leakage with
    DataSAIL.”** Formalizes multi-entity splits and demonstrates the harder
    regime where both protein and ligand are dissimilar to training examples.
    [Nature Communications](https://doi.org/10.1038/s41467-025-58606-8)
22. **Mattsson and Walters (2026), “Identifying and Addressing Systematic Data
    Leakage in Protein-Ligand Affinity Benchmarks.”** A recent preprint directly
    examining leakage in cofolding affinity benchmarks, including FEP+/OpenFE;
    its conclusions should be treated as pre-peer-review evidence.
    [doi:10.64898/2026.06.29.735309](https://doi.org/10.64898/2026.06.29.735309)

## Reading policy for this project

The Boltz-2 report tells us what the authors intended and measured; it is not an
independent validation. The benchmark-practice, pose-validity, experimental-
reproducibility, and leakage papers determine how we design our own tests and
limit our conclusions. A one-case smoke test is the beginning of that ladder,
not the end.
