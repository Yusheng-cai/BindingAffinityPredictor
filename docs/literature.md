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

## The exp001 biological and experimental reference

5. **Liang et al. (2013), “Lead Identification of Novel and Selective TYK2
   Inhibitors.”** Primary source for the ligand-46 (K_i) used here (Table 4,
   entry 46).
   [doi:10.1016/j.ejmech.2013.03.070](https://doi.org/10.1016/j.ejmech.2013.03.070)
6. **Liang et al. (2013), “Lead Optimization of a 4-Aminopyridine Benzamide
   Scaffold To Identify Potent, Selective, and Orally Bioavailable TYK2
   Inhibitors.”** Companion medicinal-chemistry series and assay context.
   [doi:10.1021/jm400266t](https://doi.org/10.1021/jm400266t)
7. **PDB 4GIH.** The 2.00 Å TYK2 JH1–inhibitor 46 coordinate reference used only
   after prediction. [RCSB](https://www.rcsb.org/structure/4GIH) and
   [doi:10.2210/pdb4GIH/pdb](https://doi.org/10.2210/pdb4GIH/pdb)
8. **OpenFF Protein-Ligand Benchmark.** Curated target, ligand, and measurement
   metadata from which the one-case manifest is frozen.
   [Repository](https://github.com/openforcefield/protein-ligand-benchmark)

## Benchmark design and statistics

9. **Hahn et al. (2022), “Best Practices for Constructing, Preparing, and
   Evaluating Protein-Ligand Binding Affinity Benchmarks.”** The main guide for
   assay consistency, chemical curation, uncertainty, metric choice, and
   reproducible reporting.
   [doi:10.33011/livecoms.4.1.1497](https://doi.org/10.33011/livecoms.4.1.1497)
10. **Schindler et al. (2020), “Large-Scale Assessment of Binding Free Energy
    Calculations in Active Drug Discovery Projects.”** Origin and context for a
    widely reused multi-target FEP benchmark that includes TYK2.
    [doi:10.1021/acs.jcim.0c00900](https://doi.org/10.1021/acs.jcim.0c00900)
11. **Ross et al. (2023), “The maximal and current accuracy of rigorous
    protein-ligand binding free energy calculations.”** Important context on
    experimental reproducibility and the ceiling it places on apparent method
    accuracy.
    [doi:10.1038/s42004-023-01019-9](https://doi.org/10.1038/s42004-023-01019-9)
12. **Buttenschoen et al. (2024), “PoseBusters: AI-based docking methods fail to
    generate physically valid poses or generalise to novel sequences.”** Shows
    why ligand RMSD alone is insufficient and motivates chemical/physical pose
    checks.
    [doi:10.1039/D3SC04185A](https://doi.org/10.1039/D3SC04185A)

## Leakage and generalization cautions

13. **Li et al. (2023), “Leak Proof PDBBind.”** Introduces splits designed to
    reduce protein- and ligand-similarity leakage; useful when we graduate from
    smoke testing to claims about generalization.
    [arXiv:2308.09639](https://arxiv.org/abs/2308.09639)
14. **Joeres et al. (2025), “Data splitting to avoid information leakage with
    DataSAIL.”** Formalizes multi-entity splits and demonstrates the harder
    regime where both protein and ligand are dissimilar to training examples.
    [Nature Communications](https://doi.org/10.1038/s41467-025-58606-8)
15. **Mattsson and Walters (2026), “Identifying and Addressing Systematic Data
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
