# Project Scope

## Objective

Understand the behavior, accuracy, limitations, and computational cost of protein-ligand affinity models through small reproducible experiments and progressively stronger benchmarks.

## Initial scope

- Start with off-the-shelf Boltz-2 inference.
- Test structure prediction, stochasticity, affinity ranking, binder classification, simple bias controls, and resource use.
- Prefer homogeneous single-assay measurements for quantitative ranking.
- Begin with small panels before scaling to high-throughput screening benchmarks.

## Non-goals for the first milestone

- Training or fine-tuning a model.
- Treating a learned affinity scalar as a rigorous binding free energy.
- Claiming broad generalization from a single target or a familiar retrospective benchmark.
- Running very large virtual screens before local runtime and memory behavior are measured.

## Planned model expansion

Adapters may later be added for Nesso-1, FlashBind/FlashAffinity, physicochemical baselines, ligand-only baselines, docking scores, and physics-based methods. Each addition should preserve the canonical input and result schemas while retaining the native raw outputs.

