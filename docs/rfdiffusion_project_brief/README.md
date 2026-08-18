# Physics-guided protein-binder design briefing

Presentation-ready overview of the proposed RFdiffusion → ProteinMPNN →
physics-guided scoring project. The same source supports two views:

- slide mode for a live talk;
- scroll mode for reading as a web page.

The diagrams are original schematics. Published results, proposed work, and
results already produced in this repository are labeled separately.

## Open locally

From the repository root:

```bash
python3 -m http.server 8000
```

Then open:

- slides: <http://localhost:8000/docs/rfdiffusion_project_brief/>
- reading view: <http://localhost:8000/docs/rfdiffusion_project_brief/?view=scroll>

In slide mode, use the arrow keys to navigate, `O` for the overview, `S` for
speaker notes, and `F` for fullscreen. To export from Chromium, append
`?print-pdf` and print with landscape orientation, no margins, and background
graphics enabled.

## Scope

This briefing answers four questions:

1. What do RFdiffusion and ProteinMPNN receive and produce?
2. How do their architectures represent three-dimensional protein geometry?
3. Where can a hydration-aware physical model enter the design loop?
4. What starter system gives a useful, falsifiable pilot experiment?

The proposed TEM-1/SHV-1 study is a project design, not a completed experiment.
The insulin-receptor structures shown as local evidence are engineering outputs
from experiments `exp003` and `exp004`; they have not been validated as binders.
