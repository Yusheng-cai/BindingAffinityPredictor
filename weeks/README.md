# Weekly research reports

Each ISO week has a self-contained report snapshot plus a complete, thin
workflow entry point. The weekly code may orchestrate source acquisition,
environment setup, downloads, inference, analysis, and report construction.
Reusable scientific implementations remain in src/, scripts/, and experiments/
rather than being copied into weekly folders.

## Layout

    weeks/
    ├── template/
    └── YYYY-Www/
        ├── manifest.yaml
        ├── report/
        │   ├── main.tex
        │   ├── figures/
        │   └── report.pdf
        ├── notebooks/
        ├── results/
        └── code/
            ├── README.md
            ├── workflow_manifest.yaml
            ├── config.sh
            ├── NN_workflow_stage.sh
            ├── build_assets.py
            └── build_report.sh

The manifest connects prose to an exact project Git revision, experiment
definitions, commands, citations, figure sources, and the compiled PDF hash.
Raw model outputs remain under ignored runs/ directories.

## GitHub weekly publishing policy

The GitHub repository should contain everything needed to understand and
reconstruct a calculation, but not multi-gigabyte external artifacts. Commit:

- the weekly report source and PDF;
- the weekly workflow entry points and pinned workflow manifest;
- canonical experiment YAML, model and benchmark configs, and input manifests;
- reusable Python implementations and focused tests;
- analysis notebooks that read canonical saved outputs;
- compact reviewed metrics, tables, and figures.

Do not commit model repositories, Python environments, checkpoints, caches,
MSAs, downloaded raw data, or raw prediction trees. The workflow scripts must
download those items from pinned sources, validate checksums when available,
and record their local run provenance.

For each new week:

1. copy `weeks/template/` to the new ISO-week directory;
2. define the week's questions and link its experiment configurations;
3. add or update safe numbered workflow entry points;
4. run focused tests and the report build;
5. record the project Git revision, run records, PDF hash, and external source
   revisions in `manifest.yaml` and `code/workflow_manifest.yaml`;
6. commit the reviewed weekly snapshot before starting the next week's work.

## Overleaf policy

Each Overleaf-connected week records its project, remote revision, and
dedicated local clone in that week's manifest. The scientific repository is
not pushed wholesale to Overleaf.

For every edit:

1. pull the latest Overleaf source before editing;
2. treat the latest Overleaf text as the writing source of truth;
3. preserve report-specific code and scientific protocols in this repository;
4. compile locally before pushing report changes;
5. mirror the reviewed source and compact figures into the weekly snapshot;
6. record both the Overleaf revision and project Git revision in manifest.yaml.

The Overleaf project contains only the editable report, compact figures, and
bibliography. Frozen PDFs and experiment provenance remain available here.
