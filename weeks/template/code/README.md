# Weekly report workflow

Use this directory as the single entry point for the week's complete workflow:
setup, pinned downloads, input preparation, inference, analysis, and report
building. Prefer numbered Bash orchestration scripts and a
`workflow_manifest.yaml` that records source revisions, checkpoints, seeds,
inputs, outputs, and canonical Python entry points.

Every week begins with `00_setup.sh`. Its default mode must be read-only. It
must provide a fresh-machine execution mode that installs pinned software under
`/home/yusheng/programs`, clones source under `/home/yusheng/source`, downloads
ignored data/model assets from pinned sources only after explicit permission,
and validates checksums before later stages are allowed to run. If data have
separate terms of use, acceptance must also be explicit.

If the week has an interactive notebook, record its canonical path in the
workflow manifest and provide a check/launch entry point here. The notebook
should read canonical saved outputs and call reusable tested code; it must not
be the only implementation of a scientific metric.

Reference reusable analysis code under `src/` or `scripts/` at the Git revision
in the weekly manifest. Do not copy model repositories, environments, model
weights, raw predictions, or reusable scientific implementations into a weekly
folder. Potentially expensive or external operations must require an explicit
execution flag.
