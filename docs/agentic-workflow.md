# Agentic Research Workflow

The repository separates scientific intent from execution artifacts so that a human or agent can reproduce, audit, and extend an experiment safely.

## Experiment lifecycle

1. **Define**
   - Create `experiments/<experiment-id>/README.md` and `experiment.yaml`.
   - State the question, hypothesis, input manifest, model configurations, metrics, exclusions, and stopping criteria.
2. **Validate inputs**
   - Validate protein constructs, ligand stereochemistry, units, assay identity, qualifiers, duplicates, and missing values.
   - Freeze a small manifest under `data/manifests/`.
3. **Resolve provenance**
   - Pin source revisions, checkpoint identifiers, environment locks, model arguments, seeds, and external structure/MSA provenance.
4. **Execute**
   - Run through a model adapter.
   - Write raw artifacts beneath `runs/<experiment-id>/<model>/<run-id>/`.
5. **Verify**
   - Check completion, schema validity, failure counts, score direction, expected output files, and obvious chemical/structural problems.
6. **Analyze**
   - Normalize outputs without discarding native values.
   - Calculate predefined metrics and uncertainty intervals.
   - Compare against simple baselines and inspect individual failures.
7. **Report**
   - Promote only compact, reviewed results to `reports/<experiment-id>/`.
   - Record deviations from the frozen protocol and any durable decision.

## Fresh-machine bootstrap

Every weekly computational workflow starts with `code/00_setup.sh`. Its
read-only mode verifies prerequisites, exact source revisions, external model
assets, and benchmark-file checksums. Its write mode must require explicit
permission before network downloads or installation and must require a
separate explicit acknowledgement when a dataset has governing terms. Source
trees belong under `/home/yusheng/source`, installed environments under
`/home/yusheng/programs`, and large data/model artifacts in Git-ignored paths.

The bootstrap downloads only the minimal external subset needed for that
week. Later workflow stages may assume Stage 00 has passed, but may not quietly
download their own datasets or checkpoints. Remote sequence submission, such
as ColabFold MSA generation, remains a distinct explicitly authorized stage.

## Canonical run layout

```text
runs/<experiment-id>/<model>/<run-id>/
├── run.json              # provenance, resolved settings, hardware, status
├── predictions.*         # canonical prediction table
├── raw/                  # unmodified model outputs
├── structures/           # generated complexes when applicable
└── logs/                 # stdout, stderr, and resource measurements
```

The complete `runs/` tree is local and ignored by Git.

## Agent handoff rules

- Leave a concise status in the experiment README or report, not only in chat history.
- Distinguish planned, running, failed, and completed work.
- Do not reinterpret a failed run as a negative scientific result.
- Never change labels, filters, metrics, or sampling settings after viewing results without documenting a new experiment revision.
- Prefer small verified steps: one input, then a tiny panel, then the full defined experiment.
