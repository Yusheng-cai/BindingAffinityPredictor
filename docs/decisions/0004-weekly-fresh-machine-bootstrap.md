# Decision 0004: every weekly workflow starts from a fresh-machine bootstrap

## Status

Accepted on 2026-08-30.

## Decision

Every new computational week will include `weeks/YYYY-Www/code/00_setup.sh` as
its first entry point. The script will have a read-only verification mode and
an explicit execution mode that can acquire all otherwise missing software,
model assets, and benchmark data needed by that week's canonical workflow.

The bootstrap must:

- pin source repositories and model/data revisions;
- validate downloaded artifacts with recorded SHA-256 checksums when available;
- require explicit permission before downloading or installing anything;
- require a distinct acceptance flag for governed data terms;
- place source under `/home/yusheng/source`, installed environments under
  `/home/yusheng/programs`, and bulk artifacts only in Git-ignored locations;
- refuse to reset or overwrite an existing source checkout at another revision;
- document unavoidable machine prerequisites and separate external-service
  submissions from local setup.

## Rationale

The previous Week 35 entry point began with already-downloaded Runs N' Poses
files and already-installed Nesso-1 and Boltz-2 environments. It could validate
the analysis on the original workstation, but it could not bootstrap a new
computer. Making setup a required, auditable stage closes that reproducibility
gap without committing checkpoints or large datasets to Git.
