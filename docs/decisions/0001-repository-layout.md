# 0001: Separate tracked experiment intent from generated artifacts

- Date: 2026-08-06
- Status: Accepted

## Context

Affinity models have large checkpoints, conflicting software dependencies, stochastic inference, heterogeneous outputs, and potentially large benchmark datasets. A repository containing copied upstream code and unstructured output directories would be difficult to audit or operate with multiple agents.

## Decision

- Keep upstream model source and installed environments external to this repository.
- Track exact external revisions and checkpoint provenance in model configurations.
- Track experiment definitions, small input manifests, reusable adapters, tests, and reviewed reports.
- Ignore downloaded data, caches, weights, MSAs, and raw run artifacts.
- Use one adapter per model and a shared canonical result representation.

## Consequences

- The Git history remains compact and scientifically interpretable.
- A run is reproducible only when its external environment and checkpoint provenance are captured correctly.
- Initial adapter and provenance plumbing requires more care than isolated shell scripts, but later comparisons become much safer.

