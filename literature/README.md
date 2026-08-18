# Literature library

This directory is the canonical citation and reading-note library for the
project. Weekly reports cite stable keys from references.bib; they should not
maintain independent, diverging bibliographies.

## Layout

    literature/
    ├── references.bib       # canonical BibTeX records
    ├── catalog.yaml         # topic, review status, and provenance
    ├── notes/               # one reviewed Markdown note per paper
    ├── reading-log/         # chronological reading and open questions
    └── papers/              # local PDFs; ignored by Git

## Source policy

- Track bibliographic metadata, DOI or authoritative URL, and our own notes.
- Do not commit publisher PDFs or copied paper figures.
- Keep author claims separate from independent evaluations.
- Mark preprints explicitly and update their record if a peer-reviewed version
  appears.
- Record which model version, dataset release, and benchmark split a claim
  concerns.
- Prefer DOI-based citation keys that remain stable across weekly reports.
- Preserve uncertainty: an abstract-only note is not a full-paper review.

The existing annotated map in docs/literature.md remains intact. It is the
source for the initial records below and can be migrated paper-by-paper after
review, avoiding a broad rewrite of existing user notes.

## Paper-note template

Each note should include:

1. full citation and persistent link;
2. source type and peer-review status;
3. question addressed;
4. data and train/test split;
5. model inputs and outputs;
6. metrics and principal results;
7. limitations and possible leakage;
8. relevance to a named project experiment;
9. exact passages, tables, or figures to revisit.

Downloaded PDFs belong in literature/papers/ and are intentionally ignored.
Ask before downloading them.
