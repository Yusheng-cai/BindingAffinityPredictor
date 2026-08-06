# Data Layout

- `manifests/`: tracked, compact records identifying the exact cases used by experiments.
- `raw/`: immutable downloaded source data; ignored by Git.
- `processed/`: reproducibly derived bulk data; ignored by Git.

Raw data should retain its original filenames and checksums. Processing code must write to `processed/` rather than modifying `raw/`.

