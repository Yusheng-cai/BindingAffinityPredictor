# Week 35 table-processing environment

This small environment reads the 376 MB Runs N' Poses Parquet similarity table
and deterministically reconstructs the frozen 50-system manifest. It is created
by `weeks/2026-W35/code/00_setup.sh` under
`/home/yusheng/programs/binding-affinity-week35-data`.

Pinned runtime:

- Python 3.9.23;
- NumPy 1.26.4;
- pandas 2.3.3;
- PyArrow 19.0.0.

The setup script validates these four scientific table-processing packages at
runtime. The first independent reconstruction produced a manifest with SHA-256
`ac98b0d4a0a632d4643beb8d999de59b2c8561c31538ab6565e95b145312e5ed`,
identical to `data/manifests/rnp_boltz2_nesso1_postcutoff50.json`.
