# Nesso-1 1.0.0 environment

## Resolved installation

- Installed: 2026-08-17
- Environment: `/home/yusheng/programs/nesso-1.0.0`
- Source: `/home/yusheng/source/nesso`
- Git revision: `f0156e9a22326448684bae09ee96f73415902dcd`
- Package version: 1.0.0
- Python: 3.11.15
- PyTorch: 2.5.1+cu121
- NumPy: 2.1.3
- Lightning: 2.6.5
- Transformers: 5.15.0
- GPU/driver: RTX 3080 10 GiB / NVIDIA 535.230.02
- Cache: `/mnt/hdd2/BindingAffinityPredictor/cache/nesso`

The environment was cloned from the verified CUDA-12.1 Boltz environment to
reuse a driver-compatible PyTorch build, then isolated by uninstalling Boltz
and installing Nesso's dependencies. NumPy is pinned to 2.1.3 because Nesso
requires NumPy 2.x while the retained SciPy and Numba builds require NumPy
below 2.3 and 2.2, respectively. `pip check` reports no broken requirements.

The optional cuEquivariance kernels are not installed. All commands therefore
use `--no_kernels`; CUDA inference itself remains enabled.

## Verification

- Upstream tests: 21 passed, 7 skipped before the CCD cache existed.
- GPU tensor test: passed on the RTX 3080.
- Seed-42 tutorial prediction: completed in 23.38 seconds after assets were
  cached, with no failed examples and finite affinity outputs.
- Maximum host resident set during the smoke run: 3,927,104 KiB.

The host exports an unrelated AmberTools `PYTHONPATH`; every command must unset
it before invoking this environment. Large model, CCD, ESM, and run artifacts
remain in the ignored project cache and run directories on `/mnt/hdd2`.
