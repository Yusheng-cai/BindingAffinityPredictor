# Boltz-2 environment

## Resolved installation

- Installed: 2026-08-06
- Environment: `/home/yusheng/programs/boltz-2.2.1`
- Source: `/home/yusheng/source/boltz`
- Git revision: `b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc`
- Git description: `v2.2.1-6-gb1ebfc4`
- Package metadata version: `2.2.1`
- Python: 3.11.15
- PyTorch: 2.5.1+cu121
- GPU/driver: RTX 3080 10 GiB / NVIDIA 535.230.02
- Cache: `/mnt/hdd2/BindingAffinityPredictor/cache/boltz`

The exact commit is six commits after the upstream `v2.2.1` tag while retaining
package metadata version 2.2.1. The commit hash, rather than the display version,
is the authoritative software identity.

## Why this CUDA stack

An unconstrained installation on 2026-08-06 attempted to resolve PyTorch 2.13
and CUDA 13, which the installed driver cannot support. PyTorch 2.5.1 with its
CUDA 12.1 runtime was therefore installed explicitly and passed an on-GPU tensor
test.

The optional `cuequivariance` extra was not installed. Its current CUDA wheels
require a newer cuBLAS stack that conflicts with the driver-compatible PyTorch
build. Predictions must use Boltz's documented `--no_kernels` fallback. This
retains CUDA inference and replaces only optional optimized triangular kernels.

## Environment isolation

The host exports an AmberTools Python 3.9 path through `PYTHONPATH`. Allowing it
into this Python 3.11 environment exposes unrelated packages and produces a
false dependency failure. Every Boltz command must therefore begin with:

```bash
env -u PYTHONPATH \
  BOLTZ_CACHE=/mnt/hdd2/BindingAffinityPredictor/cache/boltz \
  /home/yusheng/programs/boltz-2.2.1/bin/boltz
```

For example, the future exp001 invocation will append `predict`, the input and
output paths, all frozen inference arguments, `--use_msa_server`, and
`--no_kernels`. It is not executed as part of environment installation.

## Verification completed

- `pip check`: no broken requirements when the external `PYTHONPATH` is unset.
- `boltz predict --help`: imports successfully and exposes the required flags.
- PyTorch sees the RTX 3080 and completed a CUDA tensor operation.
- Both checkpoints deserialize on CPU and contain state dictionaries and
  hyperparameters.

`resolved-packages.txt` records the complete Python package set.
`conda-explicit-linux-64.txt` records the base Conda artifacts.
