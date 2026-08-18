# RFdiffusion 1.1.0 local environment

## Resolved installation

- Installed: 2026-08-11
- Environment: `/mnt/hdd2/BindingAffinityPredictor/external/envs/rfdiffusion`
- RFdiffusion source: `external/src/RFdiffusion`
- RFdiffusion revision: `86507b6538f51fce57b5a72477165f03999ed7ae`
- ProteinMPNN source: `external/src/ProteinMPNN`
- ProteinMPNN revision: `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`
- Python: 3.9.23
- PyTorch: 1.9.1, CUDA 11.1 build
- DGL: 0.9.1post1, CUDA 11.1 build
- SE(3)-Transformer: 1.0.0
- Hardware verification: RTX 3080 detected and CUDA tensor operation passed

All source, environments, checkpoints, package caches, and generated artifacts
are kept beneath the repository's ignored `external/` directory at the user's
request. They are not Git content.

## Reproducibility note

The upstream `SE3nv.yml` no longer resolves reliably with a current classic
Conda solver. An unconstrained libmamba solve also selected a CPU-only PyTorch
build. The environment therefore pins the upstream-era GPU build explicitly:

```text
pytorch=1.9.1=py3.9_cuda11.1_cudnn8.0.5_0
cudatoolkit=11.1
dgl-cuda11.1=0.9.1post1
```

The host exports an unrelated AmberTools `PYTHONPATH`. Every model command must
unset it and set `DGLBACKEND=pytorch`.

## Checkpoint

`Complex_base_ckpt.pt` is stored under `external/checkpoints/RFdiffusion/` and
has SHA-256:

```text
76e4e260aefee3b582bd76b77ab95d2592e64f00c51bf344968ab9239f3250bc
```
