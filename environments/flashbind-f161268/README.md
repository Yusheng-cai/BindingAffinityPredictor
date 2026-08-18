# FlashBind f161268 inference environment

## Resolved installation

- Installed: 2026-08-17
- Stable path: `/home/yusheng/programs/flashbind-f161268`
- Physical path: `/mnt/hdd2/BindingAffinityPredictor/external/envs/flashbind-f161268`
- Source: `/home/yusheng/source/FlashBind`
- Git revision: `f161268176237ab6ce5757031a8c1b93937e0d37`
- Python: 3.11.15
- PyTorch: 2.5.1+cu121
- PyTorch Lightning: 2.5.0
- PyTorch Geometric: 2.6.1
- LMDB: 1.6.2
- GPU/driver: RTX 3080 10 GiB / NVIDIA 535.230.02

The upstream environment requests PyTorch 2.7.1 with CUDA 12.6 and a much
larger preprocessing stack. For this released-pose scoring reproduction, the
environment was isolated by cloning the verified CUDA-12.1 Nesso environment
and adding only `lmdb==1.6.2`, `timeout-decorator==0.5.0`, and
`torch-geometric==2.6.1`. The precomputed archive makes ESM3 and FABind+
execution dependencies unnecessary for this experiment. This dependency
deviation is recorded because the result is not a bitwise upstream-environment
reproduction.

The physical environment resides on `/mnt/hdd2` because the system disk had
only 7.5 GiB free. The stable program path is retained through a symbolic link.
