#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$repo_root"

model="${1:-all}"
if [[ ! "$model" =~ ^(nesso1|boltz2|all)$ ]]; then
  echo "Usage: $0 [nesso1|boltz2|all]" >&2
  exit 2
fi

experiment="experiments/exp015_boltz2_nesso1_rnp_postcutoff50/experiment.yaml"
benchmark="configs/benchmarks/rnp_boltz2_nesso1_postcutoff50.yaml"
manifest="data/manifests/rnp_boltz2_nesso1_postcutoff50.json"
runner="scripts/run_with_provenance.py"
run_root="runs/exp015_boltz2_nesso1_rnp_postcutoff50"

if [[ "$model" == "nesso1" || "$model" == "all" ]]; then
  nesso_run="$run_root/nesso1/seed42/gate5"
  env -u PYTHONPATH NESSO_CACHE="$repo_root/cache/nesso" \
    /home/yusheng/programs/nesso-1.0.0/bin/python "$runner" \
    --run-dir "$nesso_run" --model nesso1 --cwd "$repo_root" \
    --experiment-config "$experiment" --benchmark-config "$benchmark" \
    --model-config configs/models/nesso1.yaml --manifest "$manifest" -- \
    env -u PYTHONPATH NESSO_CACHE="$repo_root/cache/nesso" \
    /home/yusheng/programs/nesso-1.0.0/bin/nesso predict \
    experiments/exp015_boltz2_nesso1_rnp_postcutoff50/inputs/nesso1/gate \
    --out_dir "$nesso_run/raw" --cache "$repo_root/cache/nesso" \
    --model_revision v1.0.0 --accelerator gpu --devices 1 \
    --precision bf16-mixed --recycling_steps 5 --num_workers 1 \
    --refine_protein_inference --refine_protein_cutoff 22.0 \
    --refine_protein_tokens_budget 256 --no_kernels --seed 42 --save_metadata
fi

if [[ "$model" == "boltz2" || "$model" == "all" ]]; then
  boltz_run="$run_root/boltz2_msa1024/seed42/gate5"
  env -u PYTHONPATH BOLTZ_CACHE="$repo_root/cache/boltz" \
    /home/yusheng/programs/boltz-2.2.1/bin/python "$runner" \
    --run-dir "$boltz_run" --model boltz2_msa1024_structure_only --cwd "$repo_root" \
    --experiment-config "$experiment" --benchmark-config "$benchmark" \
    --model-config configs/models/boltz2.yaml --manifest "$manifest" -- \
    env -u PYTHONPATH BOLTZ_CACHE="$repo_root/cache/boltz" \
    /home/yusheng/programs/boltz-2.2.1/bin/boltz predict \
    experiments/exp015_boltz2_nesso1_rnp_postcutoff50/inputs/boltz2/gate \
    --out_dir "$boltz_run/raw" --cache "$repo_root/cache/boltz" --model boltz2 \
    --checkpoint /home/yusheng/programs/boltz-2.2.1/checkpoints/boltz2_conf.ckpt \
    --accelerator gpu --devices 1 --recycling_steps 3 --sampling_steps 200 \
    --diffusion_samples 1 --max_parallel_samples 1 --step_scale 1.5 \
    --output_format mmcif --num_workers 1 --seed 42 \
    --max_msa_seqs 8192 --subsample_msa --num_subsampled_msa 1024 \
    --no_kernels
fi
