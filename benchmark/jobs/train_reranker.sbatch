#!/usr/bin/env bash
#SBATCH --job-name=prism-reranker
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

repo=${PRISM_REPO:-$PWD}
candidates=${PRISM_CANDIDATES:?Set PRISM_CANDIDATES to a labeled candidate CSV}
model=${PRISM_MODEL:-$repo/results/reranker-${SLURM_JOB_ID}.pkl}
metrics=${PRISM_METRICS:-$repo/results/reranker-${SLURM_JOB_ID}.json}
python_bin=${PRISM_TRAIN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}
seed=${PRISM_SEED:-0}

cd "$repo"
mkdir -p "$(dirname "$model")" "$(dirname "$metrics")" logs
"$python_bin" benchmark/scripts/train_reranker.py "$candidates" \
  --model "$model" --metrics "$metrics" --seed "$seed"
