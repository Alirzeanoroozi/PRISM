#!/bin/bash
# Replay published filtering on retained alignment outputs; no refinement.
#SBATCH --job-name=prism-filter-replay
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail
repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
inputs="${INPUTS:?set INPUTS}"
template_list="${TEMPLATE_LIST:?set TEMPLATE_LIST}"
alignment_root="${ALIGNMENT_ROOT:?set ALIGNMENT_ROOT}"
asset_root="${ASSET_ROOT:?set ASSET_ROOT}"
output="${OUTPUT:?set OUTPUT}"
python_bin="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
cd "$repo_root"
"$python_bin" benchmark/scripts/replay_published_filter.py \
  --inputs "$inputs" --template-list "$template_list" \
  --alignment-root "$alignment_root" --asset-root "$asset_root" --output "$output"
