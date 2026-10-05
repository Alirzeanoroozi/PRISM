#!/bin/bash
# Compare current JSON filter assets against derived legacy assets.
#SBATCH --job-name=prism-filter-parity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail
repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
modern_root="${MODERN_ROOT:-$repo_root/new_template/template}"
legacy_root="${LEGACY_ROOT:-$repo_root/tmp/agent/20260718-prism-verification/template-gate/protocol-assets-ai/retry-1}"
template_list="${TEMPLATE_LIST:-$repo_root/new_template/template/final_list.txt}"
output="${OUTPUT:?set OUTPUT to a new parity TSV}"
python_bin="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
cd "$repo_root"
"$python_bin" benchmark/scripts/compare_protocol_filter_assets.py \
  --template-list "$template_list" \
  --modern-root "$modern_root" \
  --legacy-root "$legacy_root" \
  --output "$output"
