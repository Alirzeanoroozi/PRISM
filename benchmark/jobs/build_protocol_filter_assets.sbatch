#!/bin/bash
# Derive a run-local JSON filter tree from retained legacy text assets.
#SBATCH --job-name=prism-protocol-assets
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=tmp/agent/20260718-prism-verification/logs/protocol-assets-%j.out
#SBATCH --error=tmp/agent/20260718-prism-verification/logs/protocol-assets-%j.err

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
python_bin="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
output_root="${OUTPUT_ROOT:?OUTPUT_ROOT is required}"
"$python_bin" benchmark/scripts/build_protocol_filter_assets.py \
  --template-list new_template/template/final_list.txt \
  --legacy-root template_old/template \
  --output-root "$output_root"
date -Is > "$output_root/completed"
