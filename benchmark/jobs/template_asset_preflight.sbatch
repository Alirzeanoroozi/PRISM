#!/bin/bash
# Hash and validate the frozen template panel; no source assets are modified.
#SBATCH --job-name=prism-template-preflight
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=tmp/agent/20260718-prism-verification/logs/template-preflight-%j.out
#SBATCH --error=tmp/agent/20260718-prism-verification/logs/template-preflight-%j.err

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
python_bin="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
output_root="${OUTPUT_ROOT:?OUTPUT_ROOT is required}"
template_manifest="${TEMPLATE_MANIFEST:-new_template/template/final_list.txt}"
asset_root="${ASSET_ROOT:-new_template/template}"
mkdir -p "$output_root"
"$python_bin" benchmark/scripts/investigation_provenance.py \
  --repo-root . \
  --output "$output_root/template_preflight.json" \
  --template-manifest "$template_manifest" \
  --asset-root "$asset_root" \
  --template-tsv "$output_root/template_assets.tsv"
date -Is > "$output_root/completed"
