#!/bin/bash
# Minimal CPU-only completion canaries.
#SBATCH --job-name=prism-completion-canaries
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=tmp/agent/20260718-prism-verification/logs/canaries-%j.out
#SBATCH --error=tmp/agent/20260718-prism-verification/logs/canaries-%j.err

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
python_bin="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
output="${OUTPUT_ROOT:?OUTPUT_ROOT is required}"
test -x "$python_bin"
"$python_bin" benchmark/scripts/run_pipeline_completion_canaries.py --output "$output"
