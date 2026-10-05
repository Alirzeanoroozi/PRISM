#!/usr/bin/env bash
# Freeze retained verification evidence without modifying source artifacts.
# Submit from the repository root after creating the output log directory.
#SBATCH --job-name=prism-verify-baseline
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=tmp/agent/20260718-prism-verification/logs/baseline-%j.out
#SBATCH --error=tmp/agent/20260718-prism-verification/logs/baseline-%j.err

set -euo pipefail

repo_root="${PRISM_REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
pipeline_python="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
output_root="${BASELINE_OUTPUT_ROOT:-$repo_root/tmp/agent/20260718-prism-verification/baseline}"
baseline_config="${BASELINE_CONFIG:-$repo_root/benchmark/configs/prism_pipeline_verification.json}"

test -x "$pipeline_python"
test -f "$baseline_config"
test ! -e "$output_root" || test -z "$(find "$output_root" -mindepth 1 -maxdepth 1 -print -quit)"

cd "$repo_root"
"$pipeline_python" benchmark/scripts/build_pipeline_verification_baseline.py \
  --config "$baseline_config" \
  --output "$output_root"
