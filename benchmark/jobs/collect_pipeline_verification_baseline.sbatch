#!/usr/bin/env bash
# Merge independently hashed baseline shards after both have completed.
#SBATCH --job-name=prism-verify-collect
#SBATCH --partition=ai
#SBATCH --account=ai
#SBATCH --qos=ai
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=tmp/agent/20260718-prism-verification/logs/baseline-collect-%j.out
#SBATCH --error=tmp/agent/20260718-prism-verification/logs/baseline-collect-%j.err

set -euo pipefail

repo_root="${PRISM_REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
pipeline_python="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
output_root="${BASELINE_OUTPUT_ROOT:-$repo_root/tmp/agent/20260718-prism-verification/baseline}"

test -x "$pipeline_python"
test -d "$repo_root/tmp/agent/20260718-prism-verification/baseline-ai"
test -d "$repo_root/tmp/agent/20260718-prism-verification/baseline-cosbi"
test ! -e "$output_root"

cd "$repo_root"
"$pipeline_python" benchmark/scripts/collect_pipeline_verification_baseline.py \
  --shard tmp/agent/20260718-prism-verification/baseline-ai \
  --shard tmp/agent/20260718-prism-verification/baseline-cosbi \
  --output "$output_root"
