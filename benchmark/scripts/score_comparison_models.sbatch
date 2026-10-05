#!/bin/bash
#SBATCH --job-name=prism-score-compare
#SBATCH --partition=ai
#SBATCH --account=ai
#SBATCH --qos=ai
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=tmp/agent/20260712-multichain-full-comparison/score-%j.out
#SBATCH --error=tmp/agent/20260712-multichain-full-comparison/score-%j.err

set -euo pipefail
cd /scratch/rshadi25/GitHub/PRISM-prescript
export DOCKQ_N_CPU="${DOCKQ_N_CPU:-${SLURM_CPUS_PER_TASK:-1}}"

score_python="benchmark/prism_processed/env/prism_score_env/bin/python"
"$score_python" benchmark/scripts/score_comparison_models.py \
  --manifest tmp/agent/20260712-multichain-full-comparison/batches/shared_manifest.csv \
  --current-root tmp/agent/20260712-multichain-full-comparison/runs_retry2/current \
  --current-root tmp/agent/20260712-multichain-full-comparison/runs_current_ai_retry/current \
  --legacy-root tmp/agent/20260712-multichain-full-comparison/runs_legacy_ai/legacy \
  --output tmp/agent/20260712-multichain-full-comparison/model_manifest.csv \
  --score-output tmp/agent/20260712-multichain-full-comparison/scored_models.csv \
  --score-python "$score_python" \
  --n-cpu "$DOCKQ_N_CPU"
