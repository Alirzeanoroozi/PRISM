#!/bin/bash
# Submit with:
#   sbatch benchmark/scripts/submit_prism_analysis_all.sbatch

########################################
# SLURM resources
########################################
#SBATCH --job-name=prism_eval_all
#SBATCH --account=ai
#SBATCH --partition=ai
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --output=benchmark/prism_processed_results/slurm_logs/prism_eval_all-%J.log
#SBATCH --error=benchmark/prism_processed_results/slurm_logs/prism_eval_all-%J.log

set -euo pipefail

REPO_ROOT="/scratch/rshadi25/GitHub/PRISM"
cd "${REPO_ROOT}"

mkdir -p benchmark/prism_processed_results/slurm_logs

# Optional overrides at submit time:
#   sbatch --export=ALL,SETS=rigid,SCORE_TIMEOUT_SEC=10 benchmark/scripts/submit_prism_analysis_all.sbatch
SETS="${SETS:-rigid,medium,difficult}"
SCORE_TIMEOUT_SEC="${SCORE_TIMEOUT_SEC:-5}"
VERBOSE_EVERY="${VERBOSE_EVERY:-50}"

ROSETTA_ROOT="benchmark/prism_processed/rosetta_output_1"
FIXED_ROOT="benchmark/prism_processed/rosetta_output_1_chainfixed"
FIX_MAP="benchmark/prism_processed_results/chain_fix_map.csv"
SCORE_PY="benchmark/prism_processed_results/prism_score_env/bin/python"

echo "[INFO] Job started on $(hostname) at $(date)"
echo "[INFO] SETS=${SETS}"
echo "[INFO] SCORE_TIMEOUT_SEC=${SCORE_TIMEOUT_SEC}"

python3 -u benchmark/scripts/fix_model_chain_names.py \
  --rosetta-root "${ROSETTA_ROOT}" \
  --fixed-root "${FIXED_ROOT}" \
  --map-csv "${FIX_MAP}"

python3 -u benchmark/scripts/analyze_prism_all_benchmarks.py \
  --rosetta-root "${ROSETTA_ROOT}" \
  --score-python "${SCORE_PY}" \
  --metrics dockq,irmsd \
  --sets "${SETS}" \
  --verbose-every "${VERBOSE_EVERY}" \
  --model-fix-map-csv "${FIX_MAP}" \
  --score-timeout-sec "${SCORE_TIMEOUT_SEC}"

echo "[INFO] Job finished at $(date)"
