#!/bin/bash
# Score validated USalign transformation outputs with the common DockQ
# contract, then remove transformed halves and scoring scratch after compact
# checkpoints/TSV outputs are written.
# Submit after replay_compact_usalign_batches.sbatch.
#SBATCH --job-name=prism-usalign-dockq
#SBATCH --array=1-26
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

set -euo pipefail

repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
run_root="${RUN_ROOT:?RUN_ROOT is required}/current/batch_$(printf '%04d' "$SLURM_ARRAY_TASK_ID")"
score_python="${SCORE_PYTHON:-$repo_root/benchmark/prism_processed/env/prism_score_env/bin/python}"
native_root="${NATIVE_ROOT:-$repo_root/benchmark/prism_processed/results}"
workers="${DOCKQ_WORKERS:-4}"

mkdir -p "$run_root/status" "$run_root/logs"
cd "$repo_root"
test -x "$score_python"
export DOCKQ_N_CPU="${DOCKQ_N_CPU:-1}"
"$score_python" benchmark/scripts/score_transformed_usalign_batch.py \
  --run-root "$run_root" \
  --native-root "$native_root" \
  --score-python "$score_python" \
  --workers "$workers" \
  --retain-transformed-inputs \
  > "$run_root/logs/transformed_dockq.out" \
  2> "$run_root/logs/transformed_dockq.err"
