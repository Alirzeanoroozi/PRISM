#!/usr/bin/env bash
#SBATCH --job-name=prism-kuacc-dockq-batch
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-0%16
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%A_%a.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%A_%a.err

source /etc/bashrc
set -euo pipefail

RUN_ROOT=/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919
DOWNSTREAM_ROOT="$RUN_ROOT/downstream_20260927"
PIPELINE_REPO=/scratch/users/rshadi25/prism-new/code/PRISM-prescript
RUN_PYTHON=/scratch/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python3.11
DOCKQ_PYTHON=/kuacc/users/rshadi25/.conda/envs/prism_dockq_20260919/bin/python

[[ -x "$RUN_PYTHON" ]] || { echo "missing RUN_PYTHON=$RUN_PYTHON" >&2; exit 2; }
[[ -x "$DOCKQ_PYTHON" ]] || { echo "missing DOCKQ_PYTHON=$DOCKQ_PYTHON" >&2; exit 2; }
[[ -f "$DOWNSTREAM_ROOT/batches/batch_$(printf '%05d' "${SLURM_ARRAY_TASK_ID:?}").jsonl" ]] || {
  echo "missing batch for task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
}

export PYTHONNOUSERSITE=1
export DOCKQ_PYTHON="$DOCKQ_PYTHON"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

batch_id=$(printf '%05d' "$SLURM_ARRAY_TASK_ID")
exec "$RUN_PYTHON" "$DOWNSTREAM_ROOT/score_dockq_batch.py" \
  --batch-file "$DOWNSTREAM_ROOT/batches/batch_${batch_id}.jsonl" \
  --checkpoint-file "$DOWNSTREAM_ROOT/checkpoints/dockq/batch_${batch_id}.jsonl" \
  --output-csv "$DOWNSTREAM_ROOT/results/dockq_batch_${batch_id}.csv" \
  --pipeline-repo "$PIPELINE_REPO" \
  --dockq-python "$DOCKQ_PYTHON"
