#!/usr/bin/env bash
#SBATCH --job-name=prism-candidate-compact
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/prism-kuacc-candidate-compact-%j.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/prism-kuacc-candidate-compact-%j.err

set -euo pipefail

RUN_ROOT=/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919
COMPACT_ROOT="$RUN_ROOT/compact_results_20260927"
PYTHON=/kuacc/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python3.11

export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

exec "$PYTHON" "$COMPACT_ROOT/aggregate_candidate_metrics.py" \
  --submission-manifest "$RUN_ROOT/full_submission_manifest.json" \
  --output-root "$COMPACT_ROOT/output"
