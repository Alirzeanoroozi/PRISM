#!/bin/bash
# Score the completed Benchmark 5.5 current-pipeline array.
# Submit with --dependency=afterany:<generation_job_id>.
# The score stage remains fail-closed: invalid/missing models stay in the
# manifest and receive explicit score errors rather than being dropped.
#SBATCH --job-name=prism-b55-score
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=72:00:00
#SBATCH --output=tmp/agent/20260715-benchmark55-full/score-%j.out
#SBATCH --error=tmp/agent/20260715-benchmark55-full/score-%j.err

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
export DOCKQ_N_CPU="${DOCKQ_N_CPU:-${SLURM_CPUS_PER_TASK:-1}}"

batch_root="${BATCH_ROOT:-tmp/agent/20260715-benchmark55-full/batches}"
current_root="${CURRENT_ROOT:-tmp/agent/20260715-benchmark55-full/runs/current}"
legacy_root="${LEGACY_ROOT:-tmp/agent/20260715-implementation-smoke-legacy-final/runs/legacy}"
output_root="${OUTPUT_ROOT:-tmp/agent/20260715-benchmark55-full/scoring}"
score_python="${SCORE_PYTHON:-/scratch/tmp/prism-dockq-env/bin/python}"
pipeline_python="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"

mkdir -p "$output_root/dockq"
"$pipeline_python" benchmark/scripts/score_comparison_models.py \
  --manifest "$batch_root/shared_manifest.csv" \
  --current-root "$current_root" \
  --legacy-root "$legacy_root" \
  --output "$output_root/model_manifest.csv" \
  --score-output "$output_root/scored_models.csv" \
  --score-python "$score_python" \
  --dockq-json-dir "$output_root/dockq" \
  --n-cpu "$DOCKQ_N_CPU"

date -Is > "$output_root/completed"
