#!/bin/bash
# Score one completed variant stage manifest with the canonical chain-aware evaluator.
# Submit on a compute partition; this job is independent of model generation once
# the stage manifest and native root are immutable.
#SBATCH --job-name=prism-variant-score
#SBATCH --partition=ai
#SBATCH --account=ai
#SBATCH --qos=ai
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"

: "${STAGE_MANIFEST:?stage manifest is required}"
: "${NATIVE_ROOT:?native root is required}"
: "${OUTPUT:?output CSV is required}"
: "${RAW_JSON_DIR:?raw DockQ directory is required}"
score_python="${SCORE_PYTHON:-/scratch/tmp/prism-dockq-env/bin/python}"
pipeline_python="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
main_root="${MAIN_ROOT:-/scratch/rshadi25/GitHub/PRISM-main}"

mkdir -p "$(dirname "$OUTPUT")" "$RAW_JSON_DIR"
exec "$pipeline_python" benchmark/scripts/score_bijective_benchmark_models.py \
  --stage-manifest "$STAGE_MANIFEST" \
  --native-root "$NATIVE_ROOT" \
  --score-python "$score_python" \
  --irmsd-script "$main_root/benchmark/scripts/irmsd.py" \
  --output "$OUTPUT" \
  --raw-json-dir "$RAW_JSON_DIR" \
  --timeout "${SCORE_TIMEOUT:-300}"
