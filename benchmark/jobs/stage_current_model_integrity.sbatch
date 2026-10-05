#!/bin/bash
# Minimal CPU-only integrity/staging validation for one retained run.
# Override BATCH_ROOT, CURRENT_ROOT, OUTPUT_ROOT, and MANIFEST at submission.
#SBATCH --job-name=prism-stage-integrity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=tmp/agent/20260718-prism-verification/logs/stage-integrity-%j.out
#SBATCH --error=tmp/agent/20260718-prism-verification/logs/stage-integrity-%j.err

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"

pipeline_python="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
batch_root="${BATCH_ROOT:?BATCH_ROOT is required}"
current_root="${CURRENT_ROOT:?CURRENT_ROOT is required}"
output_root="${OUTPUT_ROOT:?OUTPUT_ROOT is required}"
manifest="${MANIFEST:?MANIFEST is required}"
cohort_manifest="${COHORT_MANIFEST:-}"

test -x "$pipeline_python"
test -d "$batch_root"
test -d "$current_root"
mkdir -p "$(dirname "$manifest")"

set +e
stage_args=(
  --batch-root "$batch_root"
  --current-root "$current_root"
  --output-root "$output_root"
  --manifest "$manifest"
)
if [[ -n "$cohort_manifest" ]]; then
  stage_args+=(--cohort-manifest "$cohort_manifest")
fi
"$pipeline_python" benchmark/scripts/stage_current_models_for_main_benchmark.py \
  "${stage_args[@]}"
stage_status=$?
set -e

# The stager's exit status 2 means a valid empty candidate set.  It is the
# expected result for the retained cancelled v3 negative control.
if [[ "$stage_status" -eq 2 ]]; then
  exit 0
fi
exit "$stage_status"
