#!/usr/bin/env bash
# Aggregate validated common-refinement checkpoints.  This job is read-only
# with respect to candidate scratch; cleanup remains a separate validated gate.
#SBATCH --job-name=prism-refine-aggregate
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail

: "${COMPARISON_ROOT:?set COMPARISON_ROOT}"
: "${PIPELINE_REPO:?set PIPELINE_REPO}"
: "${PIPELINE_PYTHON:?set PIPELINE_PYTHON}"
: "${EXPECTED_COUNT:?set EXPECTED_COUNT}"

mkdir -p "${COMPARISON_ROOT}/aggregate" "${COMPARISON_ROOT}/logs"
exec "${PIPELINE_PYTHON}" \
  "${PIPELINE_REPO}/benchmark/scripts/aggregate_corrected_refinement.py" \
  --comparison-root "${COMPARISON_ROOT}" \
  --output-dir "${COMPARISON_ROOT}/aggregate" \
  --expected-count "${EXPECTED_COUNT}"
