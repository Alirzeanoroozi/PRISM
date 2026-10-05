#!/usr/bin/env bash
# Merge validated per-batch USalign compact outputs without copying raw data.
# This job is read-only with respect to batch structures and alignment files.
# Destructive cleanup remains a separately validated gate.
#SBATCH --job-name=prism-usalign-aggregate
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

: "${REPO_ROOT:?set REPO_ROOT}"
: "${RUN_ROOT:?set RUN_ROOT}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${PIPELINE_PYTHON:?set PIPELINE_PYTHON}"

mkdir -p "${OUTPUT_DIR}" "${RUN_ROOT}/logs"
exec "${PIPELINE_PYTHON}" \
  "${REPO_ROOT}/benchmark/scripts/aggregate_usalign_batches.py" \
  --run-root "${RUN_ROOT}" \
  --output-dir "${OUTPUT_DIR}" \
  --expected-batches 26 \
  > "${RUN_ROOT}/logs/aggregate_usalign.out" \
  2> "${RUN_ROOT}/logs/aggregate_usalign.err"
