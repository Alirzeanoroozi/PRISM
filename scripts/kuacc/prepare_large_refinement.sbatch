#!/usr/bin/env bash
# Lightweight Slurm setup job: freeze and shard the existing candidate table.
set -euo pipefail

: "${RUN_ROOT:?set RUN_ROOT}"
: "${SOURCE_ROOT:?set SOURCE_ROOT}"
: "${SCRIPT_ROOT:?set SCRIPT_ROOT}"
: "${PYTHON:?set PYTHON}"

mkdir -p "${RUN_ROOT}/manifest/shards"
"${PYTHON}" "${SCRIPT_ROOT}/select_high_confidence_candidates.py" \
  --candidate-csv "${SOURCE_ROOT}/candidate_generated.csv" \
  --case-results-csv "${SOURCE_ROOT}/case_results.csv" \
  --output-csv "${RUN_ROOT}/manifest/selected_candidates.csv" \
  --selection-manifest "${RUN_ROOT}/manifest/selection_manifest.json" \
  --per-pipeline 100000 \
  --skip-input-hashes \
  --omit-selected-rows \
  --skip-path-validation

"${PYTHON}" "${SCRIPT_ROOT}/shard_refinement_manifest.py" \
  --input-csv "${RUN_ROOT}/manifest/selected_candidates.csv" \
  --output-dir "${RUN_ROOT}/manifest/shards" \
  --shard-size 1000

sha256sum \
  "${SOURCE_ROOT}/candidate_generated.csv" \
  "${SOURCE_ROOT}/case_results.csv" \
  "${RUN_ROOT}/manifest/selected_candidates.csv" \
  "${RUN_ROOT}/manifest/selection_manifest.json" \
  "${RUN_ROOT}/manifest/shards/shard_manifest.json" \
  > "${RUN_ROOT}/manifest/input_and_manifest.sha256"
