#!/usr/bin/env bash
# One selected transformed candidate per task; all outputs stay in the
# run-scoped comparison root and are checkpointed by the Python worker.
# Required environment: COMPARISON_ROOT, SELECTED_CSV, PIPELINE_REPO,
# FIBERDOCK_DIR, DOCKQ_PYTHON, PIPELINE_PYTHON.
set -euo pipefail

: "${COMPARISON_ROOT:?set COMPARISON_ROOT}"
: "${SELECTED_CSV:?set SELECTED_CSV}"
: "${PIPELINE_REPO:?set PIPELINE_REPO}"
: "${FIBERDOCK_DIR:?set FIBERDOCK_DIR}"
: "${DOCKQ_PYTHON:?set DOCKQ_PYTHON}"
: "${PIPELINE_PYTHON:?set PIPELINE_PYTHON}"
: "${WORKER_SCRIPT:?set WORKER_SCRIPT}"

module load rosetta/2022.42
rosetta_db="${ROSETTA_DB:-/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database}"
rosetta_prepack="${ROSETTA_PREPACK:-docking_prepack_protocol.static.linuxgccrelease}"
rosetta_dock="${ROSETTA_DOCK:-docking_protocol.static.linuxgccrelease}"

exec "${PIPELINE_PYTHON}" \
  "${WORKER_SCRIPT}" \
  --comparison-root "${COMPARISON_ROOT}" \
  --selected-csv "${SELECTED_CSV}" \
  --index "${SLURM_ARRAY_TASK_ID:-${CANDIDATE_INDEX:?set CANDIDATE_INDEX outside an array}}" \
  --pipeline-repo "${PIPELINE_REPO}" \
  --fiberdock-dir "${FIBERDOCK_DIR}" \
  --dockq-python "${DOCKQ_PYTHON}" \
  --rosetta-prepack "${rosetta_prepack}" \
  --rosetta-dock "${rosetta_dock}" \
  --rosetta-db "${rosetta_db}"
