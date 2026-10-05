#!/usr/bin/env bash
# One manifest row per task.  The candidate worker is deliberately failure
# tolerant: stage failures are checkpointed and the worker exits successfully
# after recording them, so other array elements continue independently.
set -euo pipefail

: "${COMPARISON_ROOT:?set COMPARISON_ROOT}"
: "${SHARD_DIR:?set SHARD_DIR}"
: "${SHARD_NUMBER:?set SHARD_NUMBER}"
: "${PIPELINE_REPO:?set PIPELINE_REPO}"
: "${FIBERDOCK_DIR:?set FIBERDOCK_DIR}"
: "${DOCKQ_PYTHON:?set DOCKQ_PYTHON}"
: "${PIPELINE_PYTHON:?set PIPELINE_PYTHON}"
: "${WORKER_SCRIPT:?set WORKER_SCRIPT}"

module load rosetta/2022.42
rosetta_db="${ROSETTA_DB:-/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database}"
rosetta_prepack="${ROSETTA_PREPACK:-docking_prepack_protocol.static.linuxgccrelease}"
rosetta_dock="${ROSETTA_DOCK:-docking_protocol.static.linuxgccrelease}"
shard_path="${SHARD_DIR}/shard_$(printf '%04d' "${SHARD_NUMBER}").csv"

exec "${PIPELINE_PYTHON}" \
  "${WORKER_SCRIPT}" \
  --comparison-root "${COMPARISON_ROOT}" \
  --selected-csv "${shard_path}" \
  --index "${SLURM_ARRAY_TASK_ID:?array task id is required}" \
  --pipeline-repo "${PIPELINE_REPO}" \
  --fiberdock-dir "${FIBERDOCK_DIR}" \
  --dockq-python "${DOCKQ_PYTHON}" \
  --rosetta-prepack "${rosetta_prepack}" \
  --rosetta-dock "${rosetta_dock}" \
  --rosetta-db "${rosetta_db}"
