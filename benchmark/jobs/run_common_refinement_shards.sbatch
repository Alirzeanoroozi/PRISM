#!/usr/bin/env bash
# Run a bounded shard of the validated common FiberDock/external-Rosetta/DockQ
# worker.  One Slurm task handles multiple candidates to avoid one scheduler
# job per candidate; each candidate still has an atomic checkpoint.
#SBATCH --job-name=prism-common-refine
#SBATCH --array=0-144%8
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00

set -euo pipefail

: "${COMPARISON_ROOT:?set COMPARISON_ROOT}"
: "${SELECTED_CSV:?set SELECTED_CSV}"
: "${WORKER_SCRIPT:?set WORKER_SCRIPT}"
: "${PIPELINE_REPO:?set PIPELINE_REPO}"
: "${FIBERDOCK_DIR:?set FIBERDOCK_DIR}"
: "${DOCKQ_PYTHON:?set DOCKQ_PYTHON}"
: "${PIPELINE_PYTHON:?set PIPELINE_PYTHON}"

module load rosetta/2022.42
rosetta_db="${ROSETTA_DB:-/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database}"
rosetta_prepack="${ROSETTA_PREPACK:-docking_prepack_protocol.static.linuxgccrelease}"
rosetta_dock="${ROSETTA_DOCK:-docking_protocol.static.linuxgccrelease}"
shard_size="${SHARD_SIZE:-100}"
start=$((SLURM_ARRAY_TASK_ID * shard_size))
total=$(( $(wc -l < "${SELECTED_CSV}") - 1 ))
if (( start >= total )); then
  exit 0
fi
end=$((start + shard_size))
if (( end > total )); then
  end=$total
fi

exec "${PIPELINE_PYTHON}" benchmark/scripts/run_common_refinement_adapter.py \
  --comparison-root "${COMPARISON_ROOT}" \
  --selected-csv "${SELECTED_CSV}" \
  --start "${start}" \
  --end "${end}" \
  --worker-script "${WORKER_SCRIPT}" \
  --pipeline-repo "${PIPELINE_REPO}" \
  --fiberdock-dir "${FIBERDOCK_DIR}" \
  --dockq-python "${DOCKQ_PYTHON}" \
  --rosetta-prepack "${rosetta_prepack}" \
  --rosetta-dock "${rosetta_dock}" \
  --rosetta-db "${rosetta_db}"
