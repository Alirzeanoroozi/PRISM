#!/usr/bin/env bash
# VALAR replacement for the canceled KUACC USalign full-panel array.
# Each task owns a resumable case ledger and an isolated attempt directory.
# This arm is alignment/transformation only: no ranking, PRODIGY, or refinement.
#SBATCH --job-name=prism-valar-usalign
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-256%24
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-valar-usalign-2cpu/slurm_logs/%x-%A_%a.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-valar-usalign-2cpu/slurm_logs/%x-%A_%a.err

source /etc/bashrc
set -euo pipefail

wrapper_root="${BATCH_WRAPPER_ROOT:?set BATCH_WRAPPER_ROOT to the staged wrapper directory}"
manifest="${BATCH_MANIFEST:?set BATCH_MANIFEST}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
[[ -x "$python_bin" ]] || { echo "missing RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -f "$manifest" ]] || { echo "missing BATCH_MANIFEST=$manifest" >&2; exit 2; }
[[ -f "$wrapper_root/run_cpu_case.py" ]] || { echo "missing worker under BATCH_WRAPPER_ROOT=$wrapper_root" >&2; exit 2; }

export PYTHONNOUSERSITE=1
export PYTHONPATH="/scratch/rshadi25/GitHub/PRISM-prescript/src${PYTHONPATH:+:$PYTHONPATH}"
runtime_lib="$(cd -- "$(dirname -- "$python_bin")/../lib" && pwd)"
export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

exec "$python_bin" "$wrapper_root/run_cpu_case.py" \
  --batch-manifest "$manifest" \
  --index "${SLURM_ARRAY_TASK_ID:?array task id is required}"
