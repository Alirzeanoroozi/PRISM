#!/usr/bin/env bash
# Aggregate only after every VALAR USalign case has reached a terminal ledger.
#SBATCH --job-name=prism-valar-usalign-aggregate
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-valar-usalign-2cpu/slurm_logs/%x-%j.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-valar-usalign-2cpu/slurm_logs/%x-%j.err

source /etc/bashrc
set -euo pipefail

wrapper_root="${BATCH_WRAPPER_ROOT:?set BATCH_WRAPPER_ROOT to the staged wrapper directory}"
manifest="${BATCH_MANIFEST:?set BATCH_MANIFEST}"
output_root="${BATCH_OUTPUT_ROOT:?set BATCH_OUTPUT_ROOT}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
[[ -x "$python_bin" ]] || { echo "missing RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -f "$manifest" ]] || { echo "missing BATCH_MANIFEST=$manifest" >&2; exit 2; }
[[ -f "$wrapper_root/aggregate_cpu_case_batches.py" ]] || { echo "missing aggregator under BATCH_WRAPPER_ROOT=$wrapper_root" >&2; exit 2; }

export PYTHONNOUSERSITE=1
export PYTHONPATH="/scratch/rshadi25/GitHub/PRISM-prescript/src${PYTHONPATH:+:$PYTHONPATH}"
runtime_lib="$(cd -- "$(dirname -- "$python_bin")/../lib" && pwd)"
export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

exec "$python_bin" "$wrapper_root/aggregate_cpu_case_batches.py" \
  --batch-manifest "$manifest" \
  --output-root "$output_root"
