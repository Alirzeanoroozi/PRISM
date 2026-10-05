#!/usr/bin/env bash
# Aggregate only after the case array; dependency is supplied at submission.
#SBATCH --job-name=prism-kuacc-bm55-aggregate
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%j.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%j.err

source /etc/bashrc
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
wrapper_root="${BATCH_WRAPPER_ROOT:?set BATCH_WRAPPER_ROOT to the staged wrapper directory}"
manifest="${BATCH_MANIFEST:?set BATCH_MANIFEST}"
output_root="${OUTPUT_ROOT:?set OUTPUT_ROOT}"
python_bin="${RUN_PYTHON:-/scratch/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python}"
[[ -x "$python_bin" ]] || { echo "missing RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -f "$manifest" ]] || { echo "missing BATCH_MANIFEST=$manifest" >&2; exit 2; }
[[ -f "$wrapper_root/aggregate_cpu_case_batches.py" ]] || { echo "missing aggregator under BATCH_WRAPPER_ROOT=$wrapper_root" >&2; exit 2; }
export PYTHONNOUSERSITE=1
runtime_lib="$(cd -- "$(dirname -- "$python_bin")/../lib" && pwd)"
export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

exec "$python_bin" "$wrapper_root/aggregate_cpu_case_batches.py" \
  --batch-manifest "$manifest" \
  --output-root "$output_root"
