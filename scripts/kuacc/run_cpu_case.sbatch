#!/usr/bin/env bash
# One BM55 pair per task, each using the complete calculated template panel.
# The array is resumable because every task owns an attempt directory and
# atomic task_status.json; it is intentionally ranking- and PRODIGY-free.
#SBATCH --job-name=prism-kuacc-bm55-case
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-0%16
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%A_%a.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%A_%a.err

source /etc/bashrc
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
wrapper_root="${BATCH_WRAPPER_ROOT:?set BATCH_WRAPPER_ROOT to the staged wrapper directory}"
manifest="${BATCH_MANIFEST:?set BATCH_MANIFEST}"
python_bin="${RUN_PYTHON:-/scratch/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python}"
[[ -x "$python_bin" ]] || { echo "missing RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -f "$manifest" ]] || { echo "missing BATCH_MANIFEST=$manifest" >&2; exit 2; }
[[ -f "$wrapper_root/run_cpu_case.py" ]] || { echo "missing worker under BATCH_WRAPPER_ROOT=$wrapper_root" >&2; exit 2; }
export PYTHONNOUSERSITE=1
runtime_lib="$(cd -- "$(dirname -- "$python_bin")/../lib" && pwd)"
export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

exec "$python_bin" "$wrapper_root/run_cpu_case.py" \
  --batch-manifest "$manifest" \
  --index "${SLURM_ARRAY_TASK_ID:?array task id is required}"
