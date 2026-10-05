#!/usr/bin/env bash
# KUACC portable PRISM DockQ/iRMSD scoring wrapper.
# This stage reads a completed run root and writes only below that run root.
#SBATCH --job-name=prism-kuacc-score
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%j.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%j.err

source /etc/bashrc
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
bundle_root="${PRISM_BUNDLE_ROOT:-$(cd -- "$script_dir/.." && pwd)}"
project_root="$bundle_root/code/PRISM"
pipeline_repo="$bundle_root/code/PRISM-prescript"
python_bin="${RUN_PYTHON:-/kuacc/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python}"
dockq_python="${DOCKQ_PYTHON:-/kuacc/users/rshadi25/.conda/envs/prism_dockq_20260919/bin/python}"
run_root="${RUN_ROOT:?set RUN_ROOT to an existing completed run directory}"

helper="$project_root/notebooks/tmp_prism_all_pipelines/score_transformed_models.py"
[[ -x "$python_bin" ]] || { echo "missing executable RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -x "$dockq_python" ]] || { echo "missing executable DOCKQ_PYTHON=$dockq_python" >&2; exit 2; }
[[ -d "$pipeline_repo" ]] || { echo "missing pipeline repository: $pipeline_repo" >&2; exit 2; }
[[ -f "$helper" ]] || { echo "missing scoring wrapper: $helper" >&2; exit 2; }
[[ -d "$run_root" ]] || { echo "missing run root: $run_root" >&2; exit 2; }

runtime_lib="$(cd -- "$(dirname -- "$dockq_python")/../lib" && pwd)"
export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export DOCKQ_PYTHON="$dockq_python"
exec "$python_bin" "$helper" --run-root "$run_root" --pipeline-repo "$pipeline_repo"
