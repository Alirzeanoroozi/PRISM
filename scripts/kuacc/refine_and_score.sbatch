#!/usr/bin/env bash
# Resume a completed PRISM alignment run through Rosetta and DockQ.
#SBATCH --job-name=prism-kuacc-post
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
# KUACC mid has a one-day maximum; rerun the same RUN_ROOT to resume from
# per-candidate checkpoints if this bounded post-processing job is preempted
# or reaches the walltime limit.
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%j.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%j.err

source /etc/bashrc
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
bundle_root="${PRISM_BUNDLE_ROOT:-$(cd -- "$script_dir/.." && pwd)}"
script_dir="${PRISM_WRAPPER_ROOT:-$bundle_root/kuacc_submission_20260919}"
if [[ ! -f "$script_dir/refine_and_score.py" ]]; then
  script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fi
pipeline_repo="$bundle_root/code/PRISM-prescript"
python_bin="${RUN_PYTHON:-/kuacc/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python}"
dockq_python="${DOCKQ_PYTHON:-/kuacc/users/rshadi25/.conda/envs/prism_dockq_20260919/bin/python}"
pyrosetta_python="${PYROSETTA_PYTHON:-/kuacc/users/rshadi25/.conda/envs/tmalignRosetta/bin/python}"
run_root="${RUN_ROOT:?set RUN_ROOT to an existing completed alignment run}"

module load rosetta/2022.42
rosetta_prepack="$(command -v docking_prepack_protocol.static.linuxgccrelease)"
rosetta_dock="$(command -v docking_protocol.static.linuxgccrelease)"
rosetta_db="${PRISM_ROSETTA_DB:-/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database}"

[[ -x "$python_bin" ]] || { echo "missing RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -x "$dockq_python" ]] || { echo "missing DOCKQ_PYTHON=$dockq_python" >&2; exit 2; }
[[ -d "$run_root" ]] || { echo "missing RUN_ROOT=$run_root" >&2; exit 2; }
[[ -x "$rosetta_prepack" ]] || { echo "missing Rosetta prepack executable" >&2; exit 2; }
[[ -x "$rosetta_dock" ]] || { echo "missing Rosetta docking executable" >&2; exit 2; }
[[ -d "$rosetta_db" ]] || { echo "missing PRISM_ROSETTA_DB=$rosetta_db" >&2; exit 2; }

runtime_lib="$(cd -- "$(dirname -- "$python_bin")/../lib" && pwd)"
export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PRISM_ROSETTA_PREPACK="$rosetta_prepack"
export PRISM_ROSETTA_DOCK="$rosetta_dock"
export PRISM_ROSETTA_DB="$rosetta_db"
export PYTHONNOUSERSITE=1

exec "$python_bin" "$script_dir/refine_and_score.py" \
  --run-root "$run_root" \
  --pipeline-repo "$pipeline_repo" \
  --pipeline-python "$python_bin" \
  --dockq-python "$dockq_python" \
  --pyrosetta-python "$pyrosetta_python"
