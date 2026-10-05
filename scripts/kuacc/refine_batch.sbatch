#!/usr/bin/env bash
# One candidate per task; override --array at submission after manifest creation.
#SBATCH --job-name=prism-kuacc-refine-batch
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-0%16
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%A_%a.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%A_%a.err

source /etc/bashrc
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
bundle_root="${PRISM_BUNDLE_ROOT:-$(cd -- "$script_dir/.." && pwd)}"
script_dir="${PRISM_WRAPPER_ROOT:-$bundle_root/kuacc_submission_20260919}"
if [[ ! -f "$script_dir/refine_one_external.py" ]]; then
  script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fi
pipeline_repo="$bundle_root/code/PRISM-prescript"
python_bin="${RUN_PYTHON:-/scratch/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python3.11}"
run_root="${RUN_ROOT:?set RUN_ROOT}"
batch_manifest="${BATCH_MANIFEST:-$run_root/downstream/refinement_batch_manifest.json}"

module load rosetta/2022.42
rosetta_prepack="$(command -v docking_prepack_protocol.static.linuxgccrelease)"
rosetta_dock="$(command -v docking_protocol.static.linuxgccrelease)"
rosetta_db="${PRISM_ROSETTA_DB:-/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database}"
[[ -x "$python_bin" ]] || { echo "missing RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -f "$batch_manifest" ]] || { echo "missing BATCH_MANIFEST=$batch_manifest" >&2; exit 2; }
[[ -d "$run_root" ]] || { echo "missing RUN_ROOT=$run_root" >&2; exit 2; }
[[ -x "$rosetta_prepack" && -x "$rosetta_dock" ]] || { echo "missing Rosetta executables" >&2; exit 2; }
[[ -d "$rosetta_db" ]] || { echo "missing PRISM_ROSETTA_DB=$rosetta_db" >&2; exit 2; }

runtime_lib="$(cd -- "$(dirname -- "$python_bin")/../lib" && pwd)"
export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PRISM_ROSETTA_PREPACK="$rosetta_prepack"
export PRISM_ROSETTA_DOCK="$rosetta_dock"
export PRISM_ROSETTA_DB="$rosetta_db"
export PYTHONNOUSERSITE=1

exec "$python_bin" "$script_dir/refine_one_external.py" \
  --run-root "$run_root" \
  --pipeline-repo "$pipeline_repo" \
  --pipeline-python "$python_bin" \
  --manifest "$batch_manifest" \
  --index "${SLURM_ARRAY_TASK_ID:?array task id is required}"
