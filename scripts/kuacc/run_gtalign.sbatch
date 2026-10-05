#!/usr/bin/env bash
# KUACC portable PRISM GTalign-GPU wrapper.
# Generated artifacts must be placed in a new RUN_ROOT.
#SBATCH --job-name=prism-kuacc-gtalign
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%j.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-submissions/%x-%j.err

source /etc/bashrc
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
bundle_root="${PRISM_BUNDLE_ROOT:-$(cd -- "$script_dir/.." && pwd)}"
project_root="$bundle_root/code/PRISM"
pipeline_repo="$bundle_root/code/PRISM-prescript"
python_bin="${RUN_PYTHON:-/kuacc/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python}"
gtalign_bin="${GTALIGN_BIN:-$bundle_root/optional/gtalign_gpu}"
dataset_dir="${DATASET_DIR:?set DATASET_DIR}"
template_list="${TEMPLATE_LIST:-$bundle_root/templates/template_panel_full.txt}"
template_limit="${TEMPLATE_LIMIT:?set TEMPLATE_LIMIT}"
run_root="${RUN_ROOT:?set RUN_ROOT to a new run-scoped output directory}"
scenario="${SCENARIO:?set SCENARIO}"
rank="${RANK:-false}"
rank_method="${RANK_METHOD:-baseline}"
top_k="${TOP_K:-5}"

# Use the C++ runtime from the validated portable environment for GTalign.
runtime_lib="$(cd -- "$(dirname -- "$python_bin")/../lib" && pwd)"
gtalign_runtime_lib=""
if [[ -d "$(dirname -- "$gtalign_bin")/../lib" ]]; then
  gtalign_runtime_lib="$(cd -- "$(dirname -- "$gtalign_bin")/../lib" && pwd)"
fi
if [[ -n "$gtalign_runtime_lib" && "$gtalign_runtime_lib" != "$runtime_lib" ]]; then
  export LD_LIBRARY_PATH="$gtalign_runtime_lib:$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
else
  export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

helper="$project_root/notebooks/tmp_prism_all_pipelines/run_gtalign_pipeline.py"
timed_runner="$project_root/notebooks/tmp_prism_prescript_pipeline/slurm/prism_timed_runner.py"
usalign_wrapper="$project_root/notebooks/tmp_prism_all_pipelines/usalign_prism_wrapper.py"

case "$rank" in
  true|false) ;;
  *) echo "RANK must be true or false" >&2; exit 2 ;;
esac
[[ -x "$python_bin" ]] || { echo "missing executable RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -x "$gtalign_bin" ]] || { echo "missing executable GTALIGN_BIN=$gtalign_bin" >&2; exit 2; }
[[ -d "$dataset_dir" ]] || { echo "missing dataset directory: $dataset_dir" >&2; exit 2; }
[[ -s "$template_list" ]] || { echo "missing template list: $template_list" >&2; exit 2; }
[[ -f "$helper" ]] || { echo "missing GTalign wrapper: $helper" >&2; exit 2; }
[[ -f "$timed_runner" ]] || { echo "missing timed runner: $timed_runner" >&2; exit 2; }
[[ -f "$usalign_wrapper" ]] || { echo "missing USalign wrapper: $usalign_wrapper" >&2; exit 2; }
[[ -d "$pipeline_repo" ]] || { echo "missing pipeline repository: $pipeline_repo" >&2; exit 2; }
"$python_bin" -c 'import freesasa' >/dev/null 2>&1 || {
  echo "RUN_PYTHON lacks the required freesasa module; install/select a validated environment before submission" >&2
  exit 2
}

exec "$python_bin" "$helper" \
  --dataset-dir "$dataset_dir" \
  --template-list "$template_list" \
  --template-limit "$template_limit" \
  --run-root "$run_root" \
  --scenario "$scenario" \
  --pipeline-repo "$pipeline_repo" \
  --pipeline-python "$python_bin" \
  --timed-runner "$timed_runner" \
  --usalign-wrapper "$usalign_wrapper" \
  --gtalign-path "$gtalign_bin" \
  --rank "$rank" \
  --rank-method "$rank_method" \
  --top-k "$top_k"
