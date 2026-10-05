#!/usr/bin/env bash
# KUACC portable PRISM CPU wrapper.
# The bundle and all generated artifacts are selected through explicit paths;
# the canonical PRISM checkouts are never written by this wrapper.
#SBATCH --job-name=prism-kuacc-cpu
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
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
aligner="${ALIGNER:-tmalign}"
dataset_dir="${DATASET_DIR:?set DATASET_DIR}"
template_list="${TEMPLATE_LIST:-$bundle_root/templates/template_panel_full.txt}"
template_limit="${TEMPLATE_LIMIT:?set TEMPLATE_LIMIT}"
run_root="${RUN_ROOT:?set RUN_ROOT to a new run-scoped output directory}"
scenario="${SCENARIO:?set SCENARIO}"
rank="${RANK:-false}"
rank_method="${RANK_METHOD:-baseline}"
top_k="${TOP_K:-5}"

# The KUACC image exposes an older system libstdc++. The validated portable
# environment contains the C++ runtime required by the bundled aligners.
runtime_lib="$(cd -- "$(dirname -- "$python_bin")/../lib" && pwd)"
export LD_LIBRARY_PATH="$runtime_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

helper="$project_root/notebooks/tmp_prism_all_pipelines/run_cpu_pipeline.py"
timed_runner="$project_root/notebooks/tmp_prism_prescript_pipeline/slurm/prism_timed_runner.py"
usalign_wrapper="$project_root/notebooks/tmp_prism_all_pipelines/usalign_prism_wrapper.py"

case "$aligner" in
  tmalign|multiprot|usalign) ;;
  *) echo "unsupported ALIGNER=$aligner" >&2; exit 2 ;;
esac
case "$rank" in
  true|false) ;;
  *) echo "RANK must be true or false" >&2; exit 2 ;;
esac
[[ -x "$python_bin" ]] || { echo "missing executable RUN_PYTHON=$python_bin" >&2; exit 2; }
[[ -d "$dataset_dir" ]] || { echo "missing dataset directory: $dataset_dir" >&2; exit 2; }
[[ -s "$template_list" ]] || { echo "missing template list: $template_list" >&2; exit 2; }
[[ -f "$helper" ]] || { echo "missing CPU wrapper: $helper" >&2; exit 2; }
[[ -f "$timed_runner" ]] || { echo "missing timed runner: $timed_runner" >&2; exit 2; }
[[ -f "$usalign_wrapper" ]] || { echo "missing USalign wrapper: $usalign_wrapper" >&2; exit 2; }
[[ -d "$pipeline_repo" ]] || { echo "missing pipeline repository: $pipeline_repo" >&2; exit 2; }
case "$aligner" in
  tmalign) [[ -x "$pipeline_repo/external_tools/TMalign" ]] || { echo "missing TMalign executable" >&2; exit 2; } ;;
  multiprot) [[ -x "$pipeline_repo/external_tools/multiprot.Linux" ]] || { echo "missing MultiProt executable" >&2; exit 2; } ;;
  usalign) [[ -x "$pipeline_repo/external_tools/USalign" ]] || { echo "missing USalign executable" >&2; exit 2; } ;;
esac
"$python_bin" -c 'import freesasa' >/dev/null 2>&1 || {
  echo "RUN_PYTHON lacks the required freesasa module; install/select a validated environment before submission" >&2
  exit 2
}

exec "$python_bin" "$helper" \
  --aligner "$aligner" \
  --dataset-dir "$dataset_dir" \
  --template-list "$template_list" \
  --template-limit "$template_limit" \
  --run-root "$run_root" \
  --scenario "$scenario" \
  --pipeline-repo "$pipeline_repo" \
  --pipeline-python "$python_bin" \
  --timed-runner "$timed_runner" \
  --usalign-wrapper "$usalign_wrapper" \
  --rank "$rank" \
  --rank-method "$rank_method" \
  --top-k "$top_k"
