#!/bin/bash
#SBATCH --job-name=prism-gtalign-smoke
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=tmp/agent/20260714-pyrosetta-gtalign-comparison/prism-smoke-%A.out
#SBATCH --error=tmp/agent/20260714-pyrosetta-gtalign-comparison/prism-smoke-%A.err

set -u
repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
workspace_seed="${SMOKE_WORKSPACE_SEED:-/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-pipeline-smoke-wFpW4P}"
if [[ -n "${SMOKE_WORKSPACE:-}" ]]; then
  workspace="$SMOKE_WORKSPACE"
else
  workspace="/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-pipeline-smoke-${SLURM_JOB_ID:-manual}"
  cp -a "$workspace_seed" "$workspace"
fi
output_dir="${OUTPUT_DIR:-/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260714-pyrosetta-gtalign-comparison/prism-smoke-${SLURM_JOB_ID:-manual}}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
mkdir -p "$workspace"
cd "$repo_root"
start_epoch=$(date +%s)
write_exit() {
  local code="${1:-124}"
  mkdir -p "$output_dir"
  cat > "$output_dir/exit.json" <<EOF
{
  "slurm_job_id": "${SLURM_JOB_ID:-}",
  "return_code": $code,
  "scientific_status": "$([ "$code" -eq 0 ] && echo completed || ([ "$code" -eq 124 ] && echo timeout || echo failed))",
  "start_epoch": $start_epoch,
  "end_epoch": $(date +%s),
  "elapsed_seconds": $(( $(date +%s) - start_epoch )),
  "workspace": "$workspace",
  "output_dir": "$output_dir"
}
EOF
}
trap 'write_exit 124; exit 124' TERM INT
set +e
timeout 240 "$python_bin" -u benchmark/scripts/testcase_tmalign_vs_gtalign_new_template.py \
  --repo_root "$workspace" --inputs_csv inputs.csv \
  --surface_dir processed/surface_extraction --template_root templates \
  --queries 1fgnH --templates 1kcaCH \
  --tmalign "$repo_root/external_tools/TMalign" \
  --gtalign /home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu \
  --gtalign_dev_min_length 3 --gtalign_pre_score 0.0 \
  --gtalign_speed 0 --gtalign_refinement 3 --output_dir "$output_dir"
rc=$?
set -e
write_exit "$rc"
exit "$rc"
