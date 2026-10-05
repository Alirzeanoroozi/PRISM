#!/bin/bash
# Alignment-only contract smoke over the first bounded template panel.
# This job does not run transformation, refinement, or native scoring.
#SBATCH --job-name=prism-gtalign-panel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

set -euo pipefail

repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
seed_root="${SMOKE_WORKSPACE_SEED:-/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-pipeline-smoke-wFpW4P}"
run_root="${RUN_ROOT:?set RUN_ROOT to a new output root}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
template_limit="${TEMPLATE_LIMIT:-25}"
workspace="${run_root}/workspace"
output_dir="${run_root}/results"

mkdir -p "$run_root"
cp -a "$seed_root" "$workspace"
mkdir -p "$workspace/templates"
ln -sfn "$repo_root/new_template/template/interfaces" "$workspace/templates/interfaces"
cp "$repo_root/new_template/template/final_list.txt" "$workspace/templates/calculated_templates.txt"

cd "$repo_root"
start_epoch=$(date +%s)
set +e
"$python_bin" -u benchmark/scripts/testcase_tmalign_vs_gtalign_new_template.py \
  --repo_root "$workspace" \
  --inputs_csv inputs.csv \
  --surface_dir processed/surface_extraction \
  --template_root templates \
  --max_templates "$template_limit" \
  --queries 1fgnH \
  --tmalign "$repo_root/external_tools/TMalign" \
  --gtalign /home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu \
  --gtalign_dev_min_length 3 \
  --gtalign_pre_score 0.0 \
  --gtalign_speed 0 \
  --gtalign_refinement 3 \
  --output_dir "$output_dir"
rc=$?
set -e
end_epoch=$(date +%s)
cat > "$run_root/exit.json" <<EOF
{
  "slurm_job_id": "${SLURM_JOB_ID:-}",
  "return_code": $rc,
  "scientific_status": "$([ "$rc" -eq 0 ] && echo completed || echo failed)",
  "template_limit": $template_limit,
  "workspace": "$workspace",
  "output_dir": "$output_dir",
  "start_epoch": $start_epoch,
  "end_epoch": $end_epoch,
  "elapsed_seconds": $((end_epoch-start_epoch))
}
EOF
exit "$rc"
