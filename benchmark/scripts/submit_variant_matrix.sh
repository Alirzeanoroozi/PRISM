#!/usr/bin/env bash
# Submit one reproducible PRISM current-pipeline variant.

set -euo pipefail

repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
variant="${VARIANT:?set VARIANT to tm_external, gt_external, tm_pyro, or gt_pyro}"
batch_root="${BATCH_ROOT:?absolute or repo-relative batch root is required}"
staged_pdb_dir="${STAGED_PDB_DIR:?absolute or repo-relative staged PDB root is required}"
output_root="${OUTPUT_ROOT:?absolute or repo-relative output root is required}"
partition="${PARTITION:?set PARTITION, for example kutem or v100_ai}"
account="${ACCOUNT:?set ACCOUNT, for example kutem or ai}"
qos="${QOS:?set QOS, for example kutem or v100_ai}"
array_spec="${ARRAY_SPEC:-1-1}"
pipeline_python="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
gtalign_bin_dir="$(dirname "${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}")"
template_workflow="${TEMPLATE_WORKFLOW:-json}"
tm_threshold="${PRISM_TM_SCORE_THRESHOLD:-0.4}"

case "$variant" in
  tm_external) aligner=tmalign; refiner=external_rosetta; default_gtalign=gtalign_cpu ;;
  gt_external) aligner=gtalign; refiner=external_rosetta; default_gtalign="$gtalign_bin_dir/gtalign_gpu" ;;
  tm_pyro) aligner=tmalign; refiner=pyrosetta; default_gtalign=gtalign_cpu ;;
  gt_pyro) aligner=gtalign; refiner=pyrosetta; default_gtalign="$gtalign_bin_dir/gtalign_gpu" ;;
  *) echo "unknown VARIANT=$variant" >&2; exit 2 ;;
esac

case "$tm_threshold" in
  ''|*[!0-9.]*) echo "invalid PRISM_TM_SCORE_THRESHOLD=$tm_threshold" >&2; exit 2 ;;
esac
# Reject malformed floats such as "0..4" or ".4.5"
if ! LC_NUMERIC=C awk "BEGIN{exit(!($tm_threshold+0==$tm_threshold && $tm_threshold>0))}" 2>/dev/null; then
  echo "invalid (malformed) PRISM_TM_SCORE_THRESHOLD=$tm_threshold" >&2; exit 2
fi

absolute_from_repo() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s/%s\n' "$repo_root" "$1" ;;
  esac
}

batch_root="$(absolute_from_repo "$batch_root")"
staged_pdb_dir="$(absolute_from_repo "$staged_pdb_dir")"
output_root="$(absolute_from_repo "$output_root")"
mkdir -p "$output_root"

if [[ "$aligner" == gtalign ]]; then
  export_gtalign="${GTALIGN_PATH:-$default_gtalign}"
else
  export_gtalign="${GTALIGN_PATH:-$default_gtalign}"
  # For TM-align variants, clear the GTALIGN_PATH override to avoid
  # passing a GTalign binary path that prism.py would ignore.
  unset GTALIGN_PATH
fi
gres_args=()
if [[ "$aligner" == gtalign && ( "$partition" == "ai" || "$partition" == "v100_ai" || "$partition" == "t4_ai" ) ]]; then
  gres_args+=(--gres=gpu:1)
fi

run_root="$output_root/runs/$variant"
mkdir -p "$run_root"
command_file="$output_root/submit_${variant}.command"
printf '%s\n' "variant=$variant" "aligner=$aligner" "refiner=$refiner" \
  "partition=$partition" "array_spec=$array_spec" > "$command_file"

job_id="$(sbatch --parsable \
  --job-name="prism-${variant}" \
  --array="$array_spec" \
  --partition="$partition" --account="$account" --qos="$qos" \
  "${gres_args[@]}" \
  --output="$output_root/slurm-%A_%a.out" \
  --error="$output_root/slurm-%A_%a.err" \
  --export="ALL,PIPELINE=current,BATCH_ROOT=$batch_root,STAGED_PDB_DIR=$staged_pdb_dir,RUN_ROOT=$run_root,TEMPLATE_WORKFLOW=$template_workflow,PRISM_TM_SCORE_THRESHOLD=$tm_threshold,ALIGNER=$aligner,REFINER=$refiner,PIPELINE_PYTHON=$pipeline_python,GTALIGN_PATH=$export_gtalign" \
  "$repo_root/benchmark/scripts/submit_comparison_batches.sbatch")"

cat > "$output_root/${variant}.manifest.tsv" <<EOF
field	value
variant	$variant
job_id	$job_id
pipeline	current
aligner	$aligner
refiner	$refiner
template_workflow	$template_workflow
tm_score_threshold	$tm_threshold
partition	$partition
account	$account
qos	$qos
array_spec	$array_spec
pipeline_python	$pipeline_python
gtalign_path	$export_gtalign
batch_root	$batch_root
staged_pdb_dir	$staged_pdb_dir
run_root	$run_root
EOF
printf '%s\n' "$job_id"
