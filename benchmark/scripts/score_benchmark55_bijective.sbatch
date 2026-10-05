#!/bin/bash
# Score Benchmark 5.5 using complete chain-bijective mappings.
# Generation is intentionally separate; submit this with afterany:<array_job>.
#SBATCH --job-name=prism-b55-bijective-score
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=72:00:00
#SBATCH --output=tmp/agent/20260716-benchmark55-main-contract-r2/bijective-score-%j.out
#SBATCH --error=tmp/agent/20260716-benchmark55-main-contract-r2/bijective-score-%j.err

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
export DOCKQ_N_CPU="${DOCKQ_N_CPU:-${SLURM_CPUS_PER_TASK:-1}}"

batch_root="${BATCH_ROOT:?absolute Benchmark 5.5 batch root is required}"
current_root="${CURRENT_ROOT:?absolute current-pipeline run root is required}"
output_root="${OUTPUT_ROOT:?absolute output root is required}"
score_python="${SCORE_PYTHON:-/scratch/tmp/prism-dockq-env/bin/python}"
pipeline_python="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
main_root="${MAIN_ROOT:-/scratch/rshadi25/GitHub/PRISM-main}"
stage_root="$output_root/staged_bijective"
stage_manifest="$output_root/stage_manifest.csv"

test -x "$score_python"
test -x "$pipeline_python"
test -d "$batch_root"
test -d "$current_root"
mkdir -p "$output_root"

"$pipeline_python" benchmark/scripts/stage_current_models_for_main_benchmark.py \
  --batch-root "$batch_root" \
  --current-root "$current_root" \
  --output-root "$stage_root" \
  --manifest "$stage_manifest"

score_set() {
  local set_name="$1"
  local native_root="$2"
  "$pipeline_python" benchmark/scripts/score_bijective_benchmark_models.py \
    --stage-manifest "$stage_manifest" \
    --benchmark-set "$set_name" \
    --native-root "$native_root" \
    --score-python "$score_python" \
    --irmsd-script "$main_root/benchmark/scripts/irmsd.py" \
    --output "$output_root/scored_${set_name}.csv" \
    --raw-json-dir "$output_root/raw_dockq/$set_name" \
    --timeout 300
}

score_set rigid benchmark/prism_processed/results/native_bound_complexes_t_rigid
score_set medium benchmark/prism_processed/results/native_bound_complexes_t_medium
score_set difficult benchmark/prism_processed/results/native_bound_complexes_t_difficult

sha256sum "$main_root/benchmark/scripts/irmsd.py" benchmark/scripts/score_bijective_benchmark_models.py > "$output_root/scoring_script_hashes.txt"
printf 'score_contract=bijective_model_to_native\nDockQ=python -m DockQ with complete mapping and raw interface selection\niRMSD=PRISM-main benchmark/scripts/irmsd.py with grouped chains\n' > "$output_root/scoring_contract.txt"
date -Is > "$output_root/completed"
