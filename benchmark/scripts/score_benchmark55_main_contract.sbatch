#!/bin/bash
# Score Benchmark 5.5 models with the unmodified PRISM-main benchmark scripts.
# Submit with --dependency=afterany:<current-pipeline-array-job>.
# The staging step creates symlinks only; model coordinates and chains stay intact.
#SBATCH --job-name=prism-b55-main-score
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=72:00:00
#SBATCH --output=tmp/agent/20260716-benchmark55-main-contract/score-%j.out
#SBATCH --error=tmp/agent/20260716-benchmark55-main-contract/score-%j.err

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"

main_root="${MAIN_ROOT:-/scratch/rshadi25/GitHub/PRISM-main}"
batch_root="${BATCH_ROOT:-tmp/agent/20260716-benchmark55-main-contract/batches}"
current_root="${CURRENT_ROOT:-tmp/agent/20260716-benchmark55-main-contract/runs/current}"
output_root="${OUTPUT_ROOT:-tmp/agent/20260716-benchmark55-main-contract/scoring}"
score_python="${SCORE_PYTHON:-/scratch/tmp/prism-dockq-env/bin/python}"
pipeline_python="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
stage_root="$output_root/staged_main_contract"
runner_root="$output_root/main_analyzer_runner"

test -x "$score_python"
test -x "$pipeline_python"
test -f "$main_root/benchmark/scripts/rosetta_output/analyze_prism_rigid_results.py"
test -d "$batch_root"
test -d "$current_root"
mkdir -p "$output_root"

# PRISM-main's analyzer locates dockq.py and irmsd.py beside its own
# __file__, although those files live one directory above rosetta_output.
# Copy only the unmodified analyzer bytes into a disposable runner directory,
# then link the exact PRISM-main metric scripts beside it.  No metric source is
# patched and no file in PRISM-main is modified.
mkdir "$runner_root"
cp "$main_root/benchmark/scripts/rosetta_output/analyze_prism_rigid_results.py" \
  "$runner_root/analyze_prism_rigid_results.py"
cmp -s "$main_root/benchmark/scripts/rosetta_output/analyze_prism_rigid_results.py" \
  "$runner_root/analyze_prism_rigid_results.py"
ln -s "$main_root/benchmark/scripts/dockq.py" "$runner_root/dockq.py"
ln -s "$main_root/benchmark/scripts/irmsd.py" "$runner_root/irmsd.py"
ln -s "$main_root/benchmark/scripts/irmsd_backbone.py" "$runner_root/irmsd_backbone.py"
sha256sum "$main_root/benchmark/scripts/rosetta_output/analyze_prism_rigid_results.py" \
  "$runner_root/analyze_prism_rigid_results.py" > "$output_root/main_analyzer_hashes.txt"

"$pipeline_python" benchmark/scripts/stage_current_models_for_main_benchmark.py \
  --batch-root "$batch_root" \
  --current-root "$current_root" \
  --output-root "$stage_root" \
  --manifest "$output_root/stage_manifest.csv"

run_set() {
  local set_name="$1"
  local csv_path="$2"
  local native_dir="$3"
  "$score_python" "$runner_root/analyze_prism_rigid_results.py" \
    --rosetta-root "$stage_root/$set_name" \
    --t-rigid-csv "$csv_path" \
    --native-dir "$native_dir" \
    --out-dir "$output_root/$set_name" \
    --metrics dockq,irmsd \
    --score-python "$score_python" \
    --score-timeout-sec 300 \
    --verbose-every 25
}

run_set rigid benchmark/data/T_Rigid.csv benchmark/prism_processed/results/native_bound_complexes_t_rigid
run_set medium benchmark/data/T_medium.csv benchmark/prism_processed/results/native_bound_complexes_t_medium
run_set difficult benchmark/data/T_difficult.csv benchmark/prism_processed/results/native_bound_complexes_t_difficult

printf 'main_root=%s\nscore_python=%s\nbatch_root=%s\ncurrent_root=%s\n' \
  "$main_root" "$score_python" "$batch_root" "$current_root" > "$output_root/scoring_contract.txt"
date -Is > "$output_root/completed"
