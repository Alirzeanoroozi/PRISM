#!/bin/bash
# Score a retained positive-pose staging manifest with the canonical scorer.
#SBATCH --job-name=prism-positive-score
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=tmp/agent/20260718-prism-verification/logs/positive-score-%j.out
#SBATCH --error=tmp/agent/20260718-prism-verification/logs/positive-score-%j.err

set -euo pipefail
cd "${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
python_bin="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
score_python="${SCORE_PYTHON:-/scratch/tmp/prism-dockq-env/bin/python}"
manifest="${STAGE_MANIFEST:?STAGE_MANIFEST is required}"
output_root="${OUTPUT_ROOT:?OUTPUT_ROOT is required}"
main_root="${MAIN_ROOT:-/scratch/rshadi25/GitHub/PRISM-main}"
mkdir -p "$output_root/raw_dockq"
for set_name in rigid medium difficult; do
  "$python_bin" benchmark/scripts/score_bijective_benchmark_models.py \
    --stage-manifest "$manifest" \
    --benchmark-set "$set_name" \
    --native-root "benchmark/prism_processed/results/native_bound_complexes_t_$set_name" \
    --score-python "$score_python" \
    --irmsd-script "$main_root/benchmark/scripts/irmsd.py" \
    --output "$output_root/scores_${set_name}.tsv" \
    --interfaces-output "$output_root/scores_${set_name}_interfaces.tsv" \
    --raw-json-dir "$output_root/raw_dockq/$set_name" \
    --timeout 300
done
date -Is > "$output_root/completed"
