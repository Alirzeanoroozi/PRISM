#!/usr/bin/env bash
# Chunked corrected DockQ/iRMSD scoring for one completed GTalign stage.
# Set INPUT_MANIFEST, OUTPUT_ROOT, and STAGE at submission time.
#SBATCH --job-name=prism-gtalign-score
#SBATCH --partition=ai
#SBATCH --account=ai
#SBATCH --qos=ai
#SBATCH --array=0-38%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-gtalign-score-%x/slurm-%A_%a.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-gtalign-score-%x/slurm-%A_%a.err

set -euo pipefail
python_bin="/home/rshadi25/.conda/envs/gtalign_env/bin/python"
dockq_bin="/scratch/rshadi25/GitHub/PRISM-prescript/scripts/valar/dockq_numpy1_launcher.py"
run_root="/scratch/rshadi25/GitHub/PRISM/tmp/prism_all_pipelines/runs_gtalign/20260916T034139Z/bm55_full/full/gtalign_gpu"
pipeline_repo="/scratch/rshadi25/GitHub/PRISM-prescript"
input_manifest="${INPUT_MANIFEST:?set INPUT_MANIFEST}"
output_root="${OUTPUT_ROOT:?set OUTPUT_ROOT}"
chunk_size="${CHUNK_SIZE:-400}"
start=$((SLURM_ARRAY_TASK_ID * chunk_size))

[[ -x "$python_bin" ]] || { echo "missing executable: $python_bin" >&2; exit 2; }
[[ -x "$dockq_bin" ]] || { echo "missing executable: $dockq_bin" >&2; exit 2; }
[[ -s "$input_manifest" ]] || { echo "missing input manifest: $input_manifest" >&2; exit 2; }
mkdir -p "$output_root"
export DOCKQ_BIN="$dockq_bin"
exec "$python_bin" \
  /scratch/rshadi25/GitHub/PRISM-prescript/scripts/valar/score_gtalign_batch.py \
  --run-root "$run_root" \
  --pipeline-repo "$pipeline_repo" \
  --input-manifest "$input_manifest" \
  --output-root "$output_root" \
  --start "$start" \
  --count "$chunk_size"
