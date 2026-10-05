#!/usr/bin/env bash
# Retry only missing PyRosetta refinements for the exact completed GTalign run.
# The retry root has symlinks to the original structures/inputs so prior
# outputs remain authoritative and are never overwritten or discarded.
#SBATCH --job-name=prism-gtalign-pyro-retry
#SBATCH --partition=ai
#SBATCH --account=ai
#SBATCH --qos=ai
#SBATCH --array=0-38%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260919-gtalign-pyrosetta-retry/slurm-%A_%a.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260919-gtalign-pyrosetta-retry/slurm-%A_%a.err

set -euo pipefail

# Use the already validated interpreter directly.  Batch jobs need not source
# a login-shell Conda hook, and the cluster hook is not installed below
# /home/rshadi25/.conda.
python_bin="/home/rshadi25/.conda/envs/gtalign_env/bin/python"

pipeline_root="/scratch/rshadi25/GitHub/PRISM"
run_root="$pipeline_root/tmp/prism_all_pipelines/runs_gtalign/20260916T034139Z/bm55_full/full/gtalign_gpu"
original_root="$run_root/processed/pyrosetta_refinement"
retry_root="$run_root/processed/pyrosetta_refinement_retry_20260919"
worker="$pipeline_root/tmp/refine_gtalign_pyrosetta.py"
chunk_size="${PYROSETTA_CHUNK_SIZE:-400}"
start=$((SLURM_ARRAY_TASK_ID * chunk_size))

[[ -f "$worker" ]] || { echo "missing worker: $worker" >&2; exit 2; }
[[ -x "$python_bin" ]] || { echo "missing executable: $python_bin" >&2; exit 2; }
[[ -d "$original_root/structures" ]] || { echo "missing original structures: $original_root/structures" >&2; exit 2; }
[[ -d "$original_root/assembled_inputs" ]] || { echo "missing original inputs: $original_root/assembled_inputs" >&2; exit 2; }

mkdir -p "$retry_root"
if [[ ! -e "$retry_root/structures" ]]; then
    ln -s ../pyrosetta_refinement/structures "$retry_root/structures"
fi
if [[ ! -e "$retry_root/assembled_inputs" ]]; then
    ln -s ../pyrosetta_refinement/assembled_inputs "$retry_root/assembled_inputs"
fi
[[ -d "$retry_root/structures" && -d "$retry_root/assembled_inputs" ]] || {
    echo "retry output links are invalid" >&2
    exit 2
}

exec "$python_bin" "$worker" \
    --run-root "$run_root" \
    --pipeline-root "$pipeline_root" \
    --output-root "$retry_root" \
    --start "$start" \
    --count "$chunk_size"
