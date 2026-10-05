#!/bin/bash
# Replay common transformation/audit and compact completed USalign batches.
# Submit after the provider array, for example:
# sbatch --dependency=afterany:<provider_job> --array=1-26%4 ... this-file
#
# This job never launches USalign. It consumes retained alignment JSON from a
# completed batch, writes compact audit/count tables, validates them, and then
# removes only run-scoped raw alignment files.
#SBATCH --job-name=prism-usalign-compact
#SBATCH --array=1-26
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

set -euo pipefail

repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
run_root="${RUN_ROOT:?RUN_ROOT is required}/current/batch_$(printf '%04d' "$SLURM_ARRAY_TASK_ID")"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
mkdir -p "$run_root/status" "$run_root/logs"

cd "$repo_root"
export PRISM_FILTER_MODE="${PRISM_FILTER_MODE:-geometry_only_experimental}"
"$python_bin" benchmark/scripts/replay_compact_usalign_batch.py --run-root "$run_root" \
  > "$run_root/logs/replay_compact.out" \
  2> "$run_root/logs/replay_compact.err"
