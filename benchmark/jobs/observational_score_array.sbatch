#!/bin/bash
# Re-score one independent model shard. This job never regenerates models.
# Submit with --array=1-10 and REPLAY_ROOT pointing to a prepared replay root.
#SBATCH --job-name=prism-score-replay
#SBATCH --array=1-10
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=tmp/agent/20260714-observational-score-replay/slurm-%A_%a.out
#SBATCH --error=tmp/agent/20260714-observational-score-replay/slurm-%A_%a.err

set -u
repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
replay_root="${REPLAY_ROOT:?set REPLAY_ROOT}"
run_python="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
score_python="${SCORE_PYTHON:-/scratch/tmp/prism-dockq-env/bin/python}"
shard="$replay_root/shards/shard_$(printf '%02d' "$SLURM_ARRAY_TASK_ID").csv"
task_root="$replay_root/tasks/task-${SLURM_ARRAY_TASK_ID}"
mkdir -p "$task_root"
cd "$repo_root"
start_epoch=$(date +%s)
# Requeues reuse an array index, so remove only this task's derived outputs before
# replaying it. Inputs and the parent replay manifest remain untouched.
rm -f "$task_root/scored_models.csv" "$task_root/exit.json" "$task_root/command.txt" "$task_root/parameters.tsv"
output_tmp="$task_root/scored_models.csv.partial"
rm -f "$output_tmp"
mkdir -p "$task_root/dockq_json"
printf '%s\n' "$run_python benchmark.scripts.score_comparison_models.score_manifest $shard $output_tmp score_python=$score_python dockq_align=${DOCKQ_ALIGN:-0} dockq_json_dir=$task_root/dockq_json" > "$task_root/command.txt"
cat > "$task_root/parameters.tsv" <<EOF
parameter\tvalue
repo_root\t$repo_root
replay_root\t$replay_root
shard\t$shard
run_python\t$run_python
score_python\t$score_python
dockq_align\t${DOCKQ_ALIGN:-0}
slurm_job_id\t${SLURM_JOB_ID:-}
slurm_array_job_id\t${SLURM_ARRAY_JOB_ID:-}
slurm_array_task_id\t${SLURM_ARRAY_TASK_ID:-}
partition\t${SLURM_JOB_PARTITION:-kutem}
cpus_per_task\t${SLURM_CPUS_PER_TASK:-2}
mem_per_node\t${SLURM_MEM_PER_NODE:-2G}
time_limit\t${SLURM_TIMELIMIT:-00:05:00}
EOF

"$run_python" - "$shard" "$output_tmp" "$repo_root" "$score_python" "${DOCKQ_ALIGN:-0}" "$task_root/dockq_json" <<'PY'
import sys
from pathlib import Path

from benchmark.scripts.score_comparison_models import score_manifest

source, output, root = map(Path, sys.argv[1:4])
score_python = Path(sys.argv[4])
dockq_json_dir = Path(sys.argv[6])
score_manifest(source, output, score_python, root, dockq_json_dir, dockq_no_align=(sys.argv[5] != "1"))
PY
rc=$?
status=failed
if [ "$rc" -eq 0 ] && [ -f "$task_root/scored_models.csv" ]; then
  status=completed
fi
if [ "$rc" -eq 0 ] && [ -f "$output_tmp" ]; then
  mv "$output_tmp" "$task_root/scored_models.csv"
  status=completed
fi
end_epoch=$(date +%s)
output_sha256=""
if [ -f "$task_root/scored_models.csv" ]; then
  output_sha256=$(sha256sum "$task_root/scored_models.csv" | awk '{print $1}')
fi
cat > "$task_root/exit.json" <<EOF
{
  "array_job_id": "${SLURM_ARRAY_JOB_ID:-}",
  "array_task_id": "${SLURM_ARRAY_TASK_ID:-}",
  "start_epoch": $start_epoch,
  "end_epoch": $end_epoch,
  "elapsed_seconds": $((end_epoch-start_epoch)),
  "return_code": $rc,
  "scientific_status": "$status",
  "resources": {
    "partition": "${SLURM_JOB_PARTITION:-kutem}",
    "account": "${SLURM_JOB_ACCOUNT:-kutem}",
    "qos": "${SLURM_JOB_QOS:-kutem}",
    "nodes": "${SLURM_JOB_NUM_NODES:-1}",
    "ntasks": "${SLURM_NTASKS:-1}",
    "cpus_per_task": "${SLURM_CPUS_PER_TASK:-2}",
    "mem": "${SLURM_MEM_PER_NODE:-2G}",
    "time_limit": "${SLURM_TIMELIMIT:-00:05:00}"
  },
  "retry_id": "${SLURM_RESTART_COUNT:-0}",
  "input": "$shard",
  "output": "$task_root/scored_models.csv",
  "output_sha256": "$output_sha256",
  "command_file": "$task_root/command.txt",
  "parameters_file": "$task_root/parameters.tsv"
}
EOF
exit "$rc"
