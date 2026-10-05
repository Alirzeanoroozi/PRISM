#!/bin/bash
#SBATCH --job-name=prism-pyrosetta
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --array=1-10
#SBATCH --output=tmp/agent/20260714-pyrosetta-gtalign-comparison/slurm-%A_%a.out
#SBATCH --error=tmp/agent/20260714-pyrosetta-gtalign-comparison/slurm-%A_%a.err

set -u
repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
manifest="${TASK_MANIFEST:?set TASK_MANIFEST}"
run_root="${RUN_ROOT:?set RUN_ROOT}"
python_bin="${PYROSETTA_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
init_options="${PYROSETTA_INIT_OPTIONS:--mute all -constant_seed -jran 12345}"
task_id="${SLURM_ARRAY_TASK_ID:?array task required}"
task_root="$run_root/task-$task_id"
mkdir -p "$task_root"
cd "$repo_root"
start_epoch=$(date +%s)
line=$("$python_bin" - "$manifest" "$task_id" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    rows=list(csv.DictReader(handle))
row=rows[int(sys.argv[2])-1]
print("\t".join((row["input_pdb"], row.get("partners", "A_B"))))
PY
)
IFS=$'\t' read -r input_pdb partners <<< "$line"
output_pdb="$task_root/refined.pdb"
result_json="$task_root/result.json"
printf '%s\n' "$python_bin benchmark/scripts/run_pyrosetta_refinement.py --input-pdb $input_pdb --output-pdb $output_pdb --partners $partners --init-options $init_options --json $result_json --require-available" > "$task_root/command.txt"
{
  printf 'parameter\tvalue\n'
  printf 'input_pdb\t%s\n' "$input_pdb"
  printf 'partners\t%s\n' "$partners"
  printf 'python\t%s\n' "$python_bin"
  printf 'init_options\t%s\n' "$init_options"
  printf 'output_pdb\t%s\n' "$output_pdb"
  printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-}"
} > "$task_root/parameters.tsv"
set +e
"$python_bin" benchmark/scripts/run_pyrosetta_refinement.py --input-pdb "$input_pdb" --output-pdb "$output_pdb" --partners "$partners" --init-options "$init_options" --json "$result_json" --require-available
rc=$?
set -e
end_epoch=$(date +%s)
cat > "$task_root/exit.json" <<EOF
{
  "array_task_id": $task_id,
  "return_code": $rc,
  "scientific_status": "$([ "$rc" -eq 0 ] && echo success || ([ "$rc" -eq 3 ] && echo unavailable || echo failed))",
  "start_epoch": $start_epoch,
  "end_epoch": $end_epoch,
  "elapsed_seconds": $((end_epoch-start_epoch)),
  "input_pdb": "$input_pdb",
  "output_pdb": "$output_pdb",
  "result_json": "$result_json",
  "command_file": "$task_root/command.txt",
  "parameters_file": "$task_root/parameters.tsv"
}
EOF
exit "$rc"
