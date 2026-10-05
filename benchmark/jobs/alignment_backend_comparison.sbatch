#!/bin/bash
#SBATCH --job-name=prism-align-compare
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=tmp/agent/20260714-pyrosetta-gtalign-comparison/slurm-%A.out
#SBATCH --error=tmp/agent/20260714-pyrosetta-gtalign-comparison/slurm-%A.err

set -u
repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
run_root="${RUN_ROOT:?set RUN_ROOT}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
tmalign_bin="${TMALIGN_BIN:-$repo_root/external_tools/TMalign}"
gtalign_bin="${GTALIGN_BIN:-/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu}"
mkdir -p "$(dirname "$run_root")"
cd "$repo_root"
start_epoch=$(date +%s)
set +e
"$python_bin" benchmark/scripts/run_alignment_backend_comparison.py \
  --tmalign "$tmalign_bin" --gtalign "$gtalign_bin" --output-root "$run_root" \
  --queries "${QUERIES:-10}" --refs "${REFS:-10}" --residues "${RESIDUES:-50}" \
  --tmalign-limit-pairs "${TMALIGN_LIMIT_PAIRS:-100}"
rc=$?
set -e
end_epoch=$(date +%s)
mkdir -p "$run_root"
printf '%s\n' "$python_bin benchmark/scripts/run_alignment_backend_comparison.py --tmalign $tmalign_bin --gtalign $gtalign_bin --output-root $run_root" > "$run_root/command.txt"
cat > "$run_root/parameters.tsv" <<EOF
parameter\tvalue
python\t$python_bin
tmalign\t$tmalign_bin
gtalign\t$gtalign_bin
queries\t${QUERIES:-10}
references\t${REFS:-10}
residues\t${RESIDUES:-50}
EOF
cat > "$run_root/exit.json" <<EOF
{
  "slurm_job_id": "${SLURM_JOB_ID:-}",
  "return_code": $rc,
  "scientific_status": "$([ "$rc" -eq 0 ] && echo completed || echo failed)",
  "start_epoch": $start_epoch,
  "end_epoch": $end_epoch,
  "elapsed_seconds": $((end_epoch-start_epoch)),
  "run_root": "$run_root",
  "command_file": "$run_root/command.txt",
  "parameters_file": "$run_root/parameters.tsv"
}
EOF
exit "$rc"
