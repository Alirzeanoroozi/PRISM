#!/bin/bash
# GTalign-only one-query search against the complete frozen interface tree.
# This is a backend scalability/coverage diagnostic: no TM-align, transformation,
# refinement, native scoring, or confirmatory claim is performed.
#SBATCH --job-name=prism-gtalign-full-search
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail

repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
query_pdb="${QUERY_PDB:-/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-pipeline-smoke-wFpW4P/processed/surface_extraction/1fgnH.asa.pdb}"
reference_dir="${REFERENCE_DIR:-$repo_root/new_template/template/interfaces}"
run_root="${RUN_ROOT:?set RUN_ROOT to a new output root}"
gtalign_bin="${GTALIGN_BIN:-/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu}"
mkdir -p "$run_root/query" "$run_root/output"
ln -sfn "$query_pdb" "$run_root/query/1fgnH.asa.pdb"

reference_count=$(find "$reference_dir" -maxdepth 1 -type f -name '*.pdb' | wc -l)
start_epoch=$(date +%s)
set +e
"$gtalign_bin" \
  --qrs="$run_root/query" \
  --rfs="$reference_dir" \
  -o "$run_root/output" \
  -s 0 \
  --nhits=1 \
  --nalns=1 \
  --dev-min-length=3 \
  > "$run_root/gtalign.stdout" 2> "$run_root/gtalign.stderr"
rc=$?
set -e
end_epoch=$(date +%s)
output_count=$(find "$run_root/output" -maxdepth 1 -type f -name '*.out' | wc -l)
{
  printf 'parameter\tvalue\n'
  printf 'query_pdb\t%s\n' "$query_pdb"
  printf 'reference_dir\t%s\n' "$reference_dir"
  printf 'reference_pdb_count\t%s\n' "$reference_count"
  printf 'gtalign_bin\t%s\n' "$gtalign_bin"
  printf 'nhits\t1\n'
  printf 'nalns\t1\n'
  printf 'dev_min_length\t3\n'
} > "$run_root/parameters.tsv"
cat > "$run_root/exit.json" <<EOF
{
  "slurm_job_id": "${SLURM_JOB_ID:-}",
  "return_code": $rc,
  "scientific_status": "$([ "$rc" -eq 0 ] && echo completed || echo failed)",
  "reference_pdb_count": $reference_count,
  "output_count": $output_count,
  "start_epoch": $start_epoch,
  "end_epoch": $end_epoch,
  "elapsed_seconds": $((end_epoch-start_epoch)),
  "run_root": "$run_root"
}
EOF
exit "$rc"
