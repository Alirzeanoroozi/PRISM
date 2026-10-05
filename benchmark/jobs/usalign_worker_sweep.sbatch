#!/bin/bash
#SBATCH --job-name=prism-usalign-pilot
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=00:30:00

set -euo pipefail

repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
run_root="${RUN_ROOT:?RUN_ROOT is required}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
mkdir -p "$run_root"
cd "$repo_root"

exec "$python_bin" benchmark/scripts/run_usalign_worker_sweep.py \
  --template-list "${TEMPLATE_LIST:?TEMPLATE_LIST is required}" \
  --interface-root "${INTERFACE_ROOT:?INTERFACE_ROOT is required}" \
  --query "${QUERY_PDB:?QUERY_PDB is required}" \
  --usalign "${USALIGN_BIN:-/home/rshadi25/.conda/envs/gtalign_env/bin/USalign}" \
  --output-root "$run_root" \
  --template-limit "${TEMPLATE_LIMIT:-946}" \
  --workers "${WORKERS:-1,2,4,8,16}" \
  --configs "${CONFIGS:-default,fast}"
