#!/bin/bash
#SBATCH --partition=ai
#SBATCH --qos=ai
#SBATCH --account=ai
#SBATCH --job-name=mp-tm-ledger
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-fiberdock-broader-replay/ledger_%j.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-fiberdock-broader-replay/ledger_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=00:30:00

set -euo pipefail
cd /scratch/rshadi25/GitHub/PRISM-prescript
PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python
"$PYTHON" benchmark/scripts/compare_six_pair_arms.py \
  tmp/agent/20260728-multiprot-fiberdock-broader-replay/results/multiprot_fiberdock_broader_replay.json \
  tmp/agent/20260721-full-validation/naccess_multiprot_rosetta \
  tmp/agent/20260722-pipeline-validation/tm_external \
  tmp/agent/20260728-multiprot-fiberdock-broader-replay/results/six_pair_arm_comparison.json
