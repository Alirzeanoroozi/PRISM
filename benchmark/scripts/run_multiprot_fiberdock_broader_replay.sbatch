#!/bin/bash
#SBATCH --partition=ai
#SBATCH --qos=ai
#SBATCH --account=ai
#SBATCH --job-name=mp-fd-broad
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-fiberdock-broader-replay/broader_%j.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-fiberdock-broader-replay/broader_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=02:00:00

set -euo pipefail
cd /scratch/rshadi25/GitHub/PRISM-prescript
PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python
"$PYTHON" benchmark/scripts/run_multiprot_fiberdock_broader_replay.py \
  tmp/agent/20260721-full-validation/naccess_multiprot_rosetta \
  tmp/agent/20260722-pipeline-validation/tm_external \
  tmp/agent/20260728-multiprot-fiberdock-broader-replay/results
