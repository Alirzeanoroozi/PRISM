#!/bin/bash
#SBATCH --partition=ai
#SBATCH --qos=ai
#SBATCH --account=ai
#SBATCH --job-name=fd-corrected
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-fiberdock-corrected-replay/corrected_%j.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-fiberdock-corrected-replay/corrected_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=01:00:00

set -euo pipefail
cd /scratch/rshadi25/GitHub/PRISM-prescript
PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python
"$PYTHON" benchmark/scripts/run_fiberdock_corrected_replay.py \
  tmp/agent/20260728-fiberdock-corrected-replay/results
