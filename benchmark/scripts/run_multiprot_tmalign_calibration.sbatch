#!/bin/bash
#SBATCH --partition=ai
#SBATCH --qos=ai
#SBATCH --account=ai
#SBATCH --job-name=mp-tm-calibration
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-tmalign-calibration/calibration_%j.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-tmalign-calibration/calibration_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=02:00:00

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python
"$PYTHON" /scratch/rshadi25/GitHub/PRISM-prescript/benchmark/scripts/run_multiprot_tmalign_calibration.py \
  --panel-root /scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260722-pipeline-validation/multiprot_external \
  --output-root /scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-tmalign-calibration/results \
  --workers "$SLURM_CPUS_PER_TASK"
