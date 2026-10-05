#!/bin/bash
#SBATCH --partition=ai
#SBATCH --qos=ai
#SBATCH --account=ai
#SBATCH --job-name=multiprot_diag
#SBATCH --output=multiprot_diag_%j.out
#SBATCH --error=multiprot_diag_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=02:00:00

# Diagnostic MultiProt run: 28 pairs, 8 parallel workers
# Uses PRISM_MULTIPROT_FORCE=1 to bypass seccomp
# Lowered gate thresholds: TM>=0.3, matches>=10, pct>=30%

set -euo pipefail
cd /scratch/rshadi25/GitHub/PRISM-prescript

# Load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gtalign_env

export PRISM_MULTIPROT_FORCE=1
export PRISM_TM_SCORE_THRESHOLD=0.3
export PRISM_MIN_MATCH_COUNT=10
export PRISM_MIN_MATCH_PCT=30

# Run the diagnostic alignment
python3 -c "
import sys
sys.path.insert(0, 'src')
from alignment_multiprot import align_multiprot

# BM5.5 queries: 5zngA, 4eylA
# 14 templates per query = 28 pairs
queries = ['5zngA', '4eylA']
templates = [
    '1a0cCD', '1a0hAE', '1ay7AB', '1b6cAB', '1buhAB',
    '1e6eAB', '1fq1AB', '1gpwAB', '1he1AB', '1ibrAB',
    '1k4cAB', '1k74AB', '1kacAB', '1kxqAB'
]

print(f'Running MultiProt on {len(queries)} queries x {len(templates)} templates = {len(queries)*len(templates)} pairs with 8 workers...')
align_multiprot(queries, templates, max_workers=8)
print('Completed.')
"
