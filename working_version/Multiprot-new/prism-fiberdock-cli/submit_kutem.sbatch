#!/bin/bash
#SBATCH --job-name=prism-fiberdock-10pairs
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=out-%x.%j.out
#SBATCH --error=out-%x.%j.err

set -euo pipefail

echo "=== Job started ==="
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CWD: $(pwd)"
date

# --- environment ---------------------------------------------------------
eval "$(conda shell.bash hook)"
conda activate tmalignRosetta

# --- run from the prism-fiberdock-cli directory --------------------------
cd /scratch/rshadi25/GitHub/PRISM-prescript/working_version/Multiprot-new/prism-fiberdock-cli

echo "=== Running prism.py with pairlist_10, template_default, job: 10pairs-kutem ==="
date

python2 prism.py pairlist_10 template_default 10pairs-kutem

exit_code=$?
echo "=== Job finished ==="
echo "Exit code: $exit_code"
date

# Verify outputs
echo "=== Output verification ==="
find jobs/10pairs-kutem -name '*.HB' -size 0 2>/dev/null && echo "WARNING: Zero-size HB files found!" || echo "No zero-size HB files (good)"
echo "=== Output structure ==="
find jobs/10pairs-kutem -type f 2>/dev/null | head -60
echo "=== Energy results ==="
cat jobs/10pairs-kutem/fiberdock/energies/fiberdock_energies 2>/dev/null || echo "No fiberdock energies produced"

exit $exit_code
