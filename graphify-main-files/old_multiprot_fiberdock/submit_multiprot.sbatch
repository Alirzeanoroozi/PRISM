#!/bin/bash
#SBATCH --job-name=multiprot-pipe
#SBATCH --partition=ai
#SBATCH --account=ai
#SBATCH --qos=ai
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=out-%x.%j.out
#SBATCH --error=out-%x.%j.err

set -euo pipefail

echo "=== MultiProt/FiberDock Pipeline ==="
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
date

eval "$(conda shell.bash hook)"
conda activate tmalignRosetta

cd /scratch/rshadi25/GitHub/PRISM-prescript/working_version/Multiprot-new/prism-fiberdock-cli

echo "=== Running MultiProt pipeline with pairlist_compare, template_default, job: multiprot-compare ==="
date

python2 prism.py pairlist_compare template_default multiprot-compare

exit_code=$?
echo "=== Pipeline finished ==="
echo "Exit code: $exit_code"
date

# === REPORT ===
echo ""
echo "============= RESULTS REPORT ============="
echo "Job: multiprot-compare"
echo ""

# Stage 1-2: Download + PreProcess
echo "--- Stage 1-2: PDB Download + PreProcess ---"
echo "PDB downloaded:"
ls jobs/multiprot-compare/pdb/ 2>/dev/null
echo ""
echo "Preprocessed files:"
ls jobs/multiprot-compare/preprocess/ 2>/dev/null
echo ""

# Stage 3: Surface Extraction
echo "--- Stage 3: Surface Extraction ---"
echo "ASA files:"
ls jobs/multiprot-compare/surfaceExtract/*.asa.pdb 2>/dev/null
echo "RSA files:"
ls jobs/multiprot-compare/surfaceExtract/*.rsa 2>/dev/null
echo ""

# Stage 4: Structural Alignment (MultiProt)
echo "--- Stage 4: MultiProt Alignment ---"
echo "Alignment directory:"
ls jobs/multiprot-compare/alignment/*_*_* 2>/dev/null | wc -l
echo " total alignment files"
echo ""

# Check alignment quality
echo "Sample alignment contents:"
python2 -c "
import pickle, glob, os
files = sorted(glob.glob('jobs/multiprot-compare/alignment/*_*_*'))
total = len(files)
valid = sum(1 for f in files if pickle.load(open(f)) != -1)
empty = sum(1 for f in files if pickle.load(open(f)) == -1)
print('  Total alignment pickles: %d' % total)
print('  Valid (found solutions): %d' % valid)
print('  Empty (no solution): %d' % empty)
if valid > 0:
    # Show a few samples
    import random
    samples = random.sample([f for f in files if pickle.load(open(f)) != -1], min(3, valid))
    for s in samples:
        d = pickle.load(open(s))
        basename = os.path.basename(s)
        for k in sorted(d.keys())[:3]:
            print('  %s solution %d: match_count=%d' % (basename, k, d[k][0]))
" 2>/dev/null || echo "  (check failed)"
echo ""

# Stage 5: Transformation Filtering
echo "--- Stage 5: Transformation Filtering ---"
echo "Transformation directory:"
ls jobs/multiprot-compare/transformation/ 2>/dev/null
echo ""
echo "Passed files content:"
cat jobs/multiprot-compare/transformation/passedFiles 2>/dev/null || echo "  (empty or missing)"
echo ""

# Stage 6: FiberDock Refinement
echo "--- Stage 6: FiberDock Refinement ---"
echo "FiberDock directory:"
ls jobs/multiprot-compare/fiberdock/ 2>/dev/null
echo ""
echo "Energies:"
cat jobs/multiprot-compare/fiberdock/energies/fiberdock_energies 2>/dev/null || echo "  (no energies produced)"
echo ""
echo "Structures produced:"
ls jobs/multiprot-compare/fiberdock/structures/ 2>/dev/null | head -10 || echo "  (none)"
echo ""
echo "Zero-size HB files:"
find jobs/multiprot-compare -name '*.HB' -size 0 2>/dev/null && echo "  WARNING: zero-size HB files!" || echo "  (none)"
echo ""

echo "==========================================="
exit $exit_code
