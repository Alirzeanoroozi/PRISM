#!/bin/bash
#SBATCH --job-name=tmalign-pipe
#SBATCH --partition=ai
#SBATCH --account=ai
#SBATCH --qos=ai
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=out-%x.%j.out
#SBATCH --error=out-%x.%j.err

set -euo pipefail

echo "=== TMalign/Rosetta Pipeline ==="
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
date

# First: ensure template library is complete (copy from MultiProt if needed)
echo "=== Ensuring template library ==="
TEMPLATE_SRC="/scratch/rshadi25/GitHub/PRISM-prescript/working_version/Multiprot-new/prism-fiberdock-cli/template"
TEMPLATE_DST="/scratch/rshadi25/GitHub/PRISM-prescript/working_version/TMalign/template"

if [ -d "$TEMPLATE_SRC" ]; then
    # Copy interfaces if missing
    INT_SRC="$TEMPLATE_SRC/interfaces"
    INT_DST="$TEMPLATE_DST/interfaces"
    if [ -d "$INT_SRC" ]; then
        echo "Copying template interfaces..."
        find "$INT_SRC" -maxdepth 1 -name '*.int' -type f -exec cp -n -t "$INT_DST" {} +
    fi
    
    # Copy contact files if missing
    CONT_SRC="$TEMPLATE_SRC/contact"
    CONT_DST="$TEMPLATE_DST/contact"
    if [ -d "$CONT_SRC" ]; then
        echo "Copying template contacts..."
        find "$CONT_SRC" -maxdepth 1 -name '*.txt' -type f -exec cp -n -t "$CONT_DST" {} +
    fi
    
    # Copy hotspot files if missing
    HOT_SRC="$TEMPLATE_SRC/hotspot"
    HOT_DST="$TEMPLATE_DST/hotspot"
    if [ -d "$HOT_SRC" ]; then
        echo "Copying template hotspots..."
        find "$HOT_SRC" -maxdepth 1 -type f ! -name '.*' -exec cp -n -t "$HOT_DST" {} +
    fi
fi
echo "Template interface count:" $(ls "$TEMPLATE_DST/interfaces/" 2>/dev/null | wc -l)
echo "Template contact count:" $(ls "$TEMPLATE_DST/contact/" 2>/dev/null | wc -l)
echo "Template hotspot count:" $(ls "$TEMPLATE_DST/hotspot/" 2>/dev/null | wc -l)
echo ""

# Setup environment
eval "$(conda shell.bash hook)"
conda activate tmalignRosetta
module load rosetta/2022.42

cd /scratch/rshadi25/GitHub/PRISM-prescript/working_version/TMalign

echo "=== Running TMalign pipeline with pairlist_compare, template_default, job: tmalign-compare ==="
date

python2 prism.py pairlist_compare run_files/template_default tmalign-compare

exit_code=$?
echo "=== Pipeline finished ==="
echo "Exit code: $exit_code"
date

# === REPORT ===
echo ""
echo "============= RESULTS REPORT ============="
echo "Job: tmalign-compare"
echo ""

# Stage 1-2: Download + PreProcess
echo "--- Stage 1-2: PDB Download + PreProcess ---"
echo "PDB downloaded:"
ls jobs/tmalign-compare/pdb/ 2>/dev/null
echo ""
echo "Preprocessed files:"
ls jobs/tmalign-compare/preprocess/ 2>/dev/null
echo ""

# Stage 3: Surface Extraction
echo "--- Stage 3: Surface Extraction ---"
echo "ASA files:"
ls jobs/tmalign-compare/surfaceExtract/*.asa.pdb 2>/dev/null
echo "RSA files:"
ls jobs/tmalign-compare/surfaceExtract/*.rsa 2>/dev/null
echo ""

# Stage 4: Structural Alignment (TMalign)
echo "--- Stage 4: TMalign Alignment ---"
echo "Alignment directory:"
ls jobs/tmalign-compare/alignment/*_*_* 2>/dev/null | wc -l
echo " total alignment files"
echo ""

# Check alignment quality
echo "Sample alignment contents:"
python2 -c "
import pickle, glob, os
files = sorted(glob.glob('jobs/tmalign-compare/alignment/*_*_*'))
total = len(files)
valid = sum(1 for f in files if pickle.load(open(f)) != -1)
empty = sum(1 for f in files if pickle.load(open(f)) == -1)
print('  Total alignment pickles: %d' % total)
print('  Valid (found solutions): %d' % valid)
print('  Empty (no solution): %d' % empty)
if valid > 0:
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
ls jobs/tmalign-compare/transformation/ 2>/dev/null
echo ""
echo "Passed files content:"
cat jobs/tmalign-compare/transformation/passedFiles 2>/dev/null || echo "  (empty or missing)"
echo ""

# Stage 6: Rosetta Refinement
echo "--- Stage 6: Rosetta Refinement ---"
echo "FlexibleRefinement directory:"
ls jobs/tmalign-compare/flexibleRefinement/ 2>/dev/null
echo ""
echo "Energies:"
cat jobs/tmalign-compare/flexibleRefinement/energies/refinement_energies 2>/dev/null || echo "  (no energies produced)"
echo ""
echo "Structures produced:"
ls jobs/tmalign-compare/flexibleRefinement/structures/ 2>/dev/null | head -10 || echo "  (none)"
echo ""

# Check rosetta output
echo "Rosetta output folder:"
ls rosetta_output_1/tmalign-compare/ 2>/dev/null | head -20 || echo "  (empty)"
echo ""

echo "==========================================="
exit $exit_code
