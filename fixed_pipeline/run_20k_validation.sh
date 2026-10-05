#!/usr/bin/env bash
# Run PRISM pipelines with the full 19,855-template panel and validate results.
#
# This script:
# 1. Creates isolated run directories for each variant
# 2. Links all required assets (templates, PDBs, external tools)
# 3. Runs each variant with a 1-hour timeout
# 4. Validates alignment JSONs, transform PDBs, and refinement outputs
# 5. Reports a summary table
#
# Usage: bash fixed_pipeline/run_20k_validation.sh [--quick] [--tm-only]
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RUN_ROOT="${PRISM_RUN_ROOT:-$REPO_ROOT/tmp/agent/$(date +%Y%m%d)-pipeline-20k}"
TIMEOUT_SEC="${PRISM_TIMEOUT_SEC:-3600}"  # 1 hour default
PYTHON="${PRISM_PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
TEMPLATE_LIST="$REPO_ROOT/new_template/template/final_list.txt"
TEMPLATE_COUNT=$(wc -l < "$TEMPLATE_LIST")

# Parse args
QUICK=false
TM_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --quick) QUICK=true ;;
        --tm-only) TM_ONLY=true ;;
    esac
done

echo "========================================================"
echo "  PRISM 20K Template Pipeline Validation"
echo "========================================================"
echo "Repository: $REPO_ROOT"
echo "Python:     $PYTHON"
echo "Templates:  $TEMPLATE_COUNT (from $TEMPLATE_LIST)"
echo "Run root:   $RUN_ROOT"
echo "Timeout:    ${TIMEOUT_SEC}s"
echo "Quick mode: $QUICK"
echo ""

# Limit templates in quick mode
if $QUICK; then
    TEMPLATE_LIMIT=100
    echo "⚠️  QUICK MODE: Using $TEMPLATE_LIMIT templates instead of $TEMPLATE_COUNT"
else
    TEMPLATE_LIMIT="$TEMPLATE_COUNT"
fi

setup_run() {
    local variant="$1"
    local run_dir="$RUN_ROOT/$variant"
    mkdir -p "$run_dir"/{processed/pdbs,processed/surface_extraction,processed/alignment,processed/transformation,processed/rosetta_refinement,templates,status}
    
    # Symlink assets
    ln -sfn "$REPO_ROOT/src" "$run_dir/src"
    ln -sfn "$REPO_ROOT/external_tools" "$run_dir/external_tools"
    ln -sfn "$REPO_ROOT/prism.py" "$run_dir/prism.py"
    ln -sfn "$REPO_ROOT/processed/pdbs" "$run_dir/processed/pdbs"
    ln -sfn "$REPO_ROOT/new_template/template/interfaces" "$run_dir/templates/interfaces"
    ln -sfn "$REPO_ROOT/new_template/template/interfaces_lists" "$run_dir/templates/interfaces_lists"
    ln -sfn "$REPO_ROOT/new_template/template/hotspots" "$run_dir/templates/hotspots" 2>/dev/null || true
    ln -sfn "$REPO_ROOT/new_template/template/contacts" "$run_dir/templates/contacts" 2>/dev/null || true
    ln -sfn "$REPO_ROOT/new_template/template/rsas" "$run_dir/templates/rsas" 2>/dev/null || true
    ln -sfn "$REPO_ROOT/new_template/template/pdbs" "$run_dir/templates/pdbs" 2>/dev/null || true
    
    # Write inputs.csv with a single pair
    cat > "$run_dir/inputs.csv" << EOF
Receptor,Ligand
1fgnHL,1tfhA
EOF
    
    # Write template list
    head -n "$TEMPLATE_LIMIT" "$TEMPLATE_LIST" > "$run_dir/templates/calculated_templates.txt"
    
    echo "$run_dir"
}

run_variant() {
    local variant="$1"
    local aligner="$2"
    local refiner="$3"
    local extra_args="${4:-}"
    local run_dir=$(setup_run "$variant")
    
    echo ""
    echo "--- Running $variant ($aligner + $refiner) ---"
    echo "  Dir: $run_dir"
    
    local start_time=$(date +%s)
    
    (
        cd "$run_dir"
        export PRISM_INPUTS_CSV="inputs.csv"
        export PRISM_TM_SCORE_THRESHOLD="0.4"
        export PRISM_FILTER_MODE="geometry_only_experimental"
        export PRISM_STAGE_STATUS_PATH="status/stages.jsonl"
        export PRISM_SURFACE_BACKEND="freesasa"
        export PRISM_FREESASA_PYTHON="$PYTHON"
        
        timeout "$TIMEOUT_SEC" "$PYTHON" prism.py \
            --aligner "$aligner" \
            --refiner "$refiner" \
            --surface_backend freesasa \
            $extra_args \
            2>&1
    )
    local exit_code=$?
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    
    echo "  Exit code: $exit_code (${elapsed}s)"
    
    # Collect results
    local align_count=$(/usr/bin/find "$run_dir/processed/alignment" -maxdepth 1 -name '*.json' 2>/dev/null | wc -l)
    local trans_count=$(/usr/bin/find "$run_dir/processed/transformation" -maxdepth 1 -name '*.pdb' 2>/dev/null | wc -l)
    local refine_count=$(/usr/bin/find "$run_dir/processed/rosetta_refinement/structures" -maxdepth 1 -name '*.pdb' 2>/dev/null | wc -l)
    local passed_pairs=0
    if [ -f "$run_dir/processed/rosetta_refinement/refinement_energies.txt" ]; then
        passed_pairs=$(grep -c '^' "$run_dir/processed/rosetta_refinement/refinement_energies.txt" 2>/dev/null || echo 0)
    fi
    
    # Validate alignment JSONs
    local valid_align=0
    local empty_align=0
    for f in $(/usr/bin/find "$run_dir/processed/alignment" -maxdepth 1 -name '*.json' 2>/dev/null | head -1000); do
        local mc=$(python3 -c "import json; d=json.load(open('$f')); print(d.get('match_count', -1))" 2>/dev/null || echo "-2")
        if [ "$mc" -gt 0 ] 2>/dev/null; then
            valid_align=$((valid_align + 1))
        else
            empty_align=$((empty_align + 1))
        fi
    done
    
    echo "  Alignment JSONs: $align_count total ($valid_align valid, $empty_align empty in first 1000)"
    echo "  Transform PDBs:  $trans_count"
    echo "  Refined PDBs:    $refine_count"
    echo "  Passed pairs:    $passed_pairs"
    
    # Store summary
    cat >> "$RUN_ROOT/summary.tsv" << EOF
$variant	$aligner	$refiner	$exit_code	${elapsed}s	$align_count	$valid_align	$empty_align	$trans_count	$refine_count	$passed_pairs
EOF
}

# Header for summary
mkdir -p "$RUN_ROOT"
echo -e "variant\taligner\trefiner\texit_code\telapsed\talign_total\talign_valid\talign_empty\ttransform_pdb\trefine_pdb\tpassed_pairs" > "$RUN_ROOT/summary.tsv"

# Run variants
if $TM_ONLY; then
    run_variant "tm_external" "tmalign" "external_rosetta"
    run_variant "tm_pyro" "tmalign" "pyrosetta"
else
    run_variant "tm_external" "tmalign" "external_rosetta"
    run_variant "tm_pyro" "tmalign" "pyrosetta"
    
    if ! $QUICK; then
        run_variant "multiprot_external" "multiprot" "external_rosetta" "--template-limit 100"
        run_variant "multiprot_pyro" "multiprot" "pyrosetta" "--template-limit 100"
    else
        # MultiProt with small template set even in quick mode
        run_variant "multiprot_external" "multiprot" "external_rosetta" "--template-limit 100"
        run_variant "multiprot_pyro" "multiprot" "pyrosetta" "--template-limit 100"
    fi
fi

echo ""
echo "========================================================"
echo "  RESULTS SUMMARY"
echo "========================================================"
column -t "$RUN_ROOT/summary.tsv"
echo ""
echo "Full results: $RUN_ROOT/summary.tsv"
echo "========================================================"

# Check for failures
FAILED=$(grep -c '^' "$RUN_ROOT/summary.tsv" 2>/dev/null || true)
if [ "$FAILED" -gt 0 ]; then
    echo "⚠️  Some variants had non-zero exit codes. Check individual run logs."
fi
