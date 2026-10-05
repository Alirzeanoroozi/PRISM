#!/usr/bin/env bash
# Deploy fixed pipeline modules by copying them over originals.
#
# WARNING: This overwrites files in src/ with fixed versions.
# Backup originals first if you want to preserve them.
#
# Usage:
#   ./deploy_fixes.sh              # Copy all fixed modules
#   ./deploy_fixes.sh --dry-run    # Show what would be copied
#   ./deploy_fixes.sh --backup     # Backup originals first
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FIXED_DIR="$REPO_ROOT/fixed_pipeline"
SRC_DIR="$REPO_ROOT/src"

FIXED_FILES=(
    "rosetta_refinement.py"
    "alignment.py"
    "alignment_gtalign.py"
    "alignment_multiprot.py"
    "transformation.py"
    "naccess_utils.py"
    "hotspot.py"
)

DRY_RUN=false
BACKUP=false

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
        --backup)  BACKUP=true ;;
    esac
done

echo "=== Deploy Fixed Pipeline Modules ==="
echo "Source: $FIXED_DIR/src/"
echo "Target: $SRC_DIR/"
echo ""

for fname in "${FIXED_FILES[@]}"; do
    src="$FIXED_DIR/src/$fname"
    dst="$SRC_DIR/$fname"

    if [ ! -f "$src" ]; then
        echo "⚠️   Skipping $fname (not found in $FIXED_DIR/src/)"
        continue
    fi

    if $DRY_RUN; then
        echo "🔍  Would copy: $fname"
        continue
    fi

    if $BACKUP; then
        backup="$dst.bak.$(date +%Y%m%d%H%M%S)"
        cp "$dst" "$backup"
        echo "💾  Backup: $backup"
    fi

    cp "$src" "$dst"
    echo "✅  Deployed: $fname"
done

echo ""
if $DRY_RUN; then
    echo "Dry run complete. Use without --dry-run to deploy."
else
    echo "Deployment complete. Run tests with: bash fixed_pipeline/run_tests.sh"
fi
