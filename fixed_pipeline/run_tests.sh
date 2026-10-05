#!/usr/bin/env bash
# Run all fixed-pipeline tests using proper package imports.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FIXED_DIR="$REPO_ROOT/fixed_pipeline"

# Pick a Python interpreter with numpy
PYTHON="${PRISM_TEST_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"

echo "=== Fixed Pipeline Tests ==="
echo "Python: $PYTHON"
echo "Date:   $(date)"
echo ""

# PYTHONPATH needs the parent of src/ (i.e. fixed_pipeline/) so that
# `from src.alignment import ...` resolves the package correctly.
export PYTHONPATH="$FIXED_DIR:$PYTHONPATH"

# Change to the fixed_pipeline directory so relative imports work
cd "$FIXED_DIR"
"$PYTHON" -m pytest tests/ -v --tb=short 2>&1

EXIT_CODE=$?
echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ All fixed-pipeline tests passed."
else
    echo "❌ Some tests failed (exit code $EXIT_CODE)."
fi
exit $EXIT_CODE
