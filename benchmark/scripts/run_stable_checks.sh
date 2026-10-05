#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
default_python=${PRISM_TEST_PYTHON:-$repo_dir/benchmark/prism_processed/env/prism_score_env/bin/python}
python_bin=${1:-$default_python}

if [[ ! -x "$python_bin" ]]; then
  echo "Python interpreter not found or not executable: $python_bin" >&2
  exit 1
fi

echo "[1/3] Benchmark preflight"
python3 "$repo_dir/benchmark/scripts/preflight_benchmark_check.py" --repo-root "$repo_dir"

echo
echo "[2/3] Repo unittest suite"
"$python_bin" -m unittest discover -s "$repo_dir/tests" -p 'test_*.py'

echo
echo "[3/3] PRISM pipeline helper checks"
"$python_bin" "$repo_dir/benchmark/scripts/test_prism_pipeline_helpers.py"
