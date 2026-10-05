#!/usr/bin/env bash
# Run Python 3.11 with a run-scoped PyRosetta package bundle.
# This deliberately does not create or modify a Conda environment.
set -euo pipefail

bundle_root="${PRISM_PYROSETTA_BUNDLE:?set PRISM_PYROSETTA_BUNDLE to the run-scoped bundle}"
base_python="${PRISM_PYROSETTA_BASE_PYTHON:-/scratch/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python3.11}"
site_packages="$bundle_root/lib/python3.11/site-packages"
[[ -x "$base_python" ]] || { echo "missing base Python: $base_python" >&2; exit 2; }
[[ -f "$site_packages/pyrosetta/rosetta.so" ]] || { echo "missing PyRosetta bundle: $site_packages" >&2; exit 2; }

export PYTHONNOUSERSITE=1
export PYTHONPATH="$site_packages${PYTHONPATH:+:$PYTHONPATH}"
export PYROSETTA_DATABASE="$site_packages/pyrosetta/database"
base_lib="$(cd -- "$(dirname -- "$base_python")/../lib" && pwd)"
export LD_LIBRARY_PATH="$base_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

exec "$base_python" "$@"
