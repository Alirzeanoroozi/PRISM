#!/usr/bin/env bash
# Install an authorized local PyRosetta wheel into a separate venv.
set -euo pipefail

wheel="${1:?usage: install_pyrosetta_prism.sh /path/to/authorized-pyrosetta.whl [sha256]}"
expected_sha="${2:-}"
env_dir="${PYROSETTA_ENV:-/home/rshadi25/.conda/envs/pyrosetta_prism}"
python_base="${PYTHON_BASE:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"

if [[ ! -f "$wheel" ]]; then
  echo "wheel not found: $wheel" >&2
  exit 2
fi
if [[ "$(basename "$wheel")" != *cp311-cp311-linux_x86_64.whl ]]; then
  echo "wheel is not compatible with the pinned Python 3.11 Linux x86_64 environment: $wheel" >&2
  exit 2
fi
actual_sha=$(sha256sum "$wheel" | awk '{print $1}')
if [[ -n "$expected_sha" && "$actual_sha" != "$expected_sha" ]]; then
  echo "wheel SHA256 mismatch: expected $expected_sha observed $actual_sha" >&2
  exit 3
fi
if [[ -e "$env_dir" ]]; then
  echo "refusing to overwrite existing environment: $env_dir" >&2
  exit 4
fi
"$python_base" -m venv --system-site-packages "$env_dir"
"$env_dir/bin/python" -m pip install --no-index --no-deps "$wheel"
"$env_dir/bin/python" - <<'PY'
import json
import pyrosetta
print(json.dumps({"version": getattr(pyrosetta, "__version__", None), "module": pyrosetta.__file__}))
PY
printf 'wheel=%s\nsha256=%s\nenvironment=%s\npython=%s\n' "$wheel" "$actual_sha" "$env_dir" "$env_dir/bin/python"
