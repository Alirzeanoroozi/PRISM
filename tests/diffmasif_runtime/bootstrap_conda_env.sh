#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-masif_py36}"

echo "Creating conda env: ${ENV_NAME}"
conda create -y -n "${ENV_NAME}" python=3.6 pip

echo "Installing baseline Python packages into ${ENV_NAME}"
conda run -n "${ENV_NAME}" pip install \
  biopython==1.79 \
  colour \
  ipython

cat <<EOF

Created env: ${ENV_NAME}

Next steps:
  1. Activate it:
       conda activate ${ENV_NAME}
  2. Install any additional MaSIF-compatible packages you can obtain for this machine:
       pip install open3d
  3. Install external binaries separately:
       reduce, MSMS, PDB2PQR, multivalue, APBS
  4. Fill in:
       tests/diffmasif_runtime/masif_env_template.sh
  5. Re-run:
       python3 tests/diffmasif_runtime/check_runtime.py --masif-root /scratch/rshadi25/GitHub/masif --pdb processed/pdbs/1fgn.pdb

Notes:
  - PyMesh is the hardest dependency and may require a separate build strategy.
  - The upstream MaSIF repo cloned here does not contain DiffMaSIF itself.

EOF
