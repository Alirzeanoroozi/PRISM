#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "usage: $0 <masif_root> <pdb_path> <chain_id> [structure_id]"
  exit 1
fi

MASIF_ROOT="$1"
PDB_PATH="$2"
CHAIN_ID="$3"
STRUCTURE_ID="${4:-$(basename "$PDB_PATH" .pdb)$CHAIN_ID}"

OUTPUT_DIR="tests/diffmasif_runtime/output/${STRUCTURE_ID}"
RUNTIME_REPORT="${OUTPUT_DIR}/runtime_check.json"

mkdir -p "$OUTPUT_DIR"

python3 tests/diffmasif_runtime/check_runtime.py \
  --masif-root "$MASIF_ROOT" \
  --pdb "$PDB_PATH" \
  --output "$RUNTIME_REPORT"

cat <<EOF

Runtime report written to:
  $RUNTIME_REPORT

Suggested next manual step in the external MaSIF checkout:

  cd "$MASIF_ROOT"
  # activate the MaSIF / DiffMaSIF environment first
  # run the upstream preprocessing/surface generation workflow for:
  #   pdb: $PDB_PATH
  #   chain: $CHAIN_ID

After you export surface points as CSV, compare them against Naccess with:

  python3 tests/diffmasif_replacement/compare_with_naccess.py \\
    --structure-kind target \\
    --structure-id "$STRUCTURE_ID" \\
    --surface-points /path/to/diffmasif_points.csv \\
    --assignment-cutoff 4.0 \\
    --output tests/diffmasif_replacement/output/${STRUCTURE_ID}.diffmasif.json

EOF
