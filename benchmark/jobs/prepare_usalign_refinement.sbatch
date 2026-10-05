#!/usr/bin/env bash
# Prepare the USalign common-refinement manifest only after all compact batch
# outputs validate.  This stage is read-only with respect to transformed PDBs.
# It does not submit the refinement array; the selected count is inspected
# before that separate resource placement decision.
#SBATCH --job-name=prism-usalign-refine-prep
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

: "${REPO_ROOT:?set REPO_ROOT}"
: "${RUN_ROOT:?set RUN_ROOT}"
: "${PIPELINE_PYTHON:?set PIPELINE_PYTHON}"

aggregate_dir="${RUN_ROOT}/aggregate"
handoff_dir="${RUN_ROOT}/usalign_common_refinement_19855"
mkdir -p "${handoff_dir}" "${RUN_ROOT}/logs"

"${PIPELINE_PYTHON}" - "${aggregate_dir}/aggregation_status.json" "${handoff_dir}/handoff_status.json" <<'PY'
import json
import sys
from pathlib import Path

status_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
payload = json.loads(status_path.read_text(encoding="utf-8")) if status_path.is_file() else {"status": "missing"}
if payload.get("status") != "validated_compacted":
    output_path.write_text(json.dumps({
        "schema_version": "prism-usalign-refinement-handoff/v1",
        "status": "blocked_upstream_aggregation",
        "aggregation_status": payload,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    raise SystemExit(2)
PY

"${PIPELINE_PYTHON}" \
  "${REPO_ROOT}/benchmark/scripts/prepare_usalign_refinement_manifest.py" \
  --source "${aggregate_dir}/transformed_dockq.tsv" \
  --output "${handoff_dir}/selected_candidates.csv" \
  --rejected "${handoff_dir}/rejected_candidates.csv" \
  --manifest "${handoff_dir}/manifest.json" \
  > "${RUN_ROOT}/logs/prepare_usalign_refinement.out" \
  2> "${RUN_ROOT}/logs/prepare_usalign_refinement.err"

"${PIPELINE_PYTHON}" - "${handoff_dir}/manifest.json" "${handoff_dir}/handoff_status.json" <<'PY'
import json
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
status_path = Path(sys.argv[2])
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
manifest["status"] = "validated_manifest_ready_for_resource_review"
manifest["refinement_array"] = "not_submitted"
status_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
