#!/bin/bash
# Run the derived legacy PRISM controller under two explicitly different
# NACCESS profiles. Each task has its own workspace and output tree.
#SBATCH --job-name=array-kutem
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --array=1-2

set -euo pipefail

REPO_ROOT="${MP_REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
OUTPUT_ROOT="${MP_PIPELINE_SMOKE_ROOT:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/legacy-pipeline-smoke-v10}"
PYTHON3="${MP_PYTHON3:-/opt/ohpc/pub/compiler/conda3/latest/bin/python3}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
CURRENT_ENV_ROOT="${MP_CURRENT_ENV_ROOT:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/legacy-tool-environment-v5}"
CURRENT_WORKSPACE="${MP_CURRENT_WORKSPACE:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/legacy-pipeline-workspace-current-v8}"
HISTORICAL_ENV_ROOT="${MP_HISTORICAL_ENV_ROOT:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/legacy-tool-environment-historical-v5}"
HISTORICAL_WORKSPACE="${MP_HISTORICAL_WORKSPACE:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/legacy-pipeline-workspace-historical-v9}"

case "$TASK_ID" in
    1)
        PROFILE="compatibility_naccess"
        ENV_ROOT="$CURRENT_ENV_ROOT"
        WORKSPACE="$CURRENT_WORKSPACE"
        ;;
    2)
        PROFILE="historical_naccess"
        ENV_ROOT="$HISTORICAL_ENV_ROOT"
        WORKSPACE="$HISTORICAL_WORKSPACE"
        ;;
    *)
        PROFILE="unsupported"
        ENV_ROOT="$REPO_ROOT/.unsupported-environment"
        WORKSPACE="$REPO_ROOT/.unsupported-workspace"
        ;;
esac

SOURCE_WORKSPACE="$WORKSPACE"
JOB_ID="${SLURM_JOB_ID:-local-${BASHPID}}"
REQUEUE_COUNT="${SLURM_RESTART_COUNT:-0}"
RETRY_ID="${JOB_ID}-r${REQUEUE_COUNT}"
TASK_ROOT="$OUTPUT_ROOT/task-${TASK_ID}/retry-${RETRY_ID}"
EXIT_WRITTEN=0
export TASK_ROOT TASK_ID RETRY_ID PROFILE SOURCE_WORKSPACE ENV_ROOT JOB_ID REQUEUE_COUNT
write_failure_exit() {
    local status="$1"
    if [[ -x "$PYTHON3" ]]; then
        FAILURE_STATUS="$status" "$PYTHON3" - <<'PY' || true
import json
import os
from pathlib import Path
root = Path(os.environ["TASK_ROOT"])
value = {
    "job_id": os.environ.get("SLURM_JOB_ID", ""),
    "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
    "array_task_id": os.environ.get("TASK_ID", ""),
    "retry_id": os.environ.get("RETRY_ID", ""),
    "requeue_count": int(os.environ.get("REQUEUE_COUNT", "0")),
    "profile": os.environ.get("PROFILE", ""),
    "environment_root": os.environ.get("ENV_ROOT", ""),
    "source_workspace": os.environ.get("SOURCE_WORKSPACE", ""),
    "return_code": int(os.environ.get("FAILURE_STATUS", "1")),
    "scientific_status": "launcher_failed_before_controller",
}
(root / "exit.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
PY
    else
        printf '{"return_code":%s,"scientific_status":"launcher_failed_before_controller"}\n' "$status" > "$TASK_ROOT/exit.json" || true
    fi
}
on_exit() {
    local status=$?
    if [[ "$EXIT_WRITTEN" -eq 0 ]]; then
        write_failure_exit "$status"
    fi
    trap - EXIT
    exit "$status"
}
trap on_exit EXIT
mkdir -p "$OUTPUT_ROOT/task-${TASK_ID}"
if [[ -e "$TASK_ROOT" ]] && [[ -n "$(find "$TASK_ROOT" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "refusing to reuse non-empty retry directory: $TASK_ROOT" >&2
    TASK_ROOT="$OUTPUT_ROOT/task-${TASK_ID}/launcher-failure-${RETRY_ID}-${BASHPID}"
    mkdir -p "$TASK_ROOT"
    export TASK_ROOT
    exit 3
fi
mkdir -p "$TASK_ROOT"
test -f "$SOURCE_WORKSPACE/smoke_manifest.json"
test -f "$ENV_ROOT/environment_manifest.json"
test -x "$ENV_ROOT/bin/python2"
test -x "$PYTHON3"

cp -a "$SOURCE_WORKSPACE" "$TASK_ROOT/workspace"
WORKSPACE="$TASK_ROOT/workspace"
cp "$WORKSPACE/smoke_manifest.json" "$TASK_ROOT/smoke_manifest.json"
if [[ -f "$WORKSPACE/distinct_chain_probe_manifest.json" ]]; then
    cp "$WORKSPACE/distinct_chain_probe_manifest.json" "$TASK_ROOT/distinct_chain_probe_manifest.json"
fi
# The historical controller calls os.mkdir on the final job directory but
# does not create its parent.  Initialize that expected parent here; this is
# launcher plumbing and does not change FiberDock inputs or algorithms.
mkdir -p "$WORKSPACE/jobs"
mkdir -p "$WORKSPACE/fiberdock_output/smoke"
SMOKE_MATCH_PERCENTAGE="${MP_SMOKE_MINIMUM_RESIDUE_MATCH_PERCENTAGE:-}"
if [[ -n "$SMOKE_MATCH_PERCENTAGE" ]]; then
    if [[ ! "$SMOKE_MATCH_PERCENTAGE" =~ ^(0\.[0-9]+|[1-9][0-9]?(\.[0-9]+)?|100(\.0*)?)$ ]]; then
        echo "invalid diagnostic minimum-residue-match percentage: $SMOKE_MATCH_PERCENTAGE" >&2
        exit 2
    fi
    if [[ "$(grep -Ec '^minimum_residue_match_percentage[[:space:]]*=' "$WORKSPACE/prism.ini")" -ne 1 ]]; then
        echo "expected exactly one minimum_residue_match_percentage setting" >&2
        exit 2
    fi
    sed -i -E "s/^minimum_residue_match_percentage[[:space:]]*=.*/minimum_residue_match_percentage = ${SMOKE_MATCH_PERCENTAGE}/" "$WORKSPACE/prism.ini"
fi
PAIR_SPEC=$(tr '\n' ' ' < "$WORKSPACE/input/pair_list" | sed 's/[[:space:]]*$//')
TEMPLATE_SPEC=$(tr '\n' ' ' < "$WORKSPACE/input/template_list" | sed 's/[[:space:]]*$//')
printf 'profile\tenvironment_root\tsource_workspace\ttask_workspace\tpair\ttemplate\tminimum_residue_match_percentage\n' > "$TASK_ROOT/parameters.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$PROFILE" "$ENV_ROOT" "$SOURCE_WORKSPACE" "$WORKSPACE" "$PAIR_SPEC" "$TEMPLATE_SPEC" "${SMOKE_MATCH_PERCENTAGE:-from-source-workspace}" >> "$TASK_ROOT/parameters.tsv"
printf '%s\n' "source $ENV_ROOT/activate.sh && cd $WORKSPACE && $ENV_ROOT/bin/python2 prism.py input/pair_list input/template_list smoke" > "$TASK_ROOT/command.txt"

START_EPOCH=$(date +%s)
export PRISM_COMPAT_EVENT_LOG="$TASK_ROOT/compat_events.jsonl"
source "$ENV_ROOT/activate.sh"
cd "$WORKSPACE"
set +e
"$ENV_ROOT/bin/python2" prism.py input/pair_list input/template_list smoke > "$TASK_ROOT/stdout.log" 2> "$TASK_ROOT/stderr.log"
RC=$?
set -e
END_EPOCH=$(date +%s)

export TASK_ROOT TASK_ID START_EPOCH END_EPOCH RC PROFILE ENV_ROOT WORKSPACE
set +e
"$PYTHON3" - <<'PY'
import hashlib
import json
import os
import sys
from pathlib import Path

root = Path(os.environ["TASK_ROOT"])
workspace = Path(os.environ["WORKSPACE"])
job = workspace / "jobs" / "smoke"
environment_manifest = json.loads((Path(os.environ["ENV_ROOT"]) / "environment_manifest.json").read_text())
capabilities = environment_manifest.get("capabilities", {})
full_refinement_ready = bool(capabilities.get("fiberdock_full_refinement", False))

def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

outputs = {}
for artifact_root in (job, workspace / "fiberdock_output"):
    if artifact_root.exists():
        for path in sorted(artifact_root.rglob("*")):
            if path.is_file():
                outputs[str(path.relative_to(workspace))] = {"size": path.stat().st_size, "sha256": sha(path)}
for path in (root / "stdout.log", root / "stderr.log", root / "compat_events.jsonl"):
    if path.is_file():
        outputs[str(path.name)] = {"size": path.stat().st_size, "sha256": sha(path)}

provenance = {}
for name in ("smoke_manifest.json", "distinct_chain_probe_manifest.json", "parameters.tsv", "command.txt", "input/pair_list", "input/template_list", "workspace/prism.ini", "workspace/config.inc"):
    if name.startswith("input/"):
        path = workspace / name
    elif name.startswith("workspace/"):
        path = workspace / name.split("/", 1)[1]
    else:
        path = root / name
    provenance[name] = {"path": str(path), "sha256": sha(path)} if path.is_file() else {"path": str(path), "sha256": ""}
for name, path in (
    ("environment_manifest.json", Path(os.environ["ENV_ROOT"]) / "environment_manifest.json"),
    ("compatibility_manifest.json", workspace / "compatibility_manifest.json"),
):
    provenance[name] = {"path": str(path), "sha256": sha(path)} if path.is_file() else {"path": str(path), "sha256": ""}

stage_dirs = {name: (job / name).is_dir() for name in ("preprocess", "surfaceExtract", "alignment", "transformation", "fiberdock")}
passed_files = job / "transformation" / "passedFiles"
has_candidates = passed_files.is_file() and passed_files.stat().st_size > 0
has_fiberdock_intermediate = any(path.suffix in {".fib", ".ref", ".ref.pdb"} for path in (job / "fiberdock").rglob("*") if path.is_file()) if (job / "fiberdock").exists() else False
stdout_text = (root / "stdout.log").read_text(errors="replace") if (root / "stdout.log").is_file() else ""
refinement_started = "Flexible Refinement stage started" in stdout_text
capability_boundary_observed = has_candidates and not full_refinement_ready and refinement_started
final_output = workspace / "fiberdock_output" / "smoke"
has_fiberdock_final = bool(list(final_output.glob("*.fiberdock.pdb"))) and bool(list(final_output.glob("*.intRes.txt"))) if final_output.exists() else False
surface_outputs = list((job / "surfaceExtract").glob("*.rsa")) + list((job / "surfaceExtract").glob("*.asa.pdb")) if (job / "surfaceExtract").exists() else []
alignment_outputs = list((job / "alignment").iterdir()) if (job / "alignment").exists() else []
has_required_intermediates = bool(surface_outputs) and bool(alignment_outputs) and passed_files.is_file()
if all(stage_dirs.values()) and capability_boundary_observed:
    status = "pipeline_blocked_full_refinement_capability"
elif int(os.environ["RC"]) != 0:
    status = "controller_failed"
elif all(stage_dirs.values()) and has_candidates and has_fiberdock_intermediate and has_fiberdock_final:
    status = "pipeline_smoke_with_refinement"
elif all(stage_dirs.values()) and has_required_intermediates:
    status = "pipeline_plumbing_complete_no_candidates"
elif all(stage_dirs.values()):
    status = "pipeline_blocked_missing_intermediate"
elif any(stage_dirs.values()):
    status = "pipeline_smoke_partial"
else:
    status = "pipeline_smoke_no_intermediates"
value = {
    "job_id": os.environ.get("SLURM_JOB_ID", ""),
    "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
    "array_task_id": os.environ["TASK_ID"],
    "retry_id": os.environ.get("RETRY_ID", ""),
    "requeue_count": int(os.environ.get("REQUEUE_COUNT", "0")),
    "profile": os.environ["PROFILE"],
    "environment_root": os.environ["ENV_ROOT"],
    "workspace": str(workspace),
    "start_epoch": int(os.environ["START_EPOCH"]),
    "end_epoch": int(os.environ["END_EPOCH"]),
    "return_code": 0 if status in {"pipeline_smoke_with_refinement", "pipeline_plumbing_complete_no_candidates"} else 1,
    "controller_return_code": int(os.environ["RC"]),
    "stage_directories": stage_dirs,
    "outputs": outputs,
    "provenance_files": provenance,
    "requested_resources": {"nodes": "1", "ntasks": "1", "cpus_per_task": "2", "mem": "2G", "time": "00:05:00"},
    "scientific_status": status,
    "candidate_count": sum(1 for _ in passed_files.open()) if passed_files.is_file() else 0,
    "full_refinement_ready": full_refinement_ready,
    "refinement_started": refinement_started,
    "scientific_success": status == "pipeline_smoke_with_refinement",
    "full_refinement_blockers": capabilities.get("fiberdock_full_refinement_blockers", []),
    "refinement_contract": {
        "has_fiberdock_intermediate": has_fiberdock_intermediate,
        "has_fiberdock_final": has_fiberdock_final,
        "final_output_root": str(final_output),
    },
}
(root / "exit.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
raise SystemExit(0 if status in {"pipeline_smoke_with_refinement", "pipeline_plumbing_complete_no_candidates"} else 1)
PY
REPORT_RC=$?
set -e
if [[ -s "$TASK_ROOT/exit.json" ]]; then
    EXIT_WRITTEN=1
fi
exit "$REPORT_RC"
