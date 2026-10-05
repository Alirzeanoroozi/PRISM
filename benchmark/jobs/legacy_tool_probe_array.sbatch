#!/bin/bash
# Probe each staged native tool in an isolated KUTEM array task.
#SBATCH --job-name=array-kutem
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --array=1-4

set -euo pipefail

REPO_ROOT="${MP_REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
ENV_ROOT="${MP_LEGACY_ENV_ROOT:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/legacy-tool-environment-v5}"
OUTPUT_ROOT="${MP_LEGACY_PROBE_ROOT:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/legacy-tool-probes-v5}"
PYTHON3="${MP_PYTHON3:-/opt/ohpc/pub/compiler/conda3/latest/bin/python3}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
JOB_ID="${SLURM_JOB_ID:-local-${BASHPID}}"
REQUEUE_COUNT="${SLURM_RESTART_COUNT:-0}"
RETRY_ID="${JOB_ID}-r${REQUEUE_COUNT}"
TASK_ROOT="$OUTPUT_ROOT/task-${TASK_ID}/retry-${RETRY_ID}"
START_EPOCH=$(date +%s)

EXIT_WRITTEN=0
export TASK_ROOT TASK_ID RETRY_ID ENV_ROOT JOB_ID REQUEUE_COUNT
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
    "environment_root": os.environ.get("ENV_ROOT", ""),
    "return_code": int(os.environ.get("FAILURE_STATUS", "1")),
    "scientific_status": "launcher_failed_before_tool",
}
(root / "exit.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
PY
    else
        printf '{"return_code":%s,"scientific_status":"launcher_failed_before_tool"}\n' "$status" > "$TASK_ROOT/exit.json" || true
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
    echo "refusing to reuse non-empty task directory: $TASK_ROOT" >&2
    TASK_ROOT="$OUTPUT_ROOT/task-${TASK_ID}/launcher-failure-${RETRY_ID}-${BASHPID}"
    mkdir -p "$TASK_ROOT"
    export TASK_ROOT
    exit 3
fi
mkdir -p "$TASK_ROOT"

INPUTS=()
COMMAND=()
TOOL=""
EXPECTED_OUTPUTS=()
case "$TASK_ID" in
    1)
        TOOL="naccess"
        INPUTS=("$REPO_ROOT/benchmark/data/pdbs/rigid/1fgn.pdb")
        COMMAND=("$ENV_ROOT/external_tools/naccess/naccess" "input.pdb" "-f")
        EXPECTED_OUTPUTS=("input.rsa")
        cp "${INPUTS[0]}" "$TASK_ROOT/input.pdb"
        ;;
    2)
        TOOL="pops"
        INPUTS=("$REPO_ROOT/benchmark/data/pdbs/rigid/1fgn.pdb")
        COMMAND=("$ENV_ROOT/external_tools/pops/bin/pops" "--pdb" "input.pdb" "--residueOut" "--popsOut" "pops.rsa")
        EXPECTED_OUTPUTS=("pops.rsa")
        cp "${INPUTS[0]}" "$TASK_ROOT/input.pdb"
        ;;
    3)
        TOOL="multiprot"
        INPUTS=("$REPO_ROOT/benchmark/data/pdbs/rigid/1jpt.pdb" "$REPO_ROOT/benchmark/data/pdbs/rigid/1mkf.pdb")
        COMMAND=("$ENV_ROOT/bin/multiprot" "receptor.pdb" "ligand.pdb")
        EXPECTED_OUTPUTS=("2_sol.res" "2_sets.res")
        cp "${INPUTS[0]}" "$TASK_ROOT/receptor.pdb"
        cp "${INPUTS[1]}" "$TASK_ROOT/ligand.pdb"
        cp "$ENV_ROOT/external_tools/multiprot/params.txt" "$TASK_ROOT/params.txt"
        ;;
    4)
        TOOL="fiberdock"
        EXAMPLE="$ENV_ROOT/external_tools/fiberdock/example"
        INPUTS=("$EXAMPLE/1t6g_r.pdb" "$EXAMPLE/1t6g_l.pdb" "$EXAMPLE/1t6g_trans.txt")
        cp "${INPUTS[0]}" "$TASK_ROOT/1t6g_r.pdb"
        cp "${INPUTS[1]}" "$TASK_ROOT/1t6g_l.pdb"
        cp "${INPUTS[2]}" "$TASK_ROOT/1t6g_trans.txt"
        cp "$EXAMPLE/paramsFile" "$TASK_ROOT/paramsFile"
        sed -i \
            -e 's#^receptorPDBFileName .*#receptorPDBFileName 1t6g_r.pdb#' \
            -e 's#^ligandPDBFileName .*#ligandPDBFileName 1t6g_l.pdb#' \
            -e 's#^transFileName .*#transFileName 1t6g_trans.txt#' \
            -e 's#^rotamerLibFile .*#rotamerLibFile ROTAMER_LIB#' \
            -e 's#^protLib .*#protLib PROT_LIB#' \
            -e 's#^exclusionsFile .*#exclusionsFile EXCLUSIONS#' \
            -e 's#^charmmParamFile .*#charmmParamFile CHARMM_PARAM#' \
            -e 's#^definitionsFile .*#definitionsFile DEFINITIONS#' \
            -e 's#^performBackboneRef .*#performBackboneRef 0#' \
            -e 's#^onlyEnergyCalculation .*#onlyEnergyCalculation 1#' \
            -e 's#^printRefinedComplexes .*#printRefinedComplexes 0#' \
            -e 's#^rigidBodyMinimizationCycles .*#rigidBodyMinimizationCycles 0#' \
            -e 's#^energiesOutFileName .*#energiesOutFileName resFile#' \
            "$TASK_ROOT/paramsFile"
        sed -i \
            -e "s#ROTAMER_LIB#$ENV_ROOT/external_tools/fiberdock/lib/bbdep02.May.sortlib#" \
            -e "s#PROT_LIB#$ENV_ROOT/external_tools/fiberdock/lib/chem.lib#" \
            -e "s#EXCLUSIONS#$ENV_ROOT/external_tools/fiberdock/lib/exclusions#" \
            -e "s#CHARMM_PARAM#$ENV_ROOT/external_tools/fiberdock/lib/par_all27_prot_na.prm#" \
            -e "s#DEFINITIONS#$ENV_ROOT/external_tools/fiberdock/lib/top_all27_prot_na.rtf#" \
            "$TASK_ROOT/paramsFile"
        COMMAND=("$ENV_ROOT/external_tools/fiberdock/FiberDock" "paramsFile")
        EXPECTED_OUTPUTS=("resFile.ref")
        ;;
    *)
        TOOL="unsupported"
        INPUTS=()
        COMMAND=("$REPO_ROOT/.unsupported-tool")
        EXPECTED_OUTPUTS=("unsupported-output")
        ;;
esac

for input in "${INPUTS[@]}"; do test -f "$input"; done
test -f "$ENV_ROOT/environment_manifest.json"
test -x "$PYTHON3"
source "$ENV_ROOT/activate.sh"
"$PYTHON3" - "$ENV_ROOT/environment_manifest.json" "$ENV_ROOT" <<'PY'
import json
import os
import sys
from pathlib import Path

manifest = json.loads(Path(sys.argv[1]).read_text())
if manifest["status"] not in {"ready_compatibility_naccess_toolchain", "ready_historical_legacy_toolchain"}:
    raise SystemExit("unexpected staged environment status: " + manifest["status"])
if not manifest["dependencies"]["numpy"]["available"]:
    raise SystemExit("NumPy dependency gate failed")
if not manifest["dependencies"]["mysqldb_compat"]["available"]:
    raise SystemExit("MySQL compatibility gate failed")
PY

printf 'tool\tarray_task_id\tsource_input\tstaged_input\tsha256\n' > "$TASK_ROOT/input_manifest.tsv"
for index in "${!INPUTS[@]}"; do
    staged=""
    case "$TASK_ID:$index" in
        1:0|2:0) staged="$TASK_ROOT/input.pdb" ;;
        3:0) staged="$TASK_ROOT/receptor.pdb" ;;
        3:1) staged="$TASK_ROOT/ligand.pdb" ;;
        4:0) staged="$TASK_ROOT/1t6g_r.pdb" ;;
        4:1) staged="$TASK_ROOT/1t6g_l.pdb" ;;
        4:2) staged="$TASK_ROOT/1t6g_trans.txt" ;;
    esac
    printf '%s\t%s\t%s\t%s\t%s\n' "$TOOL" "$TASK_ID" "${INPUTS[$index]}" "$staged" "$(sha256sum "$staged" | cut -d ' ' -f1)" >> "$TASK_ROOT/input_manifest.tsv"
done
printf '%s\n' "${COMMAND[*]}" > "$TASK_ROOT/command.txt"
printf 'tool\tarray_task_id\tenvironment\tcommand\n%s\t%s\t%s\t%s\n' "$TOOL" "$TASK_ID" "$ENV_ROOT" "${COMMAND[*]}" > "$TASK_ROOT/parameters.tsv"

cd "$TASK_ROOT"
set +e
"${COMMAND[@]}" > stdout.log 2> stderr.log
RC=$?
set -e
END_EPOCH=$(date +%s)

EXPECTED_OUTPUT_LIST="${EXPECTED_OUTPUTS[*]}"
export TASK_ROOT TASK_ID RETRY_ID START_EPOCH END_EPOCH RC TOOL ENV_ROOT EXPECTED_OUTPUT_LIST
set +e
"$PYTHON3" - <<'PY'
import hashlib
import json
import os
from pathlib import Path

root = Path(os.environ["TASK_ROOT"])
def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

outputs = {}
for path in sorted(root.iterdir()):
    if path.name in {"exit.json", "input_manifest.tsv", "parameters.tsv", "command.txt"} or not path.is_file():
        continue
    outputs[path.name] = {"path": str(path), "size": path.stat().st_size, "sha256": sha(path)}
expected = os.environ.get("EXPECTED_OUTPUT_LIST", "").split()
expected_ok = all((root / name).is_file() and (root / name).stat().st_size > 0 for name in expected)
provenance = {}
for name in ("input_manifest.tsv", "parameters.tsv", "command.txt"):
    path = root / name
    provenance[name] = {"path": str(path), "sha256": sha(path)} if path.is_file() else {"path": str(path), "sha256": ""}
manifest = Path(os.environ["ENV_ROOT"]) / "environment_manifest.json"
provenance["environment_manifest.json"] = {"path": str(manifest), "sha256": sha(manifest)}
value = {
    "job_id": os.environ.get("SLURM_JOB_ID", ""),
    "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
    "array_task_id": os.environ["TASK_ID"],
    "retry_id": os.environ.get("RETRY_ID", ""),
    "requeue_count": int(os.environ.get("REQUEUE_COUNT", "0")),
    "tool": os.environ["TOOL"],
    "environment_root": os.environ["ENV_ROOT"],
    "start_epoch": int(os.environ["START_EPOCH"]),
    "end_epoch": int(os.environ["END_EPOCH"]),
    "return_code": int(os.environ["RC"]),
    "outputs": outputs,
    "expected_outputs": expected,
    "provenance_files": provenance,
    "requested_resources": {"nodes": "1", "ntasks": "1", "cpus_per_task": "2", "mem": "2G", "time": "00:05:00"},
    "scientific_status": "probe_pass" if int(os.environ["RC"]) == 0 and expected_ok else "probe_failed",
}
(root / "exit.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
if value["scientific_status"] != "probe_pass":
    raise SystemExit(1)
PY
REPORT_RC=$?
set -e
if [[ -s "$TASK_ROOT/exit.json" ]]; then
    EXIT_WRITTEN=1
fi
exit "$REPORT_RC"
