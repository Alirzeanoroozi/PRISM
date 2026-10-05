#!/bin/bash
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
SMOKE_ROOT="${MP_SMOKE_ROOT:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/multiprot-smoke-array-v1}"
INSTALL_ROOT="${MP_INSTALL_ROOT:-$REPO_ROOT/tmp/agent/20260713-investigation-implementation/multiprot-environment-v1}"
PYTHON3="${MP_PYTHON3:-/opt/ohpc/pub/compiler/conda3/latest/bin/python3}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
TASK_ROOT="$SMOKE_ROOT/task-${TASK_ID}"
if [[ -e "$TASK_ROOT" ]]; then
    if [[ -n "$(find "$TASK_ROOT" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        echo "refusing to reuse non-empty task directory: $TASK_ROOT" >&2
        exit 3
    fi
else
    mkdir -p "$TASK_ROOT"
fi

INPUT1=""
INPUT2=""
BIN="$INSTALL_ROOT/bin/multiprot"
BINARY="$INSTALL_ROOT/multiprot/multiprot.Linux"
PARAMS="$INSTALL_ROOT/multiprot/params.txt"
EXIT_WRITTEN=0
export TASK_ROOT TASK_ID INPUT1 INPUT2 BIN BINARY PARAMS

write_failure_exit() {
    local status="$1"
    export EXIT_STATUS="$status"
    if [[ -x "$PYTHON3" ]]; then
        "$PYTHON3" - <<'PY' || true
import json
import os
from pathlib import Path

root = Path(os.environ["TASK_ROOT"])
value = {
    "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
    "array_task_id": os.environ.get("TASK_ID", ""),
    "job_id": os.environ.get("SLURM_JOB_ID", ""),
    "return_code": int(os.environ.get("EXIT_STATUS", "1")),
    "scientific_status": "launcher_failed_before_binary",
    "inputs": {"receptor": os.environ.get("INPUT1", ""), "ligand": os.environ.get("INPUT2", "")},
    "command": [os.environ.get("BIN", ""), "receptor.pdb", "ligand.pdb"],
}
(root / "exit.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
PY
    else
        printf '{"return_code":%s,"scientific_status":"launcher_failed_before_binary"}\n' "$status" > "$TASK_ROOT/exit.json" || true
    fi
    EXIT_WRITTEN=1
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

case "$TASK_ID" in
    1)
        INPUT1="$REPO_ROOT/benchmark/data/pdbs/rigid/1fgn.pdb"
        INPUT2="$REPO_ROOT/benchmark/data/pdbs/rigid/1fgn.pdb"
        ;;
    2)
        INPUT1="$REPO_ROOT/benchmark/data/pdbs/rigid/1jpt.pdb"
        INPUT2="$REPO_ROOT/benchmark/data/pdbs/rigid/1mkf.pdb"
        ;;
    *)
        echo "unsupported array index: $TASK_ID" >&2
        exit 2
        ;;
esac

test -f "$INPUT1"
test -f "$INPUT2"
test -x "$BIN"
test -x "$BINARY"
test -f "$PARAMS"
test -f "$INSTALL_ROOT/environment_manifest.json"
test -x "$PYTHON3"
"$PYTHON3" - "$INSTALL_ROOT/environment_manifest.json" "$BINARY" <<'PY'
import hashlib
import json
import stat
import sys
from pathlib import Path

manifest = json.loads(Path(sys.argv[1]).read_text())
binary = Path(sys.argv[2])
if not manifest["status"].startswith("ready_for_standalone_multiprot"):
    raise SystemExit("environment is not ready for standalone MultiProt: " + manifest["status"])
if not manifest["dependencies"]["pymysql"]["available"] or not manifest["dependencies"]["mysqldb_compat"]["available"]:
    raise SystemExit("MySQL compatibility dependency gate failed")
if not (stat.S_IMODE(binary.stat().st_mode) & 0o111):
    raise SystemExit("staged MultiProt binary is not executable")
digest = hashlib.sha256(binary.read_bytes()).hexdigest()
if digest != manifest["binary_sha256"]:
    raise SystemExit("staged MultiProt binary hash does not match environment manifest")
PY

STAGED1="$TASK_ROOT/receptor.pdb"
STAGED2="$TASK_ROOT/ligand.pdb"
cp "$INPUT1" "$STAGED1"
cp "$INPUT2" "$STAGED2"
cp "$PARAMS" "$TASK_ROOT/params.txt"
{
    printf 'source_path\tstaged_path\tsha256\n'
    printf '%s\t%s\t%s\n' "$INPUT1" "$STAGED1" "$(sha256sum "$STAGED1" | cut -d ' ' -f1)"
    printf '%s\t%s\t%s\n' "$INPUT2" "$STAGED2" "$(sha256sum "$STAGED2" | cut -d ' ' -f1)"
    printf '%s\t%s\t%s\n' "$PARAMS" "$TASK_ROOT/params.txt" "$(sha256sum "$TASK_ROOT/params.txt" | cut -d ' ' -f1)"
} > "$TASK_ROOT/input_manifest.tsv"
printf '%s\n' "$BIN receptor.pdb ligand.pdb" > "$TASK_ROOT/command.txt"
printf 'job_id\tarray_job_id\tarray_task_id\tinput1\tinput2\tbin\tparams\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
"${SLURM_JOB_ID:-}" "${SLURM_ARRAY_JOB_ID:-}" "$TASK_ID" "$INPUT1" "$INPUT2" "$BIN" "$PARAMS" > "$TASK_ROOT/parameters.tsv"

START_EPOCH=$(date +%s)
cd "$TASK_ROOT"
set +e
"$BIN" receptor.pdb ligand.pdb > stdout.log 2> stderr.log
RC=$?
set -e
END_EPOCH=$(date +%s)

export TASK_ROOT TASK_ID START_EPOCH END_EPOCH RC INPUT1 INPUT2 BIN BINARY PARAMS
"$PYTHON3" - <<'PY'
import json
import hashlib
import os
from pathlib import Path

root = Path(os.environ["TASK_ROOT"])

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

outputs = {}
for name in ("stdout.log", "stderr.log", "log_multiprot.txt", "2_sol.res", "2_sets.res", "3_sol.res", "4_sol.res"):
    path = root / name
    outputs[name] = {
        "path": str(path),
        "exists": path.exists(),
        "size": path.stat().st_size if path.exists() else 0,
        "sha256": sha256(path) if path.exists() else "",
    }
value = {
    "job_id": os.environ.get("SLURM_JOB_ID", ""),
    "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
    "array_task_id": os.environ["TASK_ID"],
    "start_epoch": int(os.environ["START_EPOCH"]),
    "end_epoch": int(os.environ["END_EPOCH"]),
    "return_code": int(os.environ["RC"]),
    "command": [os.environ["BIN"], "receptor.pdb", "ligand.pdb"],
    "inputs": {"receptor": os.environ["INPUT1"], "ligand": os.environ["INPUT2"]},
    "parameters": {"source": os.environ["PARAMS"], "staged": str(root / "params.txt"), "sha256": sha256(root / "params.txt")},
    "binary_sha256": sha256(Path(os.environ["BINARY"])),
    "requested_resources": {
        "nodes": os.environ.get("SLURM_JOB_NUM_NODES", "1"),
        "ntasks": os.environ.get("SLURM_NTASKS", "1"),
        "cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK", "2"),
        "mem": "2G",
        "time": "00:05:00",
    },
    "outputs": outputs,
    "scientific_status": "smoke_pass" if int(os.environ["RC"]) == 0 and outputs["log_multiprot.txt"]["exists"] else "smoke_failed",
}
(root / "exit.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
PY
EXIT_WRITTEN=1
exit "$RC"
