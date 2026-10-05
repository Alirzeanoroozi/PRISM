#!/usr/bin/env bash
#SBATCH --signal=B:USR1@1800
#SBATCH --partition=ai
#SBATCH --account=ai
#SBATCH --qos=ai
#SBATCH --gres=gpu:ampere_a40:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=01:00:00

# Generic one-worker VALAR template. The adapter substitutes @@ markers.
# It starts qwen-llm directly inside this allocation; it never calls sbatch.
# shellcheck disable=SC1091
# Source the cluster module setup before nounset: the current /etc/bashrc
# references BASHRCSOURCED while initializing non-interactive Slurm shells.
if ! type module >/dev/null 2>&1 && [[ -f /etc/bashrc ]]; then
    source /etc/bashrc
fi
set -Eeuo pipefail

RUN_DIR=/home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6
PROJECT_PATH=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-a2ac256d43744ae084719726833b3ef6
PROMPT_PATH=/home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6/prompt-continuation-03.md
PORT=18825
PROFILE_ID=a40-q6-262k-1gpu
PROFILE_GPU_TYPE=ampere_a40
PROFILE_QUANT=q6_k
PROFILE_GPU_COUNT=1
PROFILE_CONTEXT=262144
PROFILE_PARALLEL=1
PROFILE_MODEL_PATH=/home/rshadi25/llm/models/qwen3.8-27b-q6k/Qwen3.8-27B-Q6_K.gguf
PROFILE_SERVER_BINARY=/home/rshadi25/llm/source/llama.cpp-v100/build-a40-gcc12/bin/llama-server
PROFILE_BUILD_COMMIT=73a43d1f6
COPILOT_SESSION_ID=78fae6af-9164-5f65-b7e5-bd9bc8792dfd
COPILOT_COMMAND=copilot
FRAMEWORK_ROOT="${VALAR_FRAMEWORK_ROOT:-/home/rshadi25/valar-agent-framework}"
QWEN_LLM_BIN="${QWEN_LLM_BIN:-qwen-llm}"
# qwen-llm run is invoked directly in this allocation; it must not submit a job.

mkdir -p "$RUN_DIR/logs" "$RUN_DIR/evidence" "$RUN_DIR/validation"
export PYTHONPATH="$FRAMEWORK_ROOT/src${PYTHONPATH:+:${PYTHONPATH}}"
cd "$PROJECT_PATH"
WORKSPACE_MODE=editing
SERVER_LOG="$RUN_DIR/logs/qwen-server.log"
PROVIDER_SMOKE_LOG="$RUN_DIR/logs/provider-smoke.json"
COPILOT_LOG="$RUN_DIR/logs/copilot.log"
COPILOT_OUTPUT_LOG="$RUN_DIR/logs/copilot-output.jsonl"
COPILOT_RUNTIME_LOG="$RUN_DIR/evidence/copilot-runtime.json"
SLURM_LOG="$RUN_DIR/logs/slurm.log"
if ! python3 -m valar_agent.evidence copy \
    --source "$FRAMEWORK_ROOT/.github/copilot-instructions.md" \
    --destination "$RUN_DIR/evidence/framework-worker-contract.md"
then
    printf 'framework worker contract could not be copied into durable run evidence\n' >> "$SLURM_LOG"
    exit 1
fi
CONTRACT_SHA256_RECORD="$(sha256sum "$RUN_DIR/evidence/framework-worker-contract.md")"
python3 -m valar_agent.evidence write-text \
    --output "$RUN_DIR/evidence/framework-worker-contract.sha256" \
    --content "$CONTRACT_SHA256_RECORD"
COPILOT_HOME="$RUN_DIR/copilot-home"
export COPILOT_HOME
mkdir -p "$COPILOT_HOME"
python3 -m valar_agent.evidence write-text \
    --output "$COPILOT_HOME/settings.json" \
    --content '{"disableAllHooks": true}'
python3 -m valar_agent.evidence copy \
    --source "$COPILOT_HOME/settings.json" \
    --destination "$RUN_DIR/evidence/copilot-settings.json"
SERVER_PID=""
COPILOT_PID=""
CLOSEOUT_REQUESTED=0
SINGULARITY_SANDBOX_ROOT=""

write_handoff() {
    {
        printf '\n## Worker closeout (%s)\n\n' "$(date --iso-8601=seconds)"
        printf '%s\n' "profile=$PROFILE_ID"
        printf '%s\n' "slurm_job_id=${SLURM_JOB_ID:-unknown}"
        printf '%s\n' "node=${SLURMD_NODENAME:-${SLURM_JOB_NODELIST:-unknown}}"
    printf '%s\n' "port=$PORT"
        printf '%s\n' "qwen_model_id=${QWEN_MODEL_ID:-unknown}"
        printf '%s\n' "provider_base_url=${COPILOT_PROVIDER_BASE_URL:-unknown}"
        printf '%s\n' "closeout_requested=$CLOSEOUT_REQUESTED"
    } >> "$RUN_DIR/handoff.md"
    touch "$RUN_DIR/logs/worker-closeout.marker"
}

write_checkpoint() {
    local summary="$1"
    local next_action="$2"
    python3 - "$RUN_DIR" "$summary" "$next_action" <<'PY'
import sys
from valar_agent.storage import write_checkpoint

write_checkpoint(sys.argv[1], summary=sys.argv[2], next_action=sys.argv[3])
PY
}

cleanup() {
    set +e
    if [[ -n "$COPILOT_PID" ]] && kill -0 "$COPILOT_PID" 2>/dev/null; then
        kill -TERM "$COPILOT_PID" 2>/dev/null || true
        wait "$COPILOT_PID" 2>/dev/null || true
    fi
    if [[ -n "$SERVER_PID" ]] && kill -0 "$SERVER_PID" 2>/dev/null; then
        kill -TERM "$SERVER_PID" 2>/dev/null || true
        wait "$SERVER_PID" 2>/dev/null || true
    fi
    if [[ -n "$SINGULARITY_SANDBOX_ROOT" && -d "$SINGULARITY_SANDBOX_ROOT" ]]; then
        rm -rf "$SINGULARITY_SANDBOX_ROOT"
    fi
}

on_usr1() {
    CLOSEOUT_REQUESTED=1
    printf '%s walltime warning received\n' "$(date --iso-8601=seconds)" >> "$SLURM_LOG"
    write_handoff
    if [[ -n "$COPILOT_PID" ]] && kill -0 "$COPILOT_PID" 2>/dev/null; then
        kill -TERM "$COPILOT_PID" 2>/dev/null || true
    fi
}

on_term() {
    CLOSEOUT_REQUESTED=1
    write_handoff
    cleanup
    exit 143
}

trap on_usr1 USR1
trap on_term TERM INT
trap 'write_handoff; cleanup' EXIT

if ! type module >/dev/null 2>&1; then
    printf 'CUDA module command is unavailable\n' >> "$SLURM_LOG"
    exit 2
fi
module load cuda/12.6.3 >> "$SLURM_LOG" 2>&1

printf 'profile=%s gpu=%s:%s quant=%s context=%s port=%s\n' \
    "$PROFILE_ID" "$PROFILE_GPU_TYPE" "$PROFILE_GPU_COUNT" "$PROFILE_QUANT" "$PROFILE_CONTEXT" "$PORT" \
    >> "$SLURM_LOG"

if ! python3 - "$PORT" <<'PY' >> "$SLURM_LOG" 2>&1
import socket
import sys

port = int(sys.argv[1])
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
    probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    probe.bind(("127.0.0.1", port))
print(f"port preflight passed: 127.0.0.1:{port}")
PY
then
    printf 'port preflight failed: 127.0.0.1:%s is unavailable on the allocated node\n' "$PORT" >> "$SLURM_LOG"
    write_checkpoint "The coordinated localhost port was unavailable on the allocated node before Qwen startup." "Reconcile the reservation, choose a new coordinated port, and resume the same worker run."
    exit 1
fi

"$QWEN_LLM_BIN" run \
    --gpu-type "$PROFILE_GPU_TYPE" \
    --quant "$PROFILE_QUANT" \
    --gpus "$PROFILE_GPU_COUNT" \
    --context "$PROFILE_CONTEXT" \
    --parallel "$PROFILE_PARALLEL" \
    --port "$PORT" > "$SERVER_LOG" 2>&1 &
SERVER_PID=$!

for attempt in $(seq 1 120); do
    if "$QWEN_LLM_BIN" health --port "$PORT" >> "$SERVER_LOG" 2>&1; then
        break
    fi
    if ! kill -0 "$SERVER_PID" 2>/dev/null; then
        printf 'qwen server exited before health check\n' >> "$SERVER_LOG"
        exit 1
    fi
    sleep 2
    if [[ "$attempt" -eq 120 ]]; then
        printf 'qwen server health timeout\n' >> "$SERVER_LOG"
        exit 1
    fi
done

MODEL_DISCOVERY_LOG="$RUN_DIR/logs/qwen-models.json"
if ! "$QWEN_LLM_BIN" models --port "$PORT" --json > "$MODEL_DISCOVERY_LOG" 2>> "$SERVER_LOG"; then
    printf 'local provider verification failed: /v1/models unavailable\n' >> "$SLURM_LOG"
    touch "$RUN_DIR/logs/provider-verification-failed.marker"
    write_checkpoint "Local Qwen health succeeded but /v1/models provider verification failed." "Reconnect to the local Qwen server, recover the same Copilot session, and rerun provider verification before project work."
    exit 1
fi
QWEN_MODEL_ID="$(python3 - "$MODEL_DISCOVERY_LOG" <<'PY'
import json
import sys

payload = json.loads(open(sys.argv[1], encoding="utf-8").read())
records = payload.get("data") if isinstance(payload, dict) else None
if not isinstance(records, list) or not records:
    raise SystemExit("local Qwen model discovery returned no models")
model_id = records[0].get("id") if isinstance(records[0], dict) else None
if not isinstance(model_id, str) or not model_id.strip():
    raise SystemExit("local Qwen model discovery returned an invalid model ID")
print(model_id)
PY
)"
export QWEN_MODEL_ID
export COPILOT_OFFLINE="true"
export COPILOT_PROVIDER_TYPE="openai"
export COPILOT_PROVIDER_BASE_URL="http://127.0.0.1:${PORT}/v1"
export COPILOT_MODEL="$QWEN_MODEL_ID"
export COPILOT_PROVIDER_MODEL_ID="$QWEN_MODEL_ID"
export COPILOT_PROVIDER_WIRE_MODEL="$QWEN_MODEL_ID"
export COPILOT_SESSION_ID
export COPILOT_PROVIDER_MAX_PROMPT_TOKENS="${COPILOT_PROVIDER_MAX_PROMPT_TOKENS:-245760}"
export COPILOT_PROVIDER_MAX_OUTPUT_TOKENS="${COPILOT_PROVIDER_MAX_OUTPUT_TOKENS:-8192}"
export COPILOT_PROVIDER_API_KEY=""

if ! python3 -m valar_agent.provider_smoke \
    --port "$PORT" \
    --model "$QWEN_MODEL_ID" \
    --output "$PROVIDER_SMOKE_LOG"
then
    printf 'local provider smoke request failed; Copilot was not started\n' >> "$SLURM_LOG"
    touch "$RUN_DIR/logs/provider-smoke-failed.marker"
    write_checkpoint "Qwen health and model discovery passed, but the direct local /v1/chat/completions smoke request failed." "Inspect logs/provider-smoke.json, repair the Qwen endpoint or model binding, and resume this worker before launching Copilot."
    exit 1
fi

NODE_VERSION="$(node --version 2>&1 || true)"
COPILOT_PATH="$(command -v "$COPILOT_COMMAND" || true)"
if [[ -z "$COPILOT_PATH" ]]; then
    python3 -m valar_agent.copilot_runtime \
        --output "$COPILOT_RUNTIME_LOG" \
        --status MISSING \
        --configured-command "$COPILOT_COMMAND" \
        --node-version "$NODE_VERSION"
    printf 'Copilot command unavailable: %s\n' "$COPILOT_COMMAND" > "$COPILOT_LOG"
    write_checkpoint "The configured Copilot command was not available inside the Slurm allocation." "Inspect PATH and the selected profile's copilot_command, then resume the same worker after resolving the runtime installation."
    exit 127
fi
set +e
COPILOT_VERSION_OUTPUT="$("$COPILOT_COMMAND" --version 2>&1)"
COPILOT_VERSION_STATUS=$?
COPILOT_HELP_LOG="$RUN_DIR/logs/copilot-help-providers.txt"
"$COPILOT_COMMAND" help providers > "$COPILOT_HELP_LOG" 2>&1
COPILOT_HELP_STATUS=$?
COPILOT_COMMANDS_LOG="$RUN_DIR/logs/copilot-help-commands.txt"
"$COPILOT_COMMAND" help commands > "$COPILOT_COMMANDS_LOG" 2>&1
COPILOT_COMMANDS_STATUS=$?
set -e
python3 -m valar_agent.copilot_runtime \
    --output "$COPILOT_RUNTIME_LOG" \
    --status FOUND \
    --configured-command "$COPILOT_COMMAND" \
    --resolved-path "$COPILOT_PATH" \
    --version-exit-code "$COPILOT_VERSION_STATUS" \
    --version-output "$COPILOT_VERSION_OUTPUT" \
    --provider-help-exit-code "$COPILOT_HELP_STATUS" \
    --node-version "$NODE_VERSION"
if [[ "$COPILOT_VERSION_STATUS" -ne 0 || "$COPILOT_HELP_STATUS" -ne 0 || "$COPILOT_COMMANDS_STATUS" -ne 0 ]]; then
    printf 'Copilot runtime preflight failed: version_status=%s provider_help_status=%s command_help_status=%s\n' \
        "$COPILOT_VERSION_STATUS" "$COPILOT_HELP_STATUS" "$COPILOT_COMMANDS_STATUS" >> "$SLURM_LOG"
    write_checkpoint "The Copilot runtime/version preflight failed before prompt dispatch." "Inspect evidence/copilot-runtime.json and logs/copilot-help-providers.txt, then resume after correcting the compute-node Copilot installation."
    exit 1
fi

COPILOT_OBJECTIVE_COMMAND=""
COPILOT_OBJECTIVE_FLAGS=()
if grep -Eq '^[[:space:]]*/autopilot[[:space:]]' "$COPILOT_COMMANDS_LOG"; then
    COPILOT_OBJECTIVE_COMMAND="/autopilot"
    COPILOT_OBJECTIVE_FLAGS=(
        --mode autopilot
        --max-autopilot-continues "${VALAR_MAX_AUTOPILOT_CONTINUES:-5}"
        --experimental
    )
elif grep -Eq '^[[:space:]]*/goal[[:space:]]' "$COPILOT_COMMANDS_LOG"; then
    COPILOT_OBJECTIVE_COMMAND="/goal"
else
    printf 'Copilot exposes neither /autopilot nor legacy /goal objective command\n' >> "$SLURM_LOG"
    write_checkpoint "Copilot started but did not advertise a supported persistent objective command." "Inspect logs/copilot-help-commands.txt, then resume the same worker after selecting a compatible Copilot CLI runtime."
    exit 1
fi
printf 'copilot_objective_command=%s\n' "$COPILOT_OBJECTIVE_COMMAND" >> "$SLURM_LOG"
COPILOT_SANDBOX_PREFIX=()
if [[ "$WORKSPACE_MODE" == "read_only" ]]; then
    BWRAP_COMMAND="$(command -v bwrap || true)"
    if [[ -n "$BWRAP_COMMAND" ]]; then
        COPILOT_SANDBOX_PREFIX=(
            "$BWRAP_COMMAND"
            --die-with-parent
            --new-session
            --ro-bind / /
            --dev /dev
            --proc /proc
            --tmpfs /tmp
            --ro-bind "$PROJECT_PATH" "$PROJECT_PATH"
            --bind "$RUN_DIR/logs" "$RUN_DIR/logs"
            --bind "$RUN_DIR/evidence" "$RUN_DIR/evidence"
            --bind "$RUN_DIR/validation" "$RUN_DIR/validation"
            --bind "$RUN_DIR/checkpoint.md" "$RUN_DIR/checkpoint.md"
            --bind "$RUN_DIR/handoff.md" "$RUN_DIR/handoff.md"
            --bind "$RUN_DIR/report.md" "$RUN_DIR/report.md"
            --bind "$RUN_DIR/decisions.jsonl" "$RUN_DIR/decisions.jsonl"
            --bind "$COPILOT_HOME" "$COPILOT_HOME"
            --chdir "$PROJECT_PATH"
            --
        )
        printf 'read_only_boundary=bwrap project=ro run_writable_subpaths=logs,evidence,validation,checkpoint,handoff,report,decisions,copilot-home\n' >> "$SLURM_LOG"
        python3 - "$RUN_DIR/evidence/workspace-boundary.json" "$PROJECT_PATH" "$RUN_DIR" <<'PY'
import json
import sys
from pathlib import Path
from valar_agent.evidence import atomic_write_json

output, project_path, run_path = sys.argv[1:]
atomic_write_json(
    Path(output),
    {
        "mode": "read_only",
        "enforcement": "bubblewrap",
        "project_mount": "read_only",
        "run_mount": "controlled_subpaths",
        "project_path": project_path,
        "run_path": run_path,
        "writable_run_subpaths": [
            "logs",
            "evidence",
            "validation",
            "checkpoint.md",
            "handoff.md",
            "report.md",
            "decisions.jsonl",
            "copilot-home",
        ],
    },
)
PY
    else
        if ! type singularity >/dev/null 2>&1; then
            module load singularity/4.3.2 >> "$SLURM_LOG" 2>&1 || true
        fi
        SINGULARITY_COMMAND="$(command -v singularity || true)"
        if [[ -z "$SINGULARITY_COMMAND" ]]; then
            printf 'read-only worker requires bwrap or singularity, but neither is available\n' >> "$SLURM_LOG"
            write_checkpoint "The read-only worker could not establish an OS-enforced project boundary because neither Bubblewrap nor Singularity is available." "Load singularity/4.3.2 or expose bwrap on the compute node, then retry this same run without modifying the canonical project."
            exit 1
        fi
        SINGULARITY_SANDBOX_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/valar-singularity-root.XXXXXX")"
        mkdir -p "$SINGULARITY_SANDBOX_ROOT"/{bin,dev,etc,home,lib,lib64,opt,proc,scratch,tmp,usr,var/tmp}
        mkdir -p "$SINGULARITY_SANDBOX_ROOT$PROJECT_PATH" "$SINGULARITY_SANDBOX_ROOT$RUN_DIR"
        COPILOT_BIN_DIR="$(dirname "$COPILOT_PATH")"
        NODE_BIN_DIR="$(dirname "$(command -v node)")"
        for writable_dir in logs evidence validation copilot-home; do
            mkdir -p "$SINGULARITY_SANDBOX_ROOT$RUN_DIR/$writable_dir"
        done
        for writable_file in checkpoint.md handoff.md report.md decisions.jsonl; do
            touch "$SINGULARITY_SANDBOX_ROOT$RUN_DIR/$writable_file"
        done
        COPILOT_SANDBOX_PREFIX=(
            "$SINGULARITY_COMMAND" exec
            --containall
            --no-home
            --bind /usr:/usr:ro
            --bind /bin:/bin:ro
            --bind /lib:/lib:ro
            --bind /lib64:/lib64:ro
            --bind /etc:/etc:ro
            --bind /opt:/opt:ro
            --bind /home:/home:ro
            --bind /scratch:/scratch:ro
            --env "PATH=$COPILOT_BIN_DIR:$NODE_BIN_DIR:/opt/ohpc/pub/compiler/conda3/latest/bin:/usr/local/bin:/usr/bin:/bin"
            --env "COPILOT_HOME=$COPILOT_HOME"
            --env "COPILOT_OFFLINE=$COPILOT_OFFLINE"
            --env "COPILOT_PROVIDER_TYPE=$COPILOT_PROVIDER_TYPE"
            --env "COPILOT_PROVIDER_BASE_URL=$COPILOT_PROVIDER_BASE_URL"
            --env "COPILOT_MODEL=$QWEN_MODEL_ID"
            --env "COPILOT_PROVIDER_MODEL_ID=$QWEN_MODEL_ID"
            --env "COPILOT_PROVIDER_WIRE_MODEL=$QWEN_MODEL_ID"
            --env "COPILOT_SESSION_ID=$COPILOT_SESSION_ID"
            --env "COPILOT_PROVIDER_MAX_PROMPT_TOKENS=$COPILOT_PROVIDER_MAX_PROMPT_TOKENS"
            --env "COPILOT_PROVIDER_MAX_OUTPUT_TOKENS=$COPILOT_PROVIDER_MAX_OUTPUT_TOKENS"
            --bind "$PROJECT_PATH:$PROJECT_PATH:ro"
            --bind "$RUN_DIR:$RUN_DIR:ro"
            --bind "$RUN_DIR/logs:$RUN_DIR/logs:rw"
            --bind "$RUN_DIR/evidence:$RUN_DIR/evidence:rw"
            --bind "$RUN_DIR/validation:$RUN_DIR/validation:rw"
            --bind "$RUN_DIR/checkpoint.md:$RUN_DIR/checkpoint.md:rw"
            --bind "$RUN_DIR/handoff.md:$RUN_DIR/handoff.md:rw"
            --bind "$RUN_DIR/report.md:$RUN_DIR/report.md:rw"
            --bind "$RUN_DIR/decisions.jsonl:$RUN_DIR/decisions.jsonl:rw"
            --bind "$COPILOT_HOME:$COPILOT_HOME:rw"
            --pwd "$PROJECT_PATH"
            "$SINGULARITY_SANDBOX_ROOT"
        )
        printf 'read_only_boundary=singularity project=ro run_writable_subpaths=logs,evidence,validation,checkpoint,handoff,report,decisions,copilot-home\n' >> "$SLURM_LOG"
        python3 - "$RUN_DIR/evidence/workspace-boundary.json" "$PROJECT_PATH" "$RUN_DIR" "$SINGULARITY_SANDBOX_ROOT" <<'PY'
import json
import sys
from pathlib import Path
from valar_agent.evidence import atomic_write_json

output, project_path, run_path, sandbox_root = sys.argv[1:]
atomic_write_json(
    Path(output),
    {
        "mode": "read_only",
        "enforcement": "singularity",
        "project_mount": "read_only",
        "run_mount": "controlled_subpaths",
        "project_path": project_path,
        "run_path": run_path,
        "sandbox_root": sandbox_root,
        "writable_run_subpaths": [
            "logs",
            "evidence",
            "validation",
            "checkpoint.md",
            "handoff.md",
            "report.md",
            "decisions.jsonl",
            "copilot-home",
        ],
    },
)
PY
    fi
elif [[ "$WORKSPACE_MODE" == "editing" ]]; then
    printf 'read_only_boundary=not_applicable workspace=editing\n' >> "$SLURM_LOG"
else
    printf 'unknown workspace mode: %s\n' "$WORKSPACE_MODE" >> "$SLURM_LOG"
    write_checkpoint "The worker manifest contains an unknown workspace mode, so project execution was refused." "Repair the durable workspace metadata and resume the same bounded run."
    exit 1
fi
python3 -m valar_agent.evidence write-text \
    --output "$RUN_DIR/evidence/copilot-objective-dispatch.txt" \
    --content "$(printf 'command=%s\nexperimental=%s\nhelp_log=logs/copilot-help-commands.txt\n' \
        "$COPILOT_OBJECTIVE_COMMAND" \
        "$([[ ${#COPILOT_OBJECTIVE_FLAGS[@]} -gt 0 ]] && printf true || printf false)")"

python3 - "$RUN_DIR/evidence/copilot-provider-binding.json" "$QWEN_MODEL_ID" "$PORT" "$PROFILE_ID" "$PROFILE_QUANT" "$PROFILE_GPU_TYPE" "$PROFILE_GPU_COUNT" "$PROFILE_CONTEXT" "$COPILOT_HOME" <<'PY'
import sys
from pathlib import Path
from valar_agent.evidence import atomic_write_json

output, model_id, port, profile_id, quant, gpu_type, gpu_count, context, copilot_home = sys.argv[1:]
payload = {
    "copilot_home": copilot_home,
    "copilot_offline": True,
    "copilot_provider_type": "openai",
    "copilot_provider_base_url": f"http://127.0.0.1:{port}/v1",
    "copilot_provider_api_key": "",
    "copilot_provider_model_id": model_id,
    "copilot_provider_wire_model": model_id,
    "copilot_model": model_id,
    "copilot_model_argument": model_id,
    "copilot_settings_path": str(Path(copilot_home) / "settings.json"),
    "disable_all_hooks": True,
    "qwen_profile_id": profile_id,
    "qwen_quantization": quant,
    "gpu_type": gpu_type,
    "gpu_count": int(gpu_count),
    "context": int(context),
    "server_parallel": int("1"),
}
atomic_write_json(Path(output), payload)
PY

RESOURCE_INVENTORY_LOG="$RUN_DIR/logs/copilot-resources.json"
if ! "$COPILOT_COMMAND" plugins list --kind mcp --kind skill --json > "$RESOURCE_INVENTORY_LOG" 2>&1; then
    printf 'Copilot resource inventory failed; refusing to start project work\n' >> "$SLURM_LOG"
    touch "$RUN_DIR/logs/provider-verification-failed.marker"
    write_checkpoint "Local Qwen provider was verified, but Copilot skill/MCP inventory failed." "Inspect the Copilot installation and resource inventory, then resume the same bounded objective session."
    exit 1
fi

MCP_ARGS=()
if [[ "${VALAR_ENABLE_NLM_MCP:-0}" == "1" ]]; then
    while IFS= read -r discovered_mcp; do
        [[ -n "$discovered_mcp" ]] || continue
        MCP_ARGS+=(--enable-mcp-server "$discovered_mcp")
    done < <(python3 - "$RESOURCE_INVENTORY_LOG" <<'PY'
import json
import sys

payload = json.loads(open(sys.argv[1], encoding="utf-8").read())
names = set()
for plugin in payload.get("plugins", []) if isinstance(payload, dict) else []:
    if not isinstance(plugin, dict):
        continue
    kind = str(plugin.get("kind", "")).lower()
    name = plugin.get("name") or plugin.get("id") or plugin.get("server_name")
    if kind == "mcp" and isinstance(name, str) and ("nlm" in name.lower() or "notebooklm" in name.lower()):
        names.add(name)
for name in sorted(names):
    print(name)
PY
    )
fi

export VALAR_RUN_DIR="$RUN_DIR"
export VALAR_PROVIDER_TYPE="$COPILOT_PROVIDER_TYPE"
export VALAR_PROVIDER_BASE_URL="$COPILOT_PROVIDER_BASE_URL"
export VALAR_QWEN_MODEL_ID="$QWEN_MODEL_ID"
export VALAR_QWEN_SERVER_BINARY="$PROFILE_SERVER_BINARY"
export VALAR_QWEN_BUILD_COMMIT="$PROFILE_BUILD_COMMIT"
export VALAR_QUANTIZATION="$PROFILE_QUANT"
export VALAR_GPU_TYPE="$PROFILE_GPU_TYPE"
export VALAR_GPU_COUNT="$PROFILE_GPU_COUNT"
export VALAR_CONTEXT="$PROFILE_CONTEXT"
export VALAR_RESOURCE_INVENTORY_PATH="logs/copilot-resources.json"
python3 - <<'PY'
import os
from valar_agent.worker import record_provider_binding

record_provider_binding(
    os.environ["VALAR_RUN_DIR"],
    provider_type=os.environ["VALAR_PROVIDER_TYPE"],
    provider_base_url=os.environ["VALAR_PROVIDER_BASE_URL"],
    model_id=os.environ["VALAR_QWEN_MODEL_ID"],
    server_binary=os.environ["VALAR_QWEN_SERVER_BINARY"],
    build_commit=os.environ["VALAR_QWEN_BUILD_COMMIT"],
    quantization=os.environ["VALAR_QUANTIZATION"],
    gpu_type=os.environ["VALAR_GPU_TYPE"],
    gpu_count=int(os.environ["VALAR_GPU_COUNT"]),
    context=int(os.environ["VALAR_CONTEXT"]),
    resource_inventory_path=os.environ["VALAR_RESOURCE_INVENTORY_PATH"],
)
PY
touch "$RUN_DIR/logs/provider-verified.marker"

if ! command -v "$COPILOT_COMMAND" >/dev/null 2>&1; then
    printf 'Copilot command unavailable: %s\n' "$COPILOT_COMMAND" > "$COPILOT_LOG"
    exit 127
fi
COPILOT_ARGS=(
    --session-id "$COPILOT_SESSION_ID" \
    --add-dir "$PROJECT_PATH" \
    --add-dir "$RUN_DIR" \
    --log-dir "$RUN_DIR/logs/copilot" \
    --allow-all-tools \
    --no-ask-user \
    --no-remote \
    --no-remote-export \
    --no-auto-update \
    --log-level debug \
    --output-format json \
    "${COPILOT_OBJECTIVE_FLAGS[@]}" \
    --model "$QWEN_MODEL_ID" \
    "${MCP_ARGS[@]}" \
    --prompt "$COPILOT_OBJECTIVE_COMMAND $(cat "$PROMPT_PATH")"
)
"${COPILOT_SANDBOX_PREFIX[@]}" "$COPILOT_COMMAND" "${COPILOT_ARGS[@]}" > "$COPILOT_OUTPUT_LOG" 2> "$COPILOT_LOG" &
COPILOT_PID=$!
set +e
wait "$COPILOT_PID"
COPILOT_STATUS=$?
set -e
if [[ "$COPILOT_STATUS" -ne 0 ]]; then
    if [[ "$CLOSEOUT_REQUESTED" -eq 1 ]]; then
        write_checkpoint "Walltime warning interrupted Copilot with status $COPILOT_STATUS; durable worker evidence was preserved for resume." "Reconcile this job, inspect report.md, handoff.md, logs/copilot-output.jsonl, and evidence/, then resume the same bounded goal only if the success criteria are not already recorded."
        exit 0
    fi
    write_checkpoint "Copilot non-interactive objective dispatch exited with status $COPILOT_STATUS before worker completion evidence was recorded." "Inspect logs/copilot.log and logs/copilot-output.jsonl, then resume the same bounded goal only after the failure cause is understood."
    exit "$COPILOT_STATUS"
fi
GOAL_OBJECTIVE_PATH="$(find "$COPILOT_HOME" -path "*/$COPILOT_SESSION_ID/*objective*.json" -type f -print -quit 2>/dev/null || true)"
OBJECTIVE_ACTIVATION_EVIDENCE="$RUN_DIR/evidence/copilot-objective-activation.json"
if ! python3 - "$GOAL_OBJECTIVE_PATH" "$COPILOT_OUTPUT_LOG" "$PROMPT_PATH" "$COPILOT_SESSION_ID" "$COPILOT_OBJECTIVE_COMMAND" "$RUN_DIR" "$OBJECTIVE_ACTIVATION_EVIDENCE" <<'PY'
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from valar_agent.evidence import atomic_write_json

objective_path, output_path, prompt_path, session_id, command, run_path, evidence_path = sys.argv[1:]
prompt_file = Path(prompt_path)
output_file = Path(output_path)
prompt = prompt_file.read_text(encoding="utf-8")
prompt_sha256 = hashlib.sha256(prompt.encode("utf-8")).hexdigest()
run_id = Path(run_path).name

def artifact_matches() -> bool:
    if not objective_path:
        return False
    try:
        payload = json.loads(Path(objective_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    current = payload.get("current") if isinstance(payload, dict) else None
    objective = current.get("objective") if isinstance(current, dict) else None
    status = current.get("status") if isinstance(current, dict) else None
    artifact_prompt_sha256 = current.get("prompt_sha256") if isinstance(current, dict) else None
    artifact_session_id = current.get("session_id") if isinstance(current, dict) else None
    artifact_run_id = current.get("run_id") if isinstance(current, dict) else None
    return (
        status in {"active", "completed"}
        and isinstance(objective, str)
        and prompt in objective
        and run_id in objective
        and (artifact_prompt_sha256 is None or artifact_prompt_sha256 == prompt_sha256)
        and (artifact_session_id is None or artifact_session_id == session_id)
        and (artifact_run_id is None or artifact_run_id == run_id)
    )

def dispatch_event_matches() -> bool:
    if not output_file.is_file():
        return False
    for line in output_file.read_text(encoding="utf-8").splitlines():
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(value, dict) or value.get("type") != "user.message":
            continue
        data = value.get("data")
        if not isinstance(data, dict):
            continue
        content = data.get("content")
        transformed = data.get("transformedContent")
        combined = "\n".join(item for item in (content, transformed) if isinstance(item, str))
        if prompt not in combined or run_id not in combined:
            continue
        if command == "/autopilot" and "Autopilot objective:" not in combined:
            continue
        if command == "/goal" and "/goal" not in combined.lower() and "goal" not in combined.lower():
            continue
        return True
    return False

if not artifact_matches():
    raise SystemExit("no matching persisted objective artifact was found")
if not dispatch_event_matches():
    raise SystemExit("no prompt-matching objective dispatch event was found")

atomic_write_json(
    Path(evidence_path),
    {
        "mode": command.lstrip("/"),
        "objective_command": command,
        "session_id": session_id,
        "run_id": run_id,
        "prompt_path": str(prompt_file),
        "prompt_sha256": prompt_sha256,
        "objective_artifact_path": objective_path,
        "objective_artifact_verified": True,
        "dispatch_event_verified": True,
        "output_path": str(output_file),
        "recorded_at": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
    },
)
PY
then
    python3 -m valar_agent.evidence copy \
        --source "$GOAL_OBJECTIVE_PATH" \
        --destination "$RUN_DIR/evidence/copilot-goal-objective.json"
    python3 -m valar_agent.evidence copy \
        --source "$COPILOT_OUTPUT_LOG" \
        --destination "$RUN_DIR/evidence/copilot-prompt-output.jsonl"
    GOAL_ACTIVATION_PATH="evidence/copilot-objective-activation.json"
else
    printf 'Copilot did not leave objective evidence matching the persisted prompt/session\n' >> "$SLURM_LOG"
    touch "$RUN_DIR/logs/goal-activation-failed.marker"
    write_checkpoint "Copilot exited without objective evidence matching the persisted prompt and session." "Inspect logs/copilot.log and logs/copilot-output.jsonl, then recover the same Copilot session only after objective activation is understood."
    exit 1
fi
export VALAR_GOAL_ACTIVATION_PATH="$GOAL_ACTIVATION_PATH"
python3 - <<'PY'
import os
from valar_agent.worker import record_provider_binding

record_provider_binding(
    os.environ["VALAR_RUN_DIR"],
    provider_type=os.environ["VALAR_PROVIDER_TYPE"],
    provider_base_url=os.environ["VALAR_PROVIDER_BASE_URL"],
    model_id=os.environ["VALAR_QWEN_MODEL_ID"],
    server_binary=os.environ["VALAR_QWEN_SERVER_BINARY"],
    build_commit=os.environ["VALAR_QWEN_BUILD_COMMIT"],
    quantization=os.environ["VALAR_QUANTIZATION"],
    gpu_type=os.environ["VALAR_GPU_TYPE"],
    gpu_count=int(os.environ["VALAR_GPU_COUNT"]),
    context=int(os.environ["VALAR_CONTEXT"]),
    resource_inventory_path=os.environ["VALAR_RESOURCE_INVENTORY_PATH"],
    goal_activation_path=os.environ["VALAR_GOAL_ACTIVATION_PATH"],
)
PY
# A worker must leave recoverable state even when Copilot exits without
# calling a framework-specific checkpoint helper. Slurm completion is not
# validation.
if [[ ! -s "$RUN_DIR/checkpoint.md" ]]; then
    write_checkpoint \
        "Copilot objective dispatch exited after local-provider verification; inspect the recorded project evidence before treating the goal as complete." \
        "Load the manifest, checkpoint, handoff, and evidence; refresh external status and independently validate the bounded goal."
fi
write_handoff
exit "$COPILOT_STATUS"
