#!/bin/bash
# Submit one array task per ten-pair batch:
#   sbatch --export=ALL,PIPELINE=current,BATCH_ROOT=...,STAGED_PDB_DIR=... \
#     benchmark/scripts/submit_comparison_batches.sbatch
# Repeat with PIPELINE=legacy. Each task writes an explicit unavailable-input
# file and runs only pairs whose four-character PDB inputs were staged.
#SBATCH --job-name=prism-compare
#SBATCH --array=1-26
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=72:00:00
#SBATCH --output=tmp/agent/20260712-multichain-full-comparison/slurm-%A_%a.out
#SBATCH --error=tmp/agent/20260712-multichain-full-comparison/slurm-%A_%a.err

set -euo pipefail

: "${PIPELINE:?set PIPELINE=current or legacy}"
: "${BATCH_ROOT:?set BATCH_ROOT to the prepared batch directory}"
: "${STAGED_PDB_DIR:?set STAGED_PDB_DIR to login-side staged PDBs}"

template_workflow="${TEMPLATE_WORKFLOW:-$([[ "$PIPELINE" == current ]] && echo json || echo legacy)}"
tm_threshold="${PRISM_TM_SCORE_THRESHOLD:-0.4}"
case "$tm_threshold" in
  ''|*[!0-9.]* ) echo "invalid PRISM_TM_SCORE_THRESHOLD=$tm_threshold" >&2; exit 2 ;;
esac
if ! LC_NUMERIC=C awk "BEGIN{exit(!($tm_threshold+0==$tm_threshold && $tm_threshold>0))}" 2>/dev/null; then
  echo "invalid (malformed) PRISM_TM_SCORE_THRESHOLD=$tm_threshold" >&2; exit 2
fi
if [[ "$template_workflow" != json && "$template_workflow" != legacy ]]; then
  echo "TEMPLATE_WORKFLOW must be json or legacy" >&2
  exit 2
fi

repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
absolute_from_repo() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s/%s\n' "$repo_root" "$1" ;;
  esac
}
# Slurm launches this script from a spool directory, not from the repository.
# Normalize caller-supplied relative roots before any array-task I/O.
BATCH_ROOT="$(absolute_from_repo "$BATCH_ROOT")"
STAGED_PDB_DIR="$(absolute_from_repo "$STAGED_PDB_DIR")"
if [[ -n "${RUN_ROOT:-}" ]]; then
  RUN_ROOT="$(absolute_from_repo "$RUN_ROOT")"
fi
if [[ -n "${TEMPLATE_LIST:-}" ]]; then
  TEMPLATE_LIST="$(absolute_from_repo "$TEMPLATE_LIST")"
  [[ -f "$TEMPLATE_LIST" ]] || { echo "missing TEMPLATE_LIST: $TEMPLATE_LIST" >&2; exit 127; }
fi
legacy_root="${LEGACY_ROOT:-/scratch/rshadi25/GitHub/prism-oldversion}"
legacy_template_root="${LEGACY_TEMPLATE_ROOT:-$legacy_root/template_default}"
batch_id=$(printf 'batch_%04d' "$SLURM_ARRAY_TASK_ID")
batch_dir="$BATCH_ROOT/$batch_id"
run_root="${RUN_ROOT:-$BATCH_ROOT/runs}/$PIPELINE/$batch_id"
mkdir -p "$run_root/logs" "$run_root/status"
started_epoch=$(date +%s)
completion_python="${PIPELINE_COMPLETION_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
termination_signal=""
on_term() { termination_signal="SIGTERM"; exit 143; }
on_int() { termination_signal="SIGINT"; exit 130; }
finish_status() {
  rc=$?
  elapsed=$(( $(date +%s) - started_epoch ))
  if [[ "$PIPELINE" == current ]]; then
    set +e
    "$completion_python" "$repo_root/benchmark/scripts/pipeline_completion_contract.py" \
      --run-root "$run_root" \
      --return-code "$rc" \
      --termination-signal "$termination_signal" \
      --output "$run_root/status/exit.json" \
      --pipeline "$PIPELINE" \
      --aligner "${ALIGNER:-tmalign}" \
      --elapsed-seconds "$elapsed" \
      --slurm-job-id "${SLURM_JOB_ID:-}" \
      --array-task-id "${SLURM_ARRAY_TASK_ID:-}"
    completion_rc=$?
    set -e
    if [[ "$completion_rc" -ne 0 ]]; then
      printf '{"return_code":%d,"scientific_status":"incomplete","reason":"completion_classifier_failed","elapsed_seconds":%d,"pipeline":"%s","aligner":"%s","slurm_job_id":"%s","array_task_id":"%s"}\n' \
        "$rc" "$elapsed" "$PIPELINE" "${ALIGNER:-tmalign}" "${SLURM_JOB_ID:-}" "${SLURM_ARRAY_TASK_ID:-}" > "$run_root/status/exit.json"
    fi
  else
    scientific_status=completed
    [[ "$rc" -eq 0 ]] || scientific_status=failed
    if [[ "$rc" -eq 0 ]]; then
      passed_file="$run_root/jobs/$batch_id/transformation/passedFiles"
      [[ -s "$passed_file" ]] || scientific_status=completed_no_predictions
    fi
    printf '{"return_code":%d,"scientific_status":"%s","elapsed_seconds":%d,"pipeline":"%s","aligner":"%s","slurm_job_id":"%s","array_task_id":"%s"}\n' \
      "$rc" "$scientific_status" "$elapsed" "$PIPELINE" "${ALIGNER:-tmalign}" \
      "${SLURM_JOB_ID:-}" "${SLURM_ARRAY_TASK_ID:-}" > "$run_root/status/exit.json"
  fi
  exit "$rc"
}
trap on_term TERM
trap on_int INT
trap finish_status EXIT

python3 - "$batch_dir/inputs.csv" "$STAGED_PDB_DIR" "$run_root/status/available_inputs.csv" "$run_root/status/unavailable_inputs.csv" <<'PY'
import csv, shutil, sys
from pathlib import Path
inputs, staged, available, unavailable = map(Path, sys.argv[1:])
staged.mkdir(parents=True, exist_ok=True)
with inputs.open(newline='') as handle:
    rows = list(csv.DictReader(handle))
fields = list(rows[0]) if rows else ['pair_id', 'Receptor', 'Ligand']
ok, bad = [], []
for row in rows:
    missing = []
    for role in ('Receptor', 'Ligand'):
        pdb4 = row[role][:4].lower()
        source = staged / f'{pdb4}.pdb'
        if not source.exists():
            missing.append(pdb4)
    (ok if not missing else bad).append(row | {'missing_pdb_ids': ','.join(missing)})
with available.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader(); writer.writerows({k: row[k] for k in fields} for row in ok)
bad_fields = fields + ['missing_pdb_ids']
with unavailable.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=bad_fields)
    writer.writeheader(); writer.writerows(bad)
print(f'available={len(ok)} unavailable={len(bad)}')
PY

if [[ "$PIPELINE" == "current" ]]; then
  cd "$run_root"
  mkdir -p templates processed
  ln -sfn "$repo_root/src" src
  ln -sfn "$repo_root/external_tools" external_tools
  ln -sfn "$repo_root/prism.py" prism.py
  # The current alignment stage needs template coordinate files in addition
  # to interface metadata. Keep the retained template assets read-only by
  # linking the repository directories into the per-run workspace.
  ln -sfn "$repo_root/new_template/template/pdbs" templates/pdbs
  ln -sfn "$repo_root/new_template/template/rsas" templates/rsas
  ln -sfn "$repo_root/new_template/template/hotspots" templates/hotspots
  ln -sfn "$repo_root/new_template/template/contacts" templates/contacts
  ln -sfn "$repo_root/new_template/template/interfaces" templates/interfaces
  ln -sfn "$repo_root/new_template/template/interfaces_lists" templates/interfaces_lists
  if [[ -n "${TEMPLATE_LIST:-}" ]]; then
    cp "$TEMPLATE_LIST" templates/calculated_templates.txt
  else
    cp "$repo_root/new_template/template/final_list.txt" templates/calculated_templates.txt
  fi
  mkdir -p processed/pdbs
  cp "$STAGED_PDB_DIR"/*.pdb processed/pdbs/
  cp "$run_root/status/available_inputs.csv" inputs.csv
  python_bin="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
  [[ -x "$python_bin" ]] || { echo "missing pipeline interpreter: $python_bin" >&2; exit 127; }
  export PRISM_INPUTS_CSV="$run_root/inputs.csv"
  export PRISM_TM_SCORE_THRESHOLD="$tm_threshold"
  export PRISM_TEMPLATE_WORKFLOW="$template_workflow"
  export PRISM_REFINER="${REFINER:-external_rosetta}"
  export PRISM_STAGE_STATUS_PATH="$run_root/status/stages.jsonl"
  aligner="${ALIGNER:-tmalign}"
  {
    printf 'parameter\tvalue\n'
    printf 'pipeline\t%s\n' "$PIPELINE"
    printf 'template_workflow\t%s\n' "$template_workflow"
    printf 'tm_score_threshold\t%s\n' "$PRISM_TM_SCORE_THRESHOLD"
    printf 'refiner\t%s\n' "$PRISM_REFINER"
    printf 'aligner\t%s\n' "$aligner"
    printf 'alignment_workers\t%s\n' "${ALIGNMENT_WORKERS:-${SLURM_CPUS_PER_TASK:-1}}"
    printf 'usalign_path\t%s\n' "${USALIGN_PATH:-}"
    printf 'usalign_fast\t%s\n' "${USALIGN_FAST:-false}"
    printf 'no_refine\t%s\n' "${NO_REFINE:-0}"
    printf 'gtalign_path\t%s\n' "${GTALIGN_PATH:-}"
    printf 'filter_mode\t%s\n' "${PRISM_FILTER_MODE:-geometry_only_experimental}"
    printf 'filter_asset_root\t%s\n' "${PRISM_FILTER_ASSET_ROOT:-}"
    printf 'template_list\t%s\n' "${TEMPLATE_LIST:-$repo_root/new_template/template/final_list.txt}"
    printf 'template_limit\t%s\n' "${TEMPLATE_LIMIT:-}"
    printf 'minimum_residue_match_count\t%s\n' "${PRISM_MINIMUM_RESIDUE_MATCH_COUNT:-15}"
    printf 'minimum_residue_match_percentage\t%s\n' "${PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE:-50.0}"
    printf 'diff_percentage\t%s\n' "${PRISM_DIFF_PERCENTAGE:-20.0}"
    printf 'clashing_distance\t%s\n' "${PRISM_CLASHING_DISTANCE:-3}"
    printf 'max_clashing_count\t%s\n' "${PRISM_MAX_CLASHING_COUNT:-5}"
    printf 'python\t%s\n' "$python_bin"
    printf 'hostname\t%s\n' "$(hostname)"
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-}"
  } > "$run_root/status/parameters.tsv"
  export PRISM_TMALIGN_WORKERS="${ALIGNMENT_WORKERS:-${SLURM_CPUS_PER_TASK:-1}}"
  usalign_args=()
  if [[ "$aligner" == "usalign" ]]; then
    usalign_args+=(--usalign-path "${USALIGN_PATH:-/home/rshadi25/.conda/envs/gtalign_env/bin/USalign}")
    if [[ "${USALIGN_FAST:-0}" == "1" || "${USALIGN_FAST:-}" == "true" ]]; then
      usalign_args+=(--usalign-fast)
    fi
  fi
  refine_args=()
  if [[ "${NO_REFINE:-0}" == "1" || "${NO_REFINE:-}" == "true" ]]; then
    refine_args+=(--no-refine)
  fi
  if [[ "${NO_REFINE:-0}" != "1" && "${NO_REFINE:-}" != "true" ]]; then
    module load rosetta/2022.42
    export PRISM_ROSETTA_PREPACK="${PRISM_ROSETTA_PREPACK:-/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/docking_prepack_protocol.static.linuxgccrelease}"
    export PRISM_ROSETTA_DOCK="${PRISM_ROSETTA_DOCK:-/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/docking_protocol.static.linuxgccrelease}"
    export PRISM_ROSETTA_DB="${PRISM_ROSETTA_DB:-/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/}"
  fi
  template_args=()
  if [[ -n "${TEMPLATE_LIMIT:-}" ]]; then
    template_args+=(--template-limit "$TEMPLATE_LIMIT")
  fi
  "$python_bin" -u prism.py --aligner "$aligner" --refiner "$PRISM_REFINER" "${template_args[@]}" "${usalign_args[@]}" "${refine_args[@]}" --gtalign_path "${GTALIGN_PATH:-/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu}" > "$run_root/logs/pipeline.out" 2>&1
  printf '{"return_code":0}\n' > "$run_root/status/pipeline_returned.json"
elif [[ "$PIPELINE" == "legacy" ]]; then
  cd "$run_root"
  export PRISM_LOCAL_MODE=1
  [[ -d "$legacy_template_root" ]] || { echo "missing legacy template root: $legacy_template_root" >&2; exit 127; }
  cp -a "$legacy_root/prism.py" "$legacy_root/run_files" "$legacy_root/external_tools" "$legacy_root/prism.ini" .
  if [[ -d "$legacy_template_root" ]]; then
    rm -rf template
    mkdir -p template/interfaces template/contact template/hotspot
    find "$legacy_template_root/interfaces" -maxdepth 1 -type f -name '*.int' -printf '%f\n' \
      | sed -E 's/_[A-Za-z0-9]+\.int$//' > template_ids.all
    if [[ -n "${LEGACY_TEMPLATE_LIMIT:-}" ]]; then
      awk -v limit="$LEGACY_TEMPLATE_LIMIT" '!seen[$0]++ {print; if (++n >= limit) exit}' template_ids.all > template_default
    else
      sort -u template_ids.all > template_default
    fi
    shopt -s nullglob
    while IFS= read -r template_id; do
      for interface_file in "$legacy_template_root/interfaces/${template_id}"_*.int; do
        cp "$interface_file" template/interfaces/
      done
      [[ -f "$legacy_template_root/contact/${template_id}.txt" ]] && cp "$legacy_template_root/contact/${template_id}.txt" template/contact/
      [[ -f "$legacy_template_root/hotspot/hotspot${template_id}" ]] && cp "$legacy_template_root/hotspot/hotspot${template_id}" template/hotspot/
    done < template_default
    shopt -u nullglob
  else
    cp -a "$legacy_root/template" .
    cp -a "$legacy_template_root" template_default
  fi
  mkdir -p inputs jobs
  cp "$run_root/status/available_inputs.csv" /dev/null
  [[ "$template_workflow" == legacy ]] || { echo "legacy pipeline requires TEMPLATE_WORKFLOW=legacy" >&2; exit 2; }
  printf 'parameter\tvalue\n' > "$run_root/status/parameters.tsv"
  printf 'pipeline\tlegacy\ntemplate_workflow\tlegacy\ntm_score_threshold\t%s\n' "$tm_threshold" >> "$run_root/status/parameters.tsv"
  awk -F, 'NR > 1 {print $7, $8}' "$run_root/status/available_inputs.csv" > inputs/pair_list
  cp template_default inputs/template_list
  mkdir -p pdb
  cp "$STAGED_PDB_DIR"/*.pdb pdb/
  job_id="$batch_id"
  mkdir -p "jobs/$job_id/lists"
  cp inputs/pair_list "jobs/$job_id/lists/pair_list"
  cp inputs/template_list "jobs/$job_id/lists/template_list"
  cp prism.ini "jobs/$job_id/prism.ini"
  "$legacy_root/env_py27/bin/python2.7" prism.py "inputs/pair_list" "inputs/template_list" "$job_id" > "$run_root/logs/pipeline.out" 2>&1
else
  echo "unknown PIPELINE=$PIPELINE" >&2
  exit 2
fi

date -Is > "$run_root/status/completed"
