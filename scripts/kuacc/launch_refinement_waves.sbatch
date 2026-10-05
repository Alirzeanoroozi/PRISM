#!/usr/bin/env bash
# Submit bounded KUACC array waves after the manifest setup job completes.
set -u -o pipefail

: "${RUN_ROOT:?set RUN_ROOT}"
: "${SHARD_DIR:?set SHARD_DIR}"
: "${SHARD_MANIFEST:?set SHARD_MANIFEST}"
: "${REFINEMENT_SCRIPT:?set REFINEMENT_SCRIPT}"
: "${WORKER_SCRIPT:?set WORKER_SCRIPT}"
: "${PIPELINE_REPO:?set PIPELINE_REPO}"
: "${FIBERDOCK_DIR:?set FIBERDOCK_DIR}"
: "${DOCKQ_PYTHON:?set DOCKQ_PYTHON}"
: "${PIPELINE_PYTHON:?set PIPELINE_PYTHON}"
: "${EXPECTED_COUNT:?set EXPECTED_COUNT}"

wave_size="${WAVE_SIZE:-8}"
start_file="${RUN_ROOT}/submission/next_shard.txt"
events="${RUN_ROOT}/submission/events.tsv"
mkdir -p "${RUN_ROOT}/submission" "${RUN_ROOT}/logs"

read -r expected shard_count < <("${PIPELINE_PYTHON}" - "${SHARD_MANIFEST}" <<'PY'
import json, sys
d = json.load(open(sys.argv[1], encoding="utf-8"))
print(d["selected_count"], d["shard_count"])
PY
)
if [[ "${expected}" != "${EXPECTED_COUNT}" ]]; then
  echo "manifest_count_mismatch expected=${EXPECTED_COUNT} actual=${expected}" >&2
  exit 1
fi

if [[ -f "${start_file}" ]]; then
  start=$(<"${start_file}")
else
  start=0
fi
if ! [[ "${start}" =~ ^[0-9]+$ ]]; then
  echo "invalid next shard: ${start}" >&2
  exit 1
fi
if (( start >= shard_count )); then
  printf '%s\tcomplete\tnext_shard=%s\n' "$(date -Is)" "${start}" >> "${events}"
  exit 0
fi

end=$((start + wave_size - 1))
(( end >= shard_count )) && end=$((shard_count - 1))
ids=()
submitted_end=$((start - 1))
for shard in $(seq "${start}" "${end}"); do
  rows=$("${PIPELINE_PYTHON}" - "${SHARD_MANIFEST}" "${shard}" <<'PY'
import json, sys
d = json.load(open(sys.argv[1], encoding="utf-8"))
print(d["shards"][int(sys.argv[2])]["rows"])
PY
)
  last=$((rows - 1))
  job_id=$(sbatch --parsable \
    --array="0-${last}" \
    --partition=mid --account=users \
    --cpus-per-task=1 --mem=8G --time=1-00:00:00 \
    --job-name="prism_refine_${shard}" \
    --output="${RUN_ROOT}/logs/refine_%A_%a.out" \
    --error="${RUN_ROOT}/logs/refine_%A_%a.err" \
    --export="ALL,COMPARISON_ROOT=${RUN_ROOT},SHARD_DIR=${SHARD_DIR},SHARD_NUMBER=${shard},PIPELINE_REPO=${PIPELINE_REPO},FIBERDOCK_DIR=${FIBERDOCK_DIR},DOCKQ_PYTHON=${DOCKQ_PYTHON},PIPELINE_PYTHON=${PIPELINE_PYTHON},WORKER_SCRIPT=${WORKER_SCRIPT}" \
    "${REFINEMENT_SCRIPT}" 2>&1)
  rc=$?
  if (( rc != 0 )); then
    printf '%s\tsubmission_failed\tshard=%s\tmessage=%s\n' "$(date -Is)" "${shard}" "${job_id}" >> "${events}"
    exit 2
  fi
  job_id="${job_id%%;*}"
  ids+=("${job_id}")
  submitted_end="${shard}"
  printf '%s\tsubmitted\tshard=%s\tjob_id=%s\tarray=0-%s\n' "$(date -Is)" "${shard}" "${job_id}" "${last}" >> "${events}"
done

next=$((submitted_end + 1))
tmp="${start_file}.tmp.$$"
printf '%s\n' "${next}" > "${tmp}" && mv -f "${tmp}" "${start_file}"

if (( next < shard_count )); then
  dependency=$(IFS=:; echo "${ids[*]}")
  controller=$(sbatch --parsable \
    --dependency="afterany:${dependency}" \
    --partition=mid --account=users --cpus-per-task=1 --mem=1G --time=00:10:00 \
    --job-name=prism_refine_controller \
    --output="${RUN_ROOT}/logs/controller_%j.out" \
    --error="${RUN_ROOT}/logs/controller_%j.err" \
    --export="ALL,RUN_ROOT=${RUN_ROOT},SHARD_DIR=${SHARD_DIR},SHARD_MANIFEST=${SHARD_MANIFEST},REFINEMENT_SCRIPT=${REFINEMENT_SCRIPT},WORKER_SCRIPT=${WORKER_SCRIPT},PIPELINE_REPO=${PIPELINE_REPO},FIBERDOCK_DIR=${FIBERDOCK_DIR},DOCKQ_PYTHON=${DOCKQ_PYTHON},PIPELINE_PYTHON=${PIPELINE_PYTHON},EXPECTED_COUNT=${EXPECTED_COUNT},WAVE_SIZE=${wave_size}" \
    "$0" 2>&1)
  controller="${controller%%;*}"
  if ! [[ "${controller}" =~ ^[0-9]+$ ]]; then
    printf '%s\tcontinuation_submission_failed\tmessage=%s\n' "$(date -Is)" "${controller}" >> "${events}"
    exit 3
  fi
  printf '%s\tcontinuation_scheduled\tcontroller=%s\tafterany=%s\n' "$(date -Is)" "${controller}" "${dependency}" >> "${events}"
else
  printf '%s\tall_waves_submitted\tnext_shard=%s\n' "$(date -Is)" "${next}" >> "${events}"
fi
