#!/usr/bin/env bash
# Isolated execution template. Submit only through an explicit command such as:
#   python3 benchmark/scripts/isolated_kutem_runner.py --manifest <manifest.csv> \
#     --run-root <run-root> --submit
# Local validation uses --dry-run and never invokes sbatch.
#SBATCH --job-name=array-kutem
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --array=1-10
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

set -euo pipefail

: "${MANIFEST:?sbatch export must define MANIFEST=<task-manifest.csv>}"
: "${RUN_ROOT:?sbatch export must define RUN_ROOT=<isolated-run-root>}"

mkdir -p "$RUN_ROOT/launcher_failures"
launcher_failure="$RUN_ROOT/launcher_failures/array_${SLURM_ARRAY_TASK_ID:-unknown}.json"
write_launcher_failure() {
  local rc=$?
  if [[ "$rc" -ne 0 ]]; then
    printf '{"schema_version":"isolated-kutem-launcher-failure/v1","return_code":%s,"job_id":"%s","array_task_id":"%s"}\n' \
      "$rc" "${SLURM_JOB_ID:-}" "${SLURM_ARRAY_TASK_ID:-}" > "$launcher_failure"
  fi
  return "$rc"
}
trap write_launcher_failure EXIT

repo_root="${PRISM_REPO_ROOT:-$(git -C "$(dirname "$MANIFEST")" rev-parse --show-toplevel 2>/dev/null || true)}"
: "${repo_root:?unable to resolve PRISM_REPO_ROOT; export it explicitly}"
python_bin="${PRISM_RUNNER_PYTHON:-python3}"

set +e
"$python_bin" "$repo_root/benchmark/scripts/isolated_kutem_runner.py" \
  --manifest "$MANIFEST" \
  --run-root "$RUN_ROOT" \
  --run-task
runner_status=$?
exit "$runner_status"
