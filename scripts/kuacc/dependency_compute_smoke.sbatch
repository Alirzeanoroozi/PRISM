#!/usr/bin/env bash
# KUACC compute-node dependency smoke test; no production workload.
#SBATCH --job-name=prism-dependency-smoke
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/users/rshadi25/valar-remote-runs/prism-prescript-kuacc-preflight-20260927/logs/compute-smoke-%j.out
#SBATCH --error=/scratch/users/rshadi25/valar-remote-runs/prism-prescript-kuacc-preflight-20260927/logs/compute-smoke-%j.err

source /etc/bashrc
set -euo pipefail

ROOT=${PRISM_PREFLIGHT_ROOT:?set PRISM_PREFLIGHT_ROOT}
TOOLS="$ROOT/tools"
FIBER="$TOOLS/fiberdock"
LOGS="$ROOT/logs"
REPORT="$ROOT/manifest/compute_smoke_${SLURM_JOB_ID}.tsv"
mkdir -p "$LOGS" "$ROOT/manifest"
printf 'step\tstatus\telapsed_seconds\tdetail\n' > "$REPORT"

run_usage() {
  local name=$1
  local cwd=$2
  local executable=$3
  local pattern=$4
  shift 4
  local start end rc
  start=$(date +%s.%N)
  set +e
  (cd "$cwd" && timeout 20s "$executable" "$@") \
    > "$LOGS/compute_${SLURM_JOB_ID}_${name}.stdout" \
    2> "$LOGS/compute_${SLURM_JOB_ID}_${name}.stderr"
  rc=$?
  set -e
  end=$(date +%s.%N)
  if grep -Eiq "$pattern" "$LOGS/compute_${SLURM_JOB_ID}_${name}.stdout" "$LOGS/compute_${SLURM_JOB_ID}_${name}.stderr"; then
    printf '%s\tPASS\t%s\trc=%s\n' "$name" "$(awk -v s="$start" -v e="$end" 'BEGIN { print e-s }')" "$rc" >> "$REPORT"
  else
    printf '%s\tFAIL\t%s\trc=%s\n' "$name" "$(awk -v s="$start" -v e="$end" 'BEGIN { print e-s }')" "$rc" >> "$REPORT"
    return 1
  fi
}

mkdir -p "$LOGS"
printf 'host=%s\njob_id=%s\n' "$(hostname)" "$SLURM_JOB_ID" > "$ROOT/manifest/compute_smoke_${SLURM_JOB_ID}.meta"
module load rosetta/2022.42
rosetta_prepack=$(command -v docking_prepack_protocol.static.linuxgccrelease)
rosetta_dock=$(command -v docking_protocol.static.linuxgccrelease)
rosetta_db=${PRISM_ROSETTA_DB:-/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database}
[[ -d "$rosetta_db" ]]
[[ -x "$rosetta_prepack" && -x "$rosetta_dock" ]]

run_usage multiprot_usage "$ROOT" "$TOOLS/multiprot.Linux" 'Usage:' --help || exit 1
run_usage rosetta_prepack_help "$ROOT" "$rosetta_prepack" 'Rosetta|Usage|options' -help || exit 1
run_usage rosetta_dock_help "$ROOT" "$rosetta_dock" 'Rosetta|Usage|options' -help || exit 1
run_usage fiberdock_usage "$FIBER" "$FIBER/FiberDock" 'Usage:' || exit 1
run_usage nma_usage "$FIBER" "$FIBER/nma" 'Usage:|input.pdb' || exit 1
run_usage reduce_usage "$FIBER" "$FIBER/reduce" 'Usage:|arguments:' || exit 1

sha256sum "$TOOLS/multiprot.Linux" "$FIBER/FiberDock" "$FIBER/nma" "$FIBER/reduce" > "$ROOT/manifest/compute_smoke_${SLURM_JOB_ID}.sha256"
printf 'status=PASS\ncompleted_at=%s\n' "$(date -Is)" >> "$ROOT/manifest/compute_smoke_${SLURM_JOB_ID}.meta"
