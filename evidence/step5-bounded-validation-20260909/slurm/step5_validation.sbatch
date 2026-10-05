#!/usr/bin/env bash
# Bounded Step 5 validation. Generated outputs stay in /scratch/tmp.
#SBATCH --partition=kutem
#SBATCH --account=kutem
#SBATCH --qos=kutem
#SBATCH --job-name=prism-step5
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=/home/rshadi25/valar-agent-framework/evidence/prism-prescript/step5-bounded-validation-20260909/logs-slurm-%j.out
#SBATCH --error=/home/rshadi25/valar-agent-framework/evidence/prism-prescript/step5-bounded-validation-20260909/logs-slurm-%j.err

set -euo pipefail
source /etc/bashrc 2>/dev/null || true

PYTHON_BIN="/home/rshadi25/.conda/envs/gtalign_env/bin/python"
PROJECT="/scratch/rshadi25/GitHub/PRISM-prescript"
CANON="$PROJECT"
USALIGN_WT="$PROJECT/tmp/agent/worktrees/run-1e6873e93a4f4082a236f7348218ecd5"
PRODIGY_WT="$PROJECT/tmp/agent/worktrees/run-4f0e2099ca254fde8dbd32b71a96e740"
RUN_DIR="/home/rshadi25/valar-agent-framework/evidence/prism-prescript/step5-bounded-validation-20260909"
ROOT="/scratch/tmp/prism-step5-validation-${SLURM_JOB_ID}"
mkdir -p "$ROOT" "$RUN_DIR/evidence"

printf 'job_id=%s\nhostname=%s\npython=%s\nroot=%s\n' \
  "$SLURM_JOB_ID" "$(hostname)" "$PYTHON_BIN" "$ROOT" \
  > "$RUN_DIR/evidence/runtime-${SLURM_JOB_ID}.txt"

stage_workspace() {
  local name="$1"
  local prism_source="$2"
  local workspace="$ROOT/$name"
  mkdir -p "$workspace/src" "$workspace/templates" "$workspace/processed/pdbs"
  cp -a "$CANON/src/." "$workspace/src/"
  cp -a "$prism_source" "$workspace/prism.py"
  ln -s "$CANON/external_tools" "$workspace/external_tools"
  ln -s "$CANON/templates/interfaces" "$workspace/templates/interfaces"
  ln -s "$CANON/templates/interfaces_lists" "$workspace/templates/interfaces_lists"
  ln -s "$CANON/processed/pdbs/1fgnHL.pdb" "$workspace/processed/pdbs/1fgn.pdb"
  ln -s "$CANON/processed/pdbs/1tfhA.pdb" "$workspace/processed/pdbs/1tfh.pdb"
  cat > "$workspace/inputs.csv" <<'EOF'
Receptor,Ligand
1FGNHL,1TFHA
EOF
  printf '%s\n' '1kcaCH' > "$workspace/templates/calculated_templates.txt"
  printf '%s\n' "$workspace"
}

run_pipeline_arm() {
  local name="$1"
  local prism_source="$2"
  local extra_args="$3"
  local workspace
  workspace="$(stage_workspace "$name" "$prism_source")"
  set +e
  (
    cd "$workspace"
    export PRISM_NACCESS_EXECUTABLE="external_tools/naccess/naccess"
    export PRISM_TMALIGN="external_tools/TMalign"
    export PRISM_TMALIGN_WORKERS=1
    timeout 600 "$PYTHON_BIN" prism.py \
      --inputs-csv inputs.csv --template-limit 1 --no-refine \
      --surface-backend naccess --rank false $extra_args
  ) > "$workspace/run.log" 2>&1
  local rc=$?
  set -e
  printf '%s\n' "$rc" > "$workspace/return_code"
  find "$workspace/processed" -type f -name '*.json' -printf '%P\n' | sort > "$workspace/json_files.txt" || true
  find "$workspace/processed" -type f -name '*.pdb' -printf '%P\n' | sort > "$workspace/pdb_files.txt" || true
}

run_pipeline_arm canonical "$CANON/prism.py" "--aligner tmalign"
run_pipeline_arm usalign "$USALIGN_WT/prism.py" "--aligner usalign --usalign-path /home/rshadi25/.conda/envs/gtalign_env/bin/USalign"

PRODIGY_ROOT="$ROOT/prodigy"
mkdir -p "$PRODIGY_ROOT/processed/ranking/prodigy"
LEFT="$PRODIGY_ROOT/1a0cCD_5zngA_4eylA_o1_L.pdb"
RIGHT="$PRODIGY_ROOT/1a0cCD_5zngA_4eylA_o1_R.pdb"
ln -s /scratch/tmp/prism-step4-pair-fixture-20260909/left.pdb "$LEFT"
ln -s /scratch/tmp/prism-step4-pair-fixture-20260909/right.pdb "$RIGHT"
set +e
(
  cd "$PRODIGY_ROOT"
  export PATH="/home/rshadi25/.conda/envs/gtalign_env/bin:$PATH"
  PYTHONPATH="$PRODIGY_WT" "$PYTHON_BIN" \
    "$RUN_DIR/validation/probe_prodigy.py" "$LEFT" "$RIGHT" \
    "$PRODIGY_ROOT/processed/ranking/prodigy"
) > "$PRODIGY_ROOT/run.log" 2>&1
PRODIGY_RC=$?
set -e
printf '%s\n' "$PRODIGY_RC" > "$PRODIGY_ROOT/return_code"

for arm in canonical usalign prodigy; do
  printf '%s\t%s\t%s\t%s\n' \
    "$arm" \
    "$(cat "$ROOT/$arm/return_code")" \
    "$(find "$ROOT/$arm/processed" -type f -name '*.json' 2>/dev/null | wc -l)" \
    "$(find "$ROOT/$arm/processed" -type f -name '*.pdb' 2>/dev/null | wc -l)"
done > "$RUN_DIR/evidence/summary-${SLURM_JOB_ID}.tsv"
printf 'root=%s\n' "$ROOT" > "$RUN_DIR/evidence/root-${SLURM_JOB_ID}.txt"
