#!/bin/bash
# Full target/template leakage matrix. CPU-only and source-gate diagnostic.
# Required exports: SOURCE_GATE_ROOT, INPUTS, STAGED_PDB_DIR, TEMPLATE_LIST,
# INTERFACE_ROOT, SOURCE_POLICY.
#SBATCH --job-name=prism-source-gate
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --qos=cosbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=tmp/agent/20260718-template-source-gate/slurm-%A.out
#SBATCH --error=tmp/agent/20260718-template-source-gate/slurm-%A.err

set -euo pipefail

: "${SOURCE_GATE_ROOT:?set SOURCE_GATE_ROOT}"
: "${INPUTS:?set INPUTS}"
: "${STAGED_PDB_DIR:?set STAGED_PDB_DIR}"
: "${TEMPLATE_LIST:?set TEMPLATE_LIST}"
: "${INTERFACE_ROOT:?set INTERFACE_ROOT}"
: "${SOURCE_POLICY:?set SOURCE_POLICY}"

repo_root="${REPO_ROOT:-/scratch/rshadi25/GitHub/PRISM-prescript}"
case "$SOURCE_GATE_ROOT" in /*) ;; *) SOURCE_GATE_ROOT="$repo_root/$SOURCE_GATE_ROOT" ;; esac
case "$INPUTS" in /*) ;; *) INPUTS="$repo_root/$INPUTS" ;; esac
case "$STAGED_PDB_DIR" in /*) ;; *) STAGED_PDB_DIR="$repo_root/$STAGED_PDB_DIR" ;; esac
case "$TEMPLATE_LIST" in /*) ;; *) TEMPLATE_LIST="$repo_root/$TEMPLATE_LIST" ;; esac
case "$INTERFACE_ROOT" in /*) ;; *) INTERFACE_ROOT="$repo_root/$INTERFACE_ROOT" ;; esac
case "$SOURCE_POLICY" in /*) ;; *) SOURCE_POLICY="$repo_root/$SOURCE_POLICY" ;; esac

mkdir -p "$SOURCE_GATE_ROOT"
python_bin="${PIPELINE_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
printf 'parameter\tvalue\n' > "$SOURCE_GATE_ROOT/parameters.tsv"
printf 'hostname\t%s\n' "$(hostname)" >> "$SOURCE_GATE_ROOT/parameters.tsv"
printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-}" >> "$SOURCE_GATE_ROOT/parameters.tsv"
printf 'cpus_per_task\t%s\n' "${SLURM_CPUS_PER_TASK:-}" >> "$SOURCE_GATE_ROOT/parameters.tsv"
printf 'template_list\t%s\n' "$TEMPLATE_LIST" >> "$SOURCE_GATE_ROOT/parameters.tsv"
printf 'template_list_sha256\t%s\n' "$(sha256sum "$TEMPLATE_LIST" | awk '{print $1}')" >> "$SOURCE_GATE_ROOT/parameters.tsv"
printf 'started\t%s\n' "$(date -Is)" >> "$SOURCE_GATE_ROOT/parameters.tsv"

"$python_bin" "$repo_root/benchmark/scripts/run_template_source_gate.py" \
  --inputs "$INPUTS" \
  --staged-pdb-dir "$STAGED_PDB_DIR" \
  --template-list "$TEMPLATE_LIST" \
  --interface-root "$INTERFACE_ROOT" \
  --source-policy "$SOURCE_POLICY" \
  --similarity-output "$SOURCE_GATE_ROOT/target_template_similarity.tsv" \
  --eligible-output "$SOURCE_GATE_ROOT/eligible_templates_by_pair.tsv" \
  > "$SOURCE_GATE_ROOT/source_gate_summary.json"

printf 'completed\t%s\n' "$(date -Is)" >> "$SOURCE_GATE_ROOT/parameters.tsv"
