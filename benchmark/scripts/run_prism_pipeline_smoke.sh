#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
default_python=""
for candidate in \
  "${PRISM_PIPELINE_PYTHON:-}" \
  "/home/rshadi25/.conda/envs/gtalign_env/bin/python3.11" \
  "${PRISM_TEST_PYTHON:-}" \
  "$repo_dir/benchmark/prism_processed/env/prism_score_env/bin/python" \
  "python3"
do
  if [[ -z "$candidate" ]]; then
    continue
  fi
  if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c "import pandas, numpy, Bio" >/dev/null 2>&1; then
    default_python=$candidate
    break
  fi
done
python_bin=${1:-$default_python}
work_root=${2:-$(mktemp -d "$repo_dir/tmp/agent/prism-pipeline-smoke-XXXXXX")}
template_id=${PRISM_SMOKE_TEMPLATE_ID:-1kcaCH}
receptor_id=${PRISM_SMOKE_RECEPTOR_ID:-1FGNH}
ligand_id=${PRISM_SMOKE_LIGAND_ID:-1TFHA}
timeout_sec=${PRISM_SMOKE_TIMEOUT_SEC:-300}
pipeline_args=()
if [[ -n "${PRISM_SMOKE_ALIGNER:-}" ]]; then
  pipeline_args+=(--aligner "${PRISM_SMOKE_ALIGNER}")
fi
if [[ "${PRISM_SMOKE_NO_REFINE:-0}" == "1" ]]; then
  pipeline_args+=(--no-refine)
fi
if [[ -n "${PRISM_SMOKE_GTALIGN_PATH:-}" ]]; then
  pipeline_args+=(--gtalign-path "${PRISM_SMOKE_GTALIGN_PATH}")
fi

if [[ -z "$python_bin" ]]; then
  echo "No compatible Python interpreter found for the PRISM pipeline smoke test." >&2
  exit 1
fi

if ! command -v "$python_bin" >/dev/null 2>&1; then
  echo "Python interpreter not found or not executable: $python_bin" >&2
  exit 1
fi

mkdir -p "$work_root"/templates "$work_root"/processed "$work_root"/processed/rosetta_refinement
mkdir -p "$work_root"/status
ln -sfn "$repo_dir/src" "$work_root/src"
ln -sfn "$repo_dir/external_tools" "$work_root/external_tools"
ln -sfn "$repo_dir/prism.py" "$work_root/prism.py"
ln -sfn "$repo_dir/processed/pdbs" "$work_root/processed/pdbs"
ln -sfn "$repo_dir/templates/interfaces" "$work_root/templates/interfaces"
ln -sfn "$repo_dir/templates/interfaces_lists" "$work_root/templates/interfaces_lists"

cat >"$work_root/inputs.csv" <<EOF
Receptor,Ligand
$receptor_id,$ligand_id
EOF

printf '%s\n' "$template_id" >"$work_root/templates/calculated_templates.txt"

echo "Smoke workspace: $work_root"
echo "Input pair: $receptor_id vs $ligand_id"
echo "Template: $template_id"
echo "Interpreter: $python_bin"

set +e
(
  cd "$work_root"
  # Omit the template-generation flag: prism.py defaults to False, while
  # argparse's type=bool would incorrectly treat the string "false" as True.
  timeout "$timeout_sec" "$python_bin" prism.py "${pipeline_args[@]}"
) | tee "$work_root/run.log"
pipeline_status=${PIPESTATUS[0]}
set -e

echo
echo "Pipeline exit status: $pipeline_status"
echo "Log: $work_root/run.log"
echo "Alignment dir: $work_root/processed/alignment"
echo "Transformation dir: $work_root/processed/transformation"
echo "Rosetta dir: $work_root/processed/rosetta_refinement"

printf '{"return_code":%s,"work_root":"%s"}\n' \
  "$pipeline_status" "$work_root" >"$work_root/status/pipeline_returned.json"

if [[ $pipeline_status -eq 124 ]]; then
  echo "Smoke run timed out." >&2
  exit 124
fi

if [[ $pipeline_status -ne 0 ]]; then
  exit "$pipeline_status"
fi
