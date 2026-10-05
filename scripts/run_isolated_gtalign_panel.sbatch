#!/usr/bin/env bash
# Exact, alignment-only GTalign GPU panel runner.
# All generated files belong under RUN_ROOT; the selected project is read-only.
# The caller supplies --output/--error so Slurm logs never land in the cwd.
#SBATCH --job-name=prism-gtalign-panel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

source /etc/bashrc
set -euo pipefail

template_list="${TEMPLATE_LIST:?set TEMPLATE_LIST}"
template_limit="${TEMPLATE_LIMIT:?set TEMPLATE_LIMIT}"
reference_source="${REFERENCE_SOURCE:?set REFERENCE_SOURCE}"
query_pdb="${QUERY_PDB:?set QUERY_PDB}"
run_root="${RUN_ROOT:?set RUN_ROOT}"
gtalign_bin="${GTALIGN_BIN:-/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_gpu}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"

mkdir -p "$run_root/logs" "$run_root/query" "$run_root/references" "$run_root/output"
[[ -s "$query_pdb" ]] || { echo "missing query PDB: $query_pdb" >&2; exit 2; }
[[ -s "$template_list" ]] || { echo "missing template list: $template_list" >&2; exit 2; }
[[ -x "$gtalign_bin" ]] || { echo "missing executable: $gtalign_bin" >&2; exit 2; }
ln -sfn "$query_pdb" "$run_root/query/$(basename "$query_pdb")"

stage_start=$(date +%s.%N)
"$python_bin" - "$template_list" "$run_root/template_list.txt" "$template_limit" "$reference_source" "$run_root/references" <<'PY'
import hashlib
import json
import os
import pathlib
import sys

source, selected_path, limit, reference_source, reference_target = sys.argv[1:]
limit = int(limit)
lines = [line.strip() for line in pathlib.Path(source).read_text().splitlines() if line.strip()]
source_prefix = lines[:limit]
if len(source_prefix) != limit:
    raise SystemExit(f"requested {limit} templates but source has {len(lines)}")
allow_missing = os.environ.get("ALLOW_MISSING_INTERFACES", "0") == "1"
selected = []
missing = []
for template_id in source_prefix:
    missing_chains = []
    if len(template_id) < 5:
        missing_chains.append("invalid_template_id")
    for chain in template_id[4:]:
        src = pathlib.Path(reference_source) / f"{template_id}_{chain}_int.pdb"
        if not src.is_file() or src.stat().st_size == 0:
            missing_chains.append(chain)
    if missing_chains:
        missing.append({"template_id": template_id, "chains": missing_chains})
        if not allow_missing:
            raise SystemExit(f"missing/empty interface for {template_id}: {missing_chains}")
        continue
    selected.append(template_id)
pathlib.Path(selected_path).write_text("\n".join(selected) + "\n")
target = pathlib.Path(reference_target)
for template_id in selected:
    if len(template_id) < 5:
        raise SystemExit(f"invalid template id: {template_id}")
    for chain in template_id[4:]:
        src = pathlib.Path(reference_source) / f"{template_id}_{chain}_int.pdb"
        dst = target / src.name
        if not src.is_file() or src.stat().st_size == 0:
            raise SystemExit(f"missing/empty interface: {src}")
        dst.symlink_to(src)

def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

pathlib.Path(f"{selected_path}.manifest.json").write_text(json.dumps({
    "source": source,
    "source_count": len(lines),
    "requested_source_prefix_count": len(source_prefix),
    "selected_count": len(selected),
    "selected_sha256": digest(selected_path),
    "interface_count": sum(len(t[4:]) for t in selected),
    "excluded_missing_interfaces": missing,
    "allow_missing_interfaces": allow_missing,
}, indent=2, sort_keys=True) + "\n")
PY
stage_end=$(date +%s.%N)

search_start=$(date +%s.%N)
set +e
"$gtalign_bin" \
  --qrs="$run_root/query" \
  --rfs="$run_root/references" \
  -o "$run_root/output" \
  -s 0 \
  --nhits=1 \
  --nalns=1 \
  --dev-min-length=3 \
  >"$run_root/logs/gtalign.stdout" 2>"$run_root/logs/gtalign.stderr"
return_code=$?
set -e
search_end=$(date +%s.%N)

output_count=$(find "$run_root/output" -maxdepth 1 -type f -name '*.out' | wc -l)
"$python_bin" - "$run_root" "$template_list" "$template_limit" "$query_pdb" "$reference_source" "$gtalign_bin" "$return_code" "$output_count" "$stage_start" "$stage_end" "$search_start" "$search_end" <<'PY'
import json
import pathlib
import socket
import sys

root, template_list, limit, query, references, executable, rc, output_count, stage_start, stage_end, search_start, search_end = sys.argv[1:]
selected_manifest = json.loads(pathlib.Path(root, "template_list.txt.manifest.json").read_text())
data = {
    "schema_version": "prism-gtalign-panel-run-20260913",
    "status": "COMPLETED" if int(rc) == 0 else "FAILED",
    "scientific_status": "alignment_only_no_quality_claim",
    "template_list": template_list,
    "template_limit": int(limit),
    "source_prefix_count": selected_manifest["requested_source_prefix_count"],
    "selected_template_count": selected_manifest["selected_count"],
    "excluded_missing_interfaces": selected_manifest["excluded_missing_interfaces"],
    "selected_template_manifest": f"{root}/template_list.txt.manifest.json",
    "query_pdb": query,
    "reference_source": references,
    "gtalign_executable": executable,
    "gtalign_device": "GPU (requested by GTALIGN_BIN; hardware telemetry not collected)",
    "slurm_job_id": __import__('os').environ.get("SLURM_JOB_ID", ""),
    "hostname": socket.gethostname(),
    "return_code": int(rc),
    "output_count": int(output_count),
    "initialization_or_staging_seconds": float(stage_end) - float(stage_start),
    "alignment_search_seconds": float(search_end) - float(search_start),
    "total_wall_seconds": float(search_end) - float(stage_start),
    "logs": [f"{root}/logs/gtalign.stdout", f"{root}/logs/gtalign.stderr"],
}
pathlib.Path(root, "run_manifest.json").write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY
exit "$return_code"
