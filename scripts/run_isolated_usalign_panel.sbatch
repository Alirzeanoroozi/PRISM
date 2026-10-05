#!/usr/bin/env bash
# Exact alignment-only USalign panel runner. Every materialized pair receives
# one JSONL result, including failed/unparseable outcomes. No canonical file is
# written. Slurm output paths are supplied by the caller.
#SBATCH --job-name=prism-usalign-panel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00

source /etc/bashrc
set -euo pipefail

template_list="${TEMPLATE_LIST:?set TEMPLATE_LIST}"
template_limit="${TEMPLATE_LIMIT:?set TEMPLATE_LIMIT}"
reference_source="${REFERENCE_SOURCE:?set REFERENCE_SOURCE}"
query_pdb="${QUERY_PDB:?set QUERY_PDB}"
run_root="${RUN_ROOT:?set RUN_ROOT}"
usalign_bin="${USALIGN_BIN:-/home/rshadi25/.conda/envs/gtalign_env/bin/USalign}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
workers="${USALIGN_WORKERS:-8}"

mkdir -p "$run_root/logs" "$run_root/results"
[[ -s "$query_pdb" ]] || { echo "missing query PDB: $query_pdb" >&2; exit 2; }
[[ -s "$template_list" ]] || { echo "missing template list: $template_list" >&2; exit 2; }
[[ -x "$usalign_bin" ]] || { echo "missing executable: $usalign_bin" >&2; exit 2; }

stage_start=$(date +%s.%N)
"$python_bin" - "$template_list" "$template_limit" "$reference_source" "$query_pdb" "$run_root/pairs.tsv" "$run_root/panel_manifest.json" <<'PY'
import hashlib
import json
import pathlib
import sys

source, limit, ref_root, query, pair_path, manifest_path = sys.argv[1:]
limit = int(limit)
lines = [line.strip() for line in pathlib.Path(source).read_text().splitlines() if line.strip()]
prefix = lines[:limit]
if len(prefix) != limit:
    raise SystemExit(f"requested {limit} entries but source has {len(lines)}")
selected, missing = [], []
for template_id in prefix:
    missing_chains = []
    for chain in template_id[4:]:
        p = pathlib.Path(ref_root) / f"{template_id}_{chain}_int.pdb"
        if not p.is_file() or p.stat().st_size == 0:
            missing_chains.append(chain)
    if missing_chains:
        missing.append({"template_id": template_id, "chains": missing_chains})
    else:
        selected.append(template_id)

def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

pairs = []
with pathlib.Path(pair_path).open('w') as out:
    for template_id in selected:
        for chain in template_id[4:]:
            ref = pathlib.Path(ref_root) / f"{template_id}_{chain}_int.pdb"
            pairs.append((template_id, chain, str(ref)))
            out.write(f"{query}\t{template_id}\t{chain}\t{ref}\n")
pathlib.Path(manifest_path).write_text(json.dumps({
    "schema_version": "prism-usalign-panel-manifest-20260913",
    "source": source, "source_count": len(lines),
    "requested_source_prefix_count": len(prefix), "selected_template_count": len(selected),
    "selected_interface_pair_count": len(pairs),
    "selected_pair_manifest_sha256": digest(pair_path),
    "selected_template_ids_sha256": hashlib.sha256(("\n".join(selected) + "\n").encode()).hexdigest(),
    "excluded_missing_interfaces": missing,
    "invocation": ["USalign", "QUERY", "REFERENCE", "-outfmt", "-1", "-m", "-"],
}, indent=2, sort_keys=True) + "\n")
PY
stage_end=$(date +%s.%N)

search_start=$(date +%s.%N)
set +e
"$python_bin" - "$run_root/pairs.tsv" "$run_root/results/alignments.jsonl" "$usalign_bin" "$workers" <<'PY'
import hashlib
import json
import multiprocessing as mp
import pathlib
import re
import subprocess
import sys
import time

pairs_path, result_path, executable, worker_count = sys.argv[1:]
pairs = []
for line in pathlib.Path(pairs_path).read_text().splitlines():
    query, template_id, chain, ref = line.split('\t')
    pairs.append((query, template_id, chain, ref))

def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()

def run(pair):
    query, template_id, chain, ref = pair
    started = time.perf_counter()
    command = [executable, query, ref, '-outfmt', '-1', '-m', '-']
    try:
        cp = subprocess.run(command, capture_output=True, text=True, timeout=600)
        stdout = cp.stdout or ''
        aligned = re.search(r'Aligned length\s*=\s*(\d+)', stdout)
        tm = re.search(r'TM-score\s*=\s*([0-9.]+)', stdout)
        transform_rows = sum(1 for line in stdout.splitlines() if len(line.split()) == 5 and line.split()[0] in {'0', '1', '2'})
        ok = cp.returncode == 0 and aligned is not None and tm is not None and transform_rows >= 3
        status = 'success' if ok else 'unparseable' if cp.returncode == 0 else 'failed'
        return {"query": query, "template_id": template_id, "chain": chain,
                "status": status, "return_code": cp.returncode,
                "aligned_length": int(aligned.group(1)) if aligned else None,
                "tm_score": float(tm.group(1)) if tm else None,
                "transform_rows": transform_rows, "stdout_sha256": digest(stdout),
                "elapsed_seconds": time.perf_counter() - started,
                "failure_reason": '' if ok else 'missing required alignment/transform fields'}
    except subprocess.TimeoutExpired:
        return {"query": query, "template_id": template_id, "chain": chain,
                "status": 'failed', "return_code": None, "aligned_length": None,
                "tm_score": None, "transform_rows": 0, "stdout_sha256": '',
                "elapsed_seconds": time.perf_counter() - started,
                "failure_reason": 'timeout_after_600_seconds'}

with pathlib.Path(result_path).open('w') as out:
    with mp.Pool(processes=int(worker_count)) as pool:
        for result in pool.imap(run, pairs, chunksize=1):
            out.write(json.dumps(result, sort_keys=True) + '\n')
            out.flush()
PY
return_code=$?
set -e
search_end=$(date +%s.%N)

"$python_bin" - "$run_root" "$template_list" "$template_limit" "$usalign_bin" "$workers" "$return_code" "$stage_start" "$stage_end" "$search_start" "$search_end" <<'PY'
import json
import pathlib
import socket
import sys
root, source, limit, executable, workers, rc, stage_start, stage_end, search_start, search_end = sys.argv[1:]
manifest = json.loads(pathlib.Path(root, 'panel_manifest.json').read_text())
results = [json.loads(line) for line in pathlib.Path(root, 'results/alignments.jsonl').read_text().splitlines() if line.strip()]
counts = {}
for row in results:
    counts[row['status']] = counts.get(row['status'], 0) + 1
data = {
    "schema_version": "prism-usalign-panel-run-20260913",
    "status": "COMPLETED" if int(rc) == 0 and len(results) == manifest['selected_interface_pair_count'] else "FAILED",
    "scientific_status": "alignment_only_no_quality_claim",
    "template_list": source, "template_limit": int(limit),
    "selected_template_count": manifest['selected_template_count'],
    "selected_interface_pair_count": manifest['selected_interface_pair_count'],
    "excluded_missing_interfaces": manifest['excluded_missing_interfaces'],
    "query_pdb": "recorded in pairs.tsv", "usalign_executable": executable,
    "invocation": manifest['invocation'], "usalign_version": "20241108",
    "worker_count": int(workers), "slurm_job_id": __import__('os').environ.get('SLURM_JOB_ID', ''),
    "hostname": socket.gethostname(), "return_code": int(rc),
    "result_count": len(results), "result_status_counts": counts,
    "initialization_or_staging_seconds": float(stage_end) - float(stage_start),
    "alignment_search_seconds": float(search_end) - float(search_start),
    "total_wall_seconds": float(search_end) - float(stage_start),
    "panel_manifest": f"{root}/panel_manifest.json", "results": f"{root}/results/alignments.jsonl",
}
pathlib.Path(root, 'run_manifest.json').write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
PY
exit "$return_code"
