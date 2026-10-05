#!/usr/bin/env bash
# Exact alignment-only TMalign panel runner.
# Every materialized pair receives a JSONL result, including failed or
# unparseable outcomes. The canonical project and its interfaces are read-only.
# Slurm output paths are supplied by the caller.
#SBATCH --job-name=prism-tmalign-panel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00

source /etc/bashrc
set -euo pipefail

template_list="${TEMPLATE_LIST:?set TEMPLATE_LIST}"
template_limit="${TEMPLATE_LIMIT:?set TEMPLATE_LIMIT}"
reference_source="${REFERENCE_SOURCE:?set REFERENCE_SOURCE}"
query_pdb="${QUERY_PDB:?set QUERY_PDB}"
run_root="${RUN_ROOT:?set RUN_ROOT}"
tmalign_bin="${TMALIGN_BIN:-/home/rshadi25/.conda/envs/gtalign_env/bin/TMalign}"
python_bin="${RUN_PYTHON:-/home/rshadi25/.conda/envs/gtalign_env/bin/python}"
workers="${TMALIGN_WORKERS:-8}"

mkdir -p "$run_root/logs" "$run_root/results"
[[ -s "$query_pdb" ]] || { echo "missing query PDB: $query_pdb" >&2; exit 2; }
[[ -s "$template_list" ]] || { echo "missing template list: $template_list" >&2; exit 2; }
[[ -x "$tmalign_bin" ]] || { echo "missing executable: $tmalign_bin" >&2; exit 2; }

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
    if len(template_id) < 5:
        missing_chains.append("invalid_template_id")
    for chain in template_id[4:]:
        path = pathlib.Path(ref_root) / f"{template_id}_{chain}_int.pdb"
        if not path.is_file() or path.stat().st_size == 0:
            missing_chains.append(chain)
    if missing_chains:
        missing.append({"template_id": template_id, "chains": missing_chains})
    else:
        selected.append(template_id)

def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

pairs = []
with pathlib.Path(pair_path).open("w") as out:
    for template_id in selected:
        for chain in template_id[4:]:
            reference = pathlib.Path(ref_root) / f"{template_id}_{chain}_int.pdb"
            pairs.append((template_id, chain, str(reference)))
            out.write(f"{query}\t{template_id}\t{chain}\t{reference}\n")
pathlib.Path(manifest_path).write_text(json.dumps({
    "schema_version": "prism-tmalign-panel-manifest-20260913",
    "source": source,
    "source_count": len(lines),
    "requested_source_prefix_count": len(prefix),
    "selected_template_count": len(selected),
    "selected_interface_pair_count": len(pairs),
    "selected_pair_manifest_sha256": digest(pair_path),
    "selected_template_ids_sha256": hashlib.sha256(("\n".join(selected) + "\n").encode()).hexdigest(),
    "excluded_missing_interfaces": missing,
    "invocation": ["TMalign", "QUERY", "REFERENCE", "-m", "MATRIX_PATH"],
}, indent=2, sort_keys=True) + "\n")
PY
stage_end=$(date +%s.%N)

search_start=$(date +%s.%N)
set +e
"$python_bin" - "$run_root/pairs.tsv" "$run_root/results/alignments.jsonl" "$run_root/tmp" "$tmalign_bin" "$workers" <<'PY'
import hashlib
import json
import multiprocessing as mp
import pathlib
import re
import subprocess
import sys
import tempfile
import time

pairs_path, result_path, tmp_root, executable, worker_count = sys.argv[1:]
pathlib.Path(tmp_root).mkdir(parents=True, exist_ok=True)
pairs = []
for line in pathlib.Path(pairs_path).read_text().splitlines():
    query, template_id, chain, reference = line.split("\t")
    pairs.append((query, template_id, chain, reference))

def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()

def run(pair):
    query, template_id, chain, reference = pair
    started = time.perf_counter()
    try:
        with tempfile.TemporaryDirectory(prefix="tmalign-", dir=tmp_root) as scratch:
            matrix_path = pathlib.Path(scratch) / "matrix.out"
            command = [executable, query, reference, "-m", str(matrix_path)]
            completed = subprocess.run(command, capture_output=True, text=True, timeout=600)
            stdout = completed.stdout or ""
            matrix = matrix_path.read_bytes() if matrix_path.is_file() else b""
            matrix_rows = set()
            for line in matrix.decode(errors="replace").splitlines():
                fields = line.split()
                if len(fields) >= 5 and fields[0] in {"0", "1", "2"}:
                    matrix_rows.add(fields[0])
            aligned = re.search(r"Aligned length\s*=\s*(\d+)", stdout)
            scores = [float(value) for value in re.findall(r"TM-score\s*=\s*([0-9.]+)", stdout)]
            marker = any(line.startswith('(":') for line in stdout.splitlines())
            ok = completed.returncode == 0 and aligned is not None and scores and matrix_rows == {"0", "1", "2"} and marker
            status = "success" if ok else "unparseable" if completed.returncode == 0 else "failed"
            return {
                "query": query,
                "template_id": template_id,
                "chain": chain,
                "status": status,
                "return_code": completed.returncode,
                "aligned_length": int(aligned.group(1)) if aligned else None,
                "tm_score": max(scores) if scores else None,
                "transform_rows": len(matrix_rows),
                "alignment_marker": marker,
                "stdout_sha256": digest_bytes(stdout.encode()),
                "matrix_sha256": digest_bytes(matrix),
                "elapsed_seconds": time.perf_counter() - started,
                "failure_reason": "" if ok else "missing required alignment or matrix fields",
            }
    except subprocess.TimeoutExpired:
        return {
            "query": query, "template_id": template_id, "chain": chain,
            "status": "failed", "return_code": None, "aligned_length": None,
            "tm_score": None, "transform_rows": 0, "alignment_marker": False,
            "stdout_sha256": "", "matrix_sha256": "",
            "elapsed_seconds": time.perf_counter() - started,
            "failure_reason": "timeout_after_600_seconds",
        }
    except OSError as exc:
        return {
            "query": query, "template_id": template_id, "chain": chain,
            "status": "failed", "return_code": None, "aligned_length": None,
            "tm_score": None, "transform_rows": 0, "alignment_marker": False,
            "stdout_sha256": "", "matrix_sha256": "",
            "elapsed_seconds": time.perf_counter() - started,
            "failure_reason": f"os_error: {exc}",
        }

with pathlib.Path(result_path).open("w") as out:
    with mp.Pool(processes=int(worker_count)) as pool:
        for result in pool.imap(run, pairs, chunksize=1):
            out.write(json.dumps(result, sort_keys=True) + "\n")
            out.flush()
PY
return_code=$?
set -e
search_end=$(date +%s.%N)

"$python_bin" - "$run_root" "$template_list" "$template_limit" "$tmalign_bin" "$workers" "$return_code" "$stage_start" "$stage_end" "$search_start" "$search_end" <<'PY'
import json
import pathlib
import socket
import sys

root, source, limit, executable, workers, rc, stage_start, stage_end, search_start, search_end = sys.argv[1:]
panel = json.loads(pathlib.Path(root, "panel_manifest.json").read_text())
results_path = pathlib.Path(root, "results/alignments.jsonl")
results = [json.loads(line) for line in results_path.read_text().splitlines() if line.strip()]
counts = {}
for row in results:
    counts[row["status"]] = counts.get(row["status"], 0) + 1
complete = int(rc) == 0 and len(results) == panel["selected_interface_pair_count"] and counts.get("success", 0) == len(results)
data = {
    "schema_version": "prism-tmalign-panel-run-20260913",
    "status": "COMPLETED" if complete else "FAILED",
    "scientific_status": "alignment_only_no_quality_claim",
    "template_list": source,
    "template_limit": int(limit),
    "selected_template_count": panel["selected_template_count"],
    "selected_interface_pair_count": panel["selected_interface_pair_count"],
    "excluded_missing_interfaces": panel["excluded_missing_interfaces"],
    "query_pdb": "recorded in pairs.tsv",
    "tmalign_executable": executable,
    "tmalign_version": "20220412",
    "worker_count": int(workers),
    "slurm_job_id": __import__("os").environ.get("SLURM_JOB_ID", ""),
    "hostname": socket.gethostname(),
    "return_code": int(rc),
    "result_count": len(results),
    "result_status_counts": counts,
    "initialization_or_staging_seconds": float(stage_end) - float(stage_start),
    "alignment_search_seconds": float(search_end) - float(search_start),
    "total_wall_seconds": float(search_end) - float(stage_start),
    "panel_manifest": f"{root}/panel_manifest.json",
    "results": f"{root}/results/alignments.jsonl",
}
pathlib.Path(root, "run_manifest.json").write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY
exit "$return_code"
