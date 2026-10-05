
#!/usr/bin/env python3
"""Run a corrected FiberDock output-contract replay in isolated workspaces.

The input parameter files and source structures are read-only evidence. Only the
declared energiesOutFileName line is rewritten into the fresh output workspace.
FiberDock binaries are copied per candidate so the shared external_tools tree is
never mutated.
"""

from __future__ import annotations

import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.scripts.diagnose_fiberdock_output_contract import (  # noqa: E402
    parse_energy,
    pdb_integrity,
)

CANDIDATES = (
    (
        "multiprot_1ahwBC_o1",
        Path("tmp/agent/20260728-matched-align-refiner-replay/results-v2/multiprot"),
        Path("processed/fiberdock_refinement/1ahwBC_1fgnHL_1tfhA_o1/fd_params.txt"),
    ),
    (
        "tmalign_1ahwAF_o1",
        Path("tmp/agent/20260728-matched-align-refiner-replay/results-v2/tmalign"),
        Path("processed/fiberdock_refinement/1ahwAF_1fgnHL_1tfhA_o1/fd_params.txt"),
    ),
)


def rewrite_declared_stem(source_params: Path, target_params: Path, stem: Path) -> None:
    lines = []
    replaced = False
    for line in source_params.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("energiesOutFileName "):
            lines.append(f"energiesOutFileName {stem}")
            replaced = True
        else:
            lines.append(line)
    if not replaced:
        raise RuntimeError(f"missing energiesOutFileName in {source_params}")
    target_params.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_candidate(label: str, source_arm: Path, relative_params: Path, output_root: Path) -> dict:
    source_arm = source_arm.resolve()
    source_params = source_arm / relative_params
    candidate_root = output_root / label
    candidate_root.mkdir(parents=True, exist_ok=True)
    tools_root = candidate_root / "external_tools" / "fiberdock"
    shutil.copytree(
        source_arm / "external_tools" / "fiberdock",
        tools_root,
        symlinks=True,
        dirs_exist_ok=True,
    )
    params_path = candidate_root / "fd_params.txt"
    energy_stem = candidate_root / "fiberdock_energies"
    rewrite_declared_stem(source_params, params_path, energy_stem)
    start = time.perf_counter()
    result = {
        "label": label,
        "source_params": str(source_params),
        "params": str(params_path),
        "declared_energy_stem": str(energy_stem),
        "tool_root": str(tools_root),
    }
    try:
        completed = subprocess.run(
            [str(tools_root / "FiberDock"), str(params_path), "0"],
            cwd=tools_root,
            capture_output=True,
            text=True,
            timeout=900,
            check=False,
        )
        result.update({
            "return_code": completed.returncode,
            "status": "completed" if completed.returncode == 0 else "failed",
            "stdout_tail": completed.stdout[-4000:],
            "stderr_tail": completed.stderr[-4000:],
        })
    except Exception as exc:
        result.update({
            "return_code": None,
            "status": "exception",
            "error": f"{type(exc).__name__}: {exc}",
        })
    result["elapsed_seconds"] = time.perf_counter() - start
    ref_path = Path(str(energy_stem) + ".ref")
    result["declared_ref"] = str(ref_path)
    result["declared_ref_exists"] = ref_path.is_file()
    result["declared_solution"] = parse_energy(ref_path)
    pdb_paths = sorted(candidate_root.glob("fiberdock_energies*.ref.pdb"))
    result["refined_pdbs"] = [pdb_integrity(path) for path in pdb_paths]
    result["valid_refined_pdb_count"] = sum(item["valid"] for item in result["refined_pdbs"])
    result["contract_status"] = (
        "corrected_declared_output_and_valid_pdb"
        if result["declared_ref_exists"] and result["valid_refined_pdb_count"]
        else "missing_declared_output_or_invalid_pdb"
    )
    return result


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_fiberdock_corrected_replay.py OUTPUT_ROOT")
    output_root = Path(sys.argv[1]).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    tasks = [
        (label, REPO_ROOT / source_arm, relative_params, output_root)
        for label, source_arm, relative_params in CANDIDATES
    ]
    workers = min(len(tasks), max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(lambda task: run_candidate(*task), tasks))
    report = {
        "experiment": "fiberdock_corrected_output_contract_replay",
        "workers": workers,
        "records": sorted(records, key=lambda row: row["label"]),
        "summary": {
            "candidate_count": len(records),
            "completed_count": sum(row["status"] == "completed" and row["return_code"] == 0 for row in records),
            "declared_outputs": sum(row["declared_ref_exists"] for row in records),
            "valid_refined_pdbs": sum(row["valid_refined_pdb_count"] for row in records),
            "corrected_energy_values": [
                row["declared_solution"]["global_energy"]
                for row in records
                if row["declared_solution"] is not None
            ],
        },
    }
    output_path = output_root / "fiberdock_corrected_replay.json"
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["summary"]["completed_count"] == len(records) else 2


if __name__ == "__main__":
    raise SystemExit(main())
