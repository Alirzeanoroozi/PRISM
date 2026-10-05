#!/usr/bin/env python3
"""
Non-destructive preflight check for the benchmark pipeline.

This script avoids production result paths and only writes temporary files under
the chosen temp root. It verifies:
1. Python syntax for core scripts
2. Rosetta-output benchmark entry points start correctly
3. Batch script shell syntax
4. Required files and external commands exist
5. Single-file and folder scoring smoke tests produce CSV outputs
"""

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run(cmd, cwd, timeout=300):
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=timeout,
    )
    return proc.returncode, proc.stdout.strip(), proc.stderr.strip()


def check(condition, label, detail=""):
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}")
    if detail:
        print(f"       {detail}")
    return condition


def skip(label, detail=""):
    print(f"[SKIP] {label}")
    if detail:
        print(f"       {detail}")


def read_first_csv_row(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    return rows[0] if rows else None


def main():
    parser = argparse.ArgumentParser(description="Run a safe, non-destructive preflight test of the benchmark pipeline")
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root",
    )
    parser.add_argument(
        "--score-python",
        default="benchmark/prism_processed/env/prism_score_env/bin/python",
        help="Python interpreter for scoring tools",
    )
    parser.add_argument(
        "--example-model",
        default="2dvmCD_1FGNH_1TFHA_R.pdb",
        help="Example model PDB for smoke test",
    )
    parser.add_argument(
        "--example-native",
        default="1AHW_r_b.pdb",
        help="Example native PDB for smoke test",
    )
    parser.add_argument(
        "--temp-root",
        default=None,
        help="Optional temporary output root. Defaults to a new temp directory.",
    )
    parser.add_argument(
        "--timeout-sec",
        type=int,
        default=300,
        help="Timeout for each subprocess call",
    )
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    score_python = Path(args.score_python)
    if not score_python.is_absolute():
        score_python = repo_root / score_python
    example_model = Path(args.example_model)
    if not example_model.is_absolute():
        example_model = repo_root / example_model
    example_native = Path(args.example_native)
    if not example_native.is_absolute():
        example_native = repo_root / example_native

    temp_root = Path(args.temp_root) if args.temp_root else Path(tempfile.mkdtemp(prefix="prism_benchmark_preflight_"))
    temp_root.mkdir(parents=True, exist_ok=True)

    print(f"[info] repo_root={repo_root}")
    print(f"[info] score_python={score_python}")
    print(f"[info] temp_root={temp_root}")

    ok = True

    core_scripts = [
        repo_root / "benchmark/scripts/dockq.py",
        repo_root / "benchmark/scripts/irmsd.py",
        repo_root / "benchmark/scripts/score_single_prism_pair.py",
        repo_root / "benchmark/scripts/validate_prism_pipeline.py",
        repo_root / "benchmark/scripts/rosetta_output/analyze_prism_all_benchmarks.py",
        repo_root / "benchmark/scripts/rosetta_output/analyze_prism_rigid_results.py",
        repo_root / "benchmark/scripts/rosetta_output/fix_model_chain_names.py",
    ]
    rc, so, se = run(
        ["python3", "-m", "py_compile", *[str(p) for p in core_scripts]],
        cwd=repo_root,
        timeout=args.timeout_sec,
    )
    ok &= check(rc == 0, "core script syntax", se or so)

    entrypoint_checks = [
        ["python3", "benchmark/scripts/rosetta_output/fix_model_chain_names.py", "--help"],
        ["python3", "benchmark/scripts/rosetta_output/analyze_prism_all_benchmarks.py", "--help"],
        ["python3", "benchmark/scripts/rosetta_output/analyze_prism_rigid_results.py", "--help"],
    ]
    for cmd in entrypoint_checks:
        rc, so, se = run(cmd, cwd=repo_root, timeout=args.timeout_sec)
        ok &= check(rc == 0, "startup " + " ".join(cmd[1:3]), se or so)

    rc, so, se = run(
        ["bash", "-n", "benchmark/scripts/rosetta_output/submit_prism_analysis_all.sbatch"],
        cwd=repo_root,
        timeout=args.timeout_sec,
    )
    ok &= check(rc == 0, "batch script syntax", se or so)

    required_paths = [
        repo_root / "benchmark/scripts/rosetta_output/submit_prism_analysis_all.sbatch",
        repo_root / "benchmark/data/T_Rigid.csv",
        repo_root / "benchmark/data/T_medium.csv",
        repo_root / "benchmark/data/T_difficult.csv",
    ]
    missing = [str(p) for p in required_paths if not p.exists()]
    ok &= check(not missing, "required files present", ", ".join(missing))

    required_commands = ["python3", "bash", "unzip", "tar", "find", "head"]
    available = []
    missing_cmds = []
    for name in required_commands:
        if shutil.which(name):
            available.append(name)
        else:
            missing_cmds.append(name)
    ok &= check(not missing_cmds, "required external commands", f"available={','.join(available)} missing={','.join(missing_cmds)}")

    smoke_missing = [
        str(p) for p in [score_python, example_model, example_native] if not p.exists()
    ]
    if smoke_missing:
        skip("scoring smoke tests", "missing optional smoke-test inputs: " + ", ".join(smoke_missing))
    else:
        single_csv = temp_root / "single_scores.csv"
        folder_dir = temp_root / "models"
        folder_dir.mkdir(exist_ok=True)
        shutil.copy2(example_model, folder_dir / example_model.name)
        folder_csv = temp_root / "folder_scores.csv"

        single_cmd = [
            str(score_python),
            "benchmark/scripts/score_single_prism_pair.py",
            str(example_model),
            str(example_native),
            "--dockq-no-align",
            "--out-csv",
            str(single_csv),
        ]
        rc, so, se = run(single_cmd, cwd=repo_root, timeout=args.timeout_sec)
        ok &= check(rc == 0 and single_csv.exists(), "single-file scoring smoke test", se or so)

        folder_cmd = [
            str(score_python),
            "benchmark/scripts/score_single_prism_pair.py",
            str(folder_dir),
            str(example_native),
            "--dockq-no-align",
            "--out-csv",
            str(folder_csv),
        ]
        rc, so, se = run(folder_cmd, cwd=repo_root, timeout=args.timeout_sec)
        ok &= check(rc == 0 and folder_csv.exists(), "folder scoring smoke test", se or so)

        if single_csv.exists():
            row = read_first_csv_row(single_csv)
            detail = ""
            if row:
                detail = f"irmsd={row.get('irmsd')} dockq={row.get('dockq')} error_irmsd={row.get('error_irmsd')} error_dockq={row.get('error_dockq')}"
            ok &= check(bool(row is not None and row.get("irmsd") and row.get("dockq")), "single-file scoring values", detail)

        if folder_csv.exists():
            row = read_first_csv_row(folder_csv)
            detail = ""
            if row:
                detail = f"irmsd={row.get('irmsd')} dockq={row.get('dockq')} error_irmsd={row.get('error_irmsd')} error_dockq={row.get('error_dockq')}"
            ok &= check(bool(row is not None and row.get("irmsd") and row.get("dockq")), "folder scoring values", detail)

    print(f"[done] overall={'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
