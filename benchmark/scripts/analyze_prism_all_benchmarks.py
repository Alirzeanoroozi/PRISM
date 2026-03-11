#!/usr/bin/env python3
"""
Run PRISM-vs-benchmark analysis for T_Rigid, T_medium, and T_difficult.

This is a thin wrapper around:
    benchmark/scripts/analyze_prism_rigid_results.py

It reuses the same analyzer with different CSV inputs and separate output/native folders.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_one(
    analyzer_script,
    rosetta_root,
    csv_path,
    out_dir,
    native_dir,
    metrics,
    score_python,
    download_native_missing,
    verbose_every,
    model_fix_map_csv,
    score_timeout_sec,
):
    cmd = [
        sys.executable,
        str(analyzer_script),
        "--rosetta-root",
        str(rosetta_root),
        "--t-rigid-csv",  # analyzer accepts arbitrary benchmark CSV despite option name
        str(csv_path),
        "--out-dir",
        str(out_dir),
        "--native-dir",
        str(native_dir),
        "--metrics",
        metrics,
        "--score-python",
        score_python,
        "--verbose-every",
        str(verbose_every),
        "--score-timeout-sec",
        str(score_timeout_sec),
    ]
    if download_native_missing:
        cmd.append("--download-native-missing")
    if model_fix_map_csv:
        cmd.extend(["--model-fix-map-csv", str(model_fix_map_csv)])

    print("[RUN]", " ".join(cmd))
    proc = subprocess.run(cmd)
    return proc.returncode


def main():
    parser = argparse.ArgumentParser(description="Run PRISM analysis for rigid/medium/difficult benchmark CSVs")
    parser.add_argument(
        "--rosetta-root",
        default="benchmark/prism_processed/rosetta_output_1",
        help="Root directory with PRISM prediction PDBs",
    )
    parser.add_argument(
        "--score-python",
        default="benchmark/prism_processed_results/prism_score_env/bin/python",
        help="Python interpreter used by the analyzer for DockQ/iRMSD scoring",
    )
    parser.add_argument(
        "--metrics",
        default="dockq,irmsd",
        help="Metrics passed through to the analyzer (e.g., dockq,irmsd)",
    )
    parser.add_argument(
        "--download-native-missing",
        action="store_true",
        help="Download native bound complex PDBs for each benchmark set",
    )
    parser.add_argument(
        "--verbose-every",
        type=int,
        default=50,
        help="Progress frequency passed to the analyzer",
    )
    parser.add_argument(
        "--sets",
        default="rigid,medium,difficult",
        help="Comma-separated benchmark sets to run: rigid,medium,difficult",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue with remaining sets if one set fails",
    )
    parser.add_argument(
        "--model-fix-map-csv",
        default=None,
        help="Optional CSV from fix_model_chain_names.py to use corrected model files/chains during scoring",
    )
    parser.add_argument(
        "--score-timeout-sec",
        type=int,
        default=120,
        help="Timeout (seconds) for each DockQ/iRMSD subprocess call",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    analyzer_script = repo_root / "benchmark/scripts/analyze_prism_rigid_results.py"
    rosetta_root = repo_root / args.rosetta_root
    score_python = str((repo_root / args.score_python) if not Path(args.score_python).is_absolute() else Path(args.score_python))

    set_map = {
        "rigid": {
            "csv": repo_root / "benchmark/data/T_Rigid.csv",
            "out": repo_root / "benchmark/prism_processed_results/prism_rigid_analysis_all_jobs",
            "native": repo_root / "benchmark/prism_processed_results/native_bound_complexes_t_rigid",
        },
        "medium": {
            "csv": repo_root / "benchmark/data/T_medium.csv",
            "out": repo_root / "benchmark/prism_processed_results/prism_medium_analysis_all_jobs",
            "native": repo_root / "benchmark/prism_processed_results/native_bound_complexes_t_medium",
        },
        "difficult": {
            "csv": repo_root / "benchmark/data/T_difficult.csv",
            "out": repo_root / "benchmark/prism_processed_results/prism_difficult_analysis_all_jobs",
            "native": repo_root / "benchmark/prism_processed_results/native_bound_complexes_t_difficult",
        },
    }

    requested = [x.strip().lower() for x in args.sets.split(",") if x.strip()]
    invalid = [x for x in requested if x not in set_map]
    if invalid:
        print("[ERROR] Unknown set(s): %s" % ", ".join(invalid), file=sys.stderr)
        return 2

    failures = []
    for name in requested:
        cfg = set_map[name]
        print("\n[SET] %s" % name)
        rc = run_one(
            analyzer_script=analyzer_script,
            rosetta_root=rosetta_root,
            csv_path=cfg["csv"],
            out_dir=cfg["out"],
            native_dir=cfg["native"],
            metrics=args.metrics,
            score_python=score_python,
            download_native_missing=args.download_native_missing,
            verbose_every=args.verbose_every,
            model_fix_map_csv=args.model_fix_map_csv,
            score_timeout_sec=args.score_timeout_sec,
        )
        if rc != 0:
            failures.append((name, rc))
            print("[FAIL] %s exited with code %s" % (name, rc), file=sys.stderr)
            if not args.continue_on_error:
                return rc

    if failures:
        print("[DONE] Completed with failures: %s" % failures, file=sys.stderr)
        return 1

    print("\n[DONE] Completed sets: %s" % ", ".join(requested))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
