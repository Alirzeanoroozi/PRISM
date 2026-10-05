#!/usr/bin/env python3
"""
Validate PRISM benchmark analysis outputs.

Checks:
1) per_prediction -> pair_summary aggregation consistency
2) PDB ID matching consistency for matched rows
3) Optional score recomputation (DockQ/iRMSD) on a sample of matched rows
"""

import argparse
import csv
import json
import math
import os
import statistics
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple


SETS = ("rigid", "medium", "difficult")


def to_float(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if s in ("", "None", "NA"):
        return None
    try:
        return float(s)
    except Exception:
        return None


def parse_pair_key(k: str) -> Tuple[str, str]:
    parts = (k or "").split("|")
    if len(parts) != 2:
        return "", ""
    return parts[0].lower(), parts[1].lower()


def unordered_pair(a: str, b: str) -> Tuple[str, str]:
    a = (a or "").lower()
    b = (b or "").lower()
    return (a, b) if a <= b else (b, a)


def mean_var(vals: List[float]) -> Tuple[Optional[float], Optional[float]]:
    if not vals:
        return None, None
    if len(vals) == 1:
        return float(vals[0]), 0.0
    return float(sum(vals)) / float(len(vals)), float(statistics.pvariance(vals))


def run_cmd(cmd: List[str], timeout_sec: int, env: Optional[dict] = None) -> Tuple[int, str, str]:
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=timeout_sec,
        env=env,
    )
    return proc.returncode, proc.stdout or "", proc.stderr or ""


def run_dockq(score_python: str, dockq_script: Path, model_pdb: str, native_pdb: str, mapping: str, timeout_sec: int) -> Tuple[Optional[float], str]:
    py_bin = str(Path(score_python).resolve().parent)
    env = dict(os.environ)
    env["PATH"] = py_bin + os.pathsep + env.get("PATH", "")
    rc, out, err = run_cmd(
        [score_python, str(dockq_script), model_pdb, native_pdb, "--mapping", mapping],
        timeout_sec=timeout_sec,
        env=env,
    )
    if rc != 0:
        msg = (err or out).strip().replace("\n", " | ")
        return None, f"exit_{rc}:{msg}"
    val = None
    for line in out.splitlines():
        line = line.strip()
        if line.startswith("DockQ:"):
            try:
                val = float(line.split(":", 1)[1].strip())
            except Exception:
                pass
    if val is None:
        return None, "parse_error"
    return val, ""


def run_irmsd(score_python: str, irmsd_script: Path, model_pdb: str, model_r: str, model_l: str, native_pdb: str, ref_r: str, ref_l: str, timeout_sec: int) -> Tuple[Optional[float], str]:
    py_bin = str(Path(score_python).resolve().parent)
    env = dict(os.environ)
    env["PATH"] = py_bin + os.pathsep + env.get("PATH", "")
    rc, out, err = run_cmd(
        [score_python, str(irmsd_script), model_pdb, model_r, model_l, native_pdb, ref_r, ref_l],
        timeout_sec=timeout_sec,
        env=env,
    )
    if rc != 0:
        msg = (err or out).strip().replace("\n", " | ")
        return None, f"exit_{rc}:{msg}"
    lines = [ln.strip() for ln in out.splitlines() if ln.strip()]
    if not lines:
        return None, "no_output"
    try:
        return float(lines[-1]), ""
    except Exception:
        return None, f"parse_error:{lines[-1]}"


def validate_set(
    set_name: str,
    results_root: Path,
    score_python: str,
    scripts_dir: Path,
    sample_size: int,
    tolerance: float,
    timeout_sec: int,
) -> Dict:
    per_path = results_root / f"prism_{set_name}_analysis_all_jobs/per_prediction.csv"
    pair_path = results_root / f"prism_{set_name}_analysis_all_jobs/pair_summary.csv"

    with per_path.open(newline="") as fh:
        per_rows = list(csv.DictReader(fh))
    with pair_path.open(newline="") as fh:
        pair_rows = list(csv.DictReader(fh))

    errors = []
    warnings = []

    # 1) PDB matching consistency
    matched_rows = [r for r in per_rows if (r.get("pair_match_status") or "").strip() != "no_t_rigid_match"]
    for r in matched_rows:
        tpl_pair = unordered_pair(r.get("tpl1_pdb4", ""), r.get("tpl2_pdb4", ""))
        pair_key = parse_pair_key(r.get("pdb_pair_key", ""))
        if tpl_pair != pair_key:
            errors.append(f"pair_key_mismatch:{r.get('filename')} tpl_pair={tpl_pair} pair_key={pair_key}")
        p1 = (r.get("pdb1_raw") or "")[:4].lower()
        p2 = (r.get("pdb2_raw") or "")[:4].lower()
        csv_pair = unordered_pair(p1, p2)
        if csv_pair != pair_key:
            errors.append(f"csv_pair_mismatch:{r.get('filename')} csv_pair={csv_pair} pair_key={pair_key}")

    # 2) Aggregation consistency
    grouped: Dict[Tuple[str, str, str, str, str], List[dict]] = {}
    for r in per_rows:
        if not (r.get("complex_raw") or "").strip():
            continue
        key = (
            r.get("pdb_pair_key") or "",
            r.get("pdb1_raw") or "",
            r.get("pdb2_raw") or "",
            r.get("complex_raw") or "",
            r.get("bound_pdb4") or "",
        )
        grouped.setdefault(key, []).append(r)

    pair_index = {}
    for r in pair_rows:
        key = (
            r.get("pdb_pair_key") or "",
            r.get("pdb1_raw") or "",
            r.get("pdb2_raw") or "",
            r.get("complex_raw") or "",
            r.get("bound_pdb4") or "",
        )
        pair_index[key] = r

    for key, rows in grouped.items():
        if key not in pair_index:
            errors.append(f"missing_pair_summary_row:{key}")
            continue
        ps = pair_index[key]
        dockq_vals = [to_float(r.get("dockq")) for r in rows]
        dockq_vals = [v for v in dockq_vals if v is not None]
        ir_vals = [to_float(r.get("irmsd_min")) for r in rows]
        ir_vals = [v for v in ir_vals if v is not None]
        bb_vals = [to_float(r.get("irmsd_backbone_min")) for r in rows]
        bb_vals = [v for v in bb_vals if v is not None]

        checks = [
            (int(ps.get("n_predictions_total") or 0), len(rows), "n_predictions_total"),
            (int(ps.get("n_predictions_mapped") or 0), len(rows), "n_predictions_mapped"),
            (int(ps.get("n_dockq_scored") or 0), len(dockq_vals), "n_dockq_scored"),
            (int(ps.get("n_irmsd_scored") or 0), len(ir_vals), "n_irmsd_scored"),
            (int(ps.get("n_irmsd_backbone_scored") or 0), len(bb_vals), "n_irmsd_backbone_scored"),
        ]
        for got, exp, name in checks:
            if got != exp:
                errors.append(f"{name}_mismatch:{key} got={got} expected={exp}")

        for vals, pref, best_mode in (
            (dockq_vals, "dockq", "max"),
            (ir_vals, "irmsd", "min"),
            (bb_vals, "irmsd_backbone", "min"),
        ):
            m, v = mean_var(vals)
            if m is not None:
                pm = to_float(ps.get(f"{pref}_mean"))
                pv = to_float(ps.get(f"{pref}_variance"))
                pb = to_float(ps.get(f"{pref}_best" if pref == "dockq" else f"{pref}_best_min"))
                eb = max(vals) if best_mode == "max" else min(vals)
                if pm is None or abs(pm - m) > tolerance:
                    errors.append(f"{pref}_mean_mismatch:{key} got={pm} expected={m}")
                if pv is None or abs(pv - v) > tolerance:
                    errors.append(f"{pref}_var_mismatch:{key} got={pv} expected={v}")
                if pb is None or abs(pb - eb) > tolerance:
                    errors.append(f"{pref}_best_mismatch:{key} got={pb} expected={eb}")

    # 3) Recompute sample
    py_bin = str(Path(score_python).resolve().parent)
    check_env = dict(os.environ)
    check_env["PATH"] = py_bin + os.pathsep + check_env.get("PATH", "")
    dockq_available = subprocess.run(
        ["which", "DockQ"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        env=check_env,
    ).returncode == 0
    if not dockq_available:
        warnings.append("dockq_binary_not_found_in_score_env_path:DockQ recomputation checks skipped")

    sample_candidates = [
        r for r in matched_rows
        if to_float(r.get("dockq")) is not None
        and to_float(r.get("irmsd_min")) is not None
        and not (r.get("dockq_error") or "").strip()
        and not (r.get("irmsd_error") or "").strip()
    ]
    sample_rows = sorted(sample_candidates, key=lambda r: (r.get("filename") or "", r.get("complex_raw") or ""))[:sample_size]

    dockq_script = scripts_dir / "dockq.py"
    irmsd_script = scripts_dir / "irmsd.py"
    sample_checked = 0
    sample_failed = 0

    for r in sample_rows:
        model_pdb = r.get("model_path_used") or r.get("model_path")
        native_pdb = r.get("native_pdb")
        model_r = (r.get("model_receptor_chains") or "").strip()
        model_l = (r.get("model_ligand_chains") or "").strip()
        if not model_pdb or not native_pdb or not model_r or not model_l:
            warnings.append(f"sample_skipped_missing_paths_or_chains:{r.get('filename')}")
            continue
        pairs = [x for x in (r.get("native_pairs_evaluated") or "").split(",") if len(x) >= 2]
        if not pairs:
            warnings.append(f"sample_skipped_no_native_pairs:{r.get('filename')}")
            continue

        dvals = []
        ivals = []
        sample_checked += 1
        row_fail = False
        for p in pairs:
            ref_r, ref_l = p[0], p[1]
            if dockq_available:
                dval, derr = run_dockq(score_python, dockq_script, model_pdb, native_pdb, f"{model_r}{model_l}:{ref_r}{ref_l}", timeout_sec)
                if dval is not None:
                    dvals.append(dval)
                elif derr:
                    row_fail = True
            fwd, ferr = run_irmsd(score_python, irmsd_script, model_pdb, model_r, model_l, native_pdb, ref_r, ref_l, timeout_sec)
            if fwd is not None:
                ivals.append(fwd)
            elif ferr:
                row_fail = True
            rev, rerr = run_irmsd(score_python, irmsd_script, model_pdb, model_l, model_r, native_pdb, ref_l, ref_r, timeout_sec)
            if rev is not None:
                ivals.append(rev)
            elif rerr:
                row_fail = True

        expected_d = to_float(r.get("dockq"))
        expected_i = to_float(r.get("irmsd_min"))
        got_d = max(dvals) if dvals else None
        got_i = min(ivals) if ivals else None

        if dockq_available:
            if expected_d is None or got_d is None or abs(expected_d - got_d) > tolerance:
                sample_failed += 1
                errors.append(
                    f"sample_dockq_mismatch:{r.get('filename')} expected={expected_d} got={got_d} complex={r.get('complex_raw')}"
                )
        if expected_i is None or got_i is None or abs(expected_i - got_i) > tolerance:
            sample_failed += 1
            errors.append(
                f"sample_irmsd_mismatch:{r.get('filename')} expected={expected_i} got={got_i} complex={r.get('complex_raw')}"
            )
        if row_fail:
            warnings.append(f"sample_partial_runtime_failures:{r.get('filename')}")

    return {
        "set": set_name,
        "per_rows": len(per_rows),
        "pair_rows": len(pair_rows),
        "matched_rows": len(matched_rows),
        "sample_requested": sample_size,
        "sample_checked": sample_checked,
        "sample_failed_checks": sample_failed,
        "n_errors": len(errors),
        "n_warnings": len(warnings),
        "errors": errors,
        "warnings": warnings,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate PRISM pipeline outputs")
    parser.add_argument("--results-root", default="benchmark/prism_processed_results")
    parser.add_argument("--score-python", default="benchmark/prism_processed_results/prism_score_env/bin/python")
    parser.add_argument("--sample-size-per-set", type=int, default=10)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--timeout-sec", type=int, default=30)
    parser.add_argument("--out-dir", default="benchmark/prism_processed_results/validation")
    args = parser.parse_args()

    results_root = Path(args.results_root)
    scripts_dir = Path(__file__).resolve().parent
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for s in SETS:
        all_results.append(
            validate_set(
                set_name=s,
                results_root=results_root,
                score_python=args.score_python,
                scripts_dir=scripts_dir,
                sample_size=args.sample_size_per_set,
                tolerance=args.tolerance,
                timeout_sec=args.timeout_sec,
            )
        )

    total_errors = sum(r["n_errors"] for r in all_results)
    total_warnings = sum(r["n_warnings"] for r in all_results)

    json_path = out_dir / "pipeline_validation.json"
    md_path = out_dir / "pipeline_validation.md"
    with json_path.open("w") as fh:
        json.dump(
            {
                "total_errors": total_errors,
                "total_warnings": total_warnings,
                "results": all_results,
            },
            fh,
            indent=2,
        )

    lines = []
    lines.append("# PRISM Pipeline Validation")
    lines.append("")
    lines.append(f"- Total errors: `{total_errors}`")
    lines.append(f"- Total warnings: `{total_warnings}`")
    lines.append("")
    lines.append("| Set | per_rows | pair_rows | matched_rows | sample_checked | sample_failed_checks | errors | warnings |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in all_results:
        lines.append(
            f"| {r['set']} | {r['per_rows']} | {r['pair_rows']} | {r['matched_rows']} | "
            f"{r['sample_checked']} | {r['sample_failed_checks']} | {r['n_errors']} | {r['n_warnings']} |"
        )
    lines.append("")
    for r in all_results:
        lines.append(f"## {r['set'].title()}")
        lines.append("")
        lines.append("Top errors:")
        if r["errors"]:
            for e in r["errors"][:20]:
                lines.append(f"- `{e}`")
        else:
            lines.append("- None")
        lines.append("")
        lines.append("Top warnings:")
        if r["warnings"]:
            for w in r["warnings"][:20]:
                lines.append(f"- `{w}`")
        else:
            lines.append("- None")
        lines.append("")

    md_path.write_text("\n".join(lines) + "\n")

    print(f"[INFO] Wrote {json_path}")
    print(f"[INFO] Wrote {md_path}")
    print(f"[INFO] total_errors={total_errors} total_warnings={total_warnings}")
    return 0 if total_errors == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
