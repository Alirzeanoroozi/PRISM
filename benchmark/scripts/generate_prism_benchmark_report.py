#!/usr/bin/env python3
"""
Generate a consolidated cross-benchmark report (rigid/medium/difficult) with plots.

Inputs expected:
  benchmark/prism_processed_results/prism_<set>_analysis_all_jobs/per_prediction.csv
  benchmark/prism_processed_results/prism_<set>_analysis_all_jobs/pair_summary.csv

Outputs:
  benchmark/prism_processed_results/benchmark_reports/cross_benchmark_report.md
  benchmark/prism_processed_results/benchmark_reports/cross_benchmark_report.pdf
  benchmark/prism_processed_results/benchmark_reports/plots/*.png
"""

import argparse
import csv
import math
import os
from pathlib import Path
from statistics import mean


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def fnum(x):
    if x in (None, ""):
        return None
    try:
        return float(x)
    except Exception:
        return None


def count_complexes_in_csv(csv_path):
    if not csv_path.exists():
        return None
    vals = set()
    with csv_path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            c = (row.get("Complex") or row.get("complex") or "").strip()
            if c:
                vals.add(c)
    return len(vals)


def summarize_set(name, base_dir, csv_path=None):
    per_path = base_dir / ("prism_%s_analysis_all_jobs/per_prediction.csv" % name)
    pair_path = base_dir / ("prism_%s_analysis_all_jobs/pair_summary.csv" % name)
    per_rows = read_csv(per_path)
    pair_rows = read_csv(pair_path)

    matched = [r for r in per_rows if (r.get("complex_raw") or "").strip()]
    model_key = lambda r: (r.get("model_path") or r.get("model_path_used") or r.get("filename") or "")
    all_model_keys = {model_key(r) for r in per_rows if model_key(r)}
    matched_model_keys = {model_key(r) for r in matched if model_key(r)}
    total_models = len(all_model_keys)
    matched_models = len(matched_model_keys)
    unmatched_models = total_models - matched_models
    dockq_vals = [fnum(r.get("dockq")) for r in matched]
    dockq_vals = [v for v in dockq_vals if v is not None]
    irmsd_vals = [fnum(r.get("irmsd_min")) for r in matched]
    irmsd_vals = [v for v in irmsd_vals if v is not None]

    capri_counts = {}
    for r in matched:
        c = (r.get("dockq_capri") or "").strip()
        if c:
            capri_counts[c] = capri_counts.get(c, 0) + 1

    def categorize_err(msg, metric):
        msg = (msg or "").strip()
        if not msg:
            return None
        low = msg.lower()
        if metric == "dockq":
            if "could not find interfaces in the native model" in low:
                return "No interface found in native for selected mapping"
            if "dockq not found" in low:
                return "DockQ executable not found"
            if "dockq failed (exit 1)" in low:
                return "DockQ runtime failure (exit 1)"
            if "dockq_parse_error" in low:
                return "DockQ output parse error"
            return "Other DockQ error"
        if metric == "irmsd":
            if "module named 'numpy'" in low:
                return "Missing numpy dependency"
            if "module named 'bio'" in low:
                return "Missing biopython dependency"
            if "missing_model_chain_assignment" in low:
                return "Missing model chain assignment"
            if "native_bound_complex_not_found" in low:
                return "Missing native bound complex PDB"
            if "irmsd_failed_exit_" in low:
                return "iRMSD script runtime failure"
            if "irmsd_parse_error" in low:
                return "iRMSD output parse error"
            return "Other iRMSD error"
        return "Other error"

    dockq_error_counts = {}
    irmsd_error_counts = {}
    for r in matched:
        dk = categorize_err(r.get("dockq_error"), "dockq")
        ir = categorize_err(r.get("irmsd_error"), "irmsd")
        if dk:
            dockq_error_counts[dk] = dockq_error_counts.get(dk, 0) + 1
        if ir:
            irmsd_error_counts[ir] = irmsd_error_counts.get(ir, 0) + 1

    top_dockq_pairs = []
    top_irmsd_pairs = []
    for r in pair_rows:
        dv = fnum(r.get("dockq_best"))
        iv = fnum(r.get("irmsd_best_min"))
        if dv is not None:
            top_dockq_pairs.append((dv, r))
        if iv is not None:
            top_irmsd_pairs.append((iv, r))
    top_dockq_pairs = sorted(top_dockq_pairs, key=lambda x: x[0], reverse=True)[:5]
    top_irmsd_pairs = sorted(top_irmsd_pairs, key=lambda x: x[0])[:5]

    scored_both = []
    for r in matched:
        dv = fnum(r.get("dockq"))
        iv = fnum(r.get("irmsd_min"))
        if dv is not None and iv is not None:
            scored_both.append((dv, iv, r))

    complexes_with_models = sum(1 for r in pair_rows if fnum(r.get("n_predictions_total")) and fnum(r.get("n_predictions_total")) > 0)
    complexes_with_dockq = sum(1 for r in pair_rows if fnum(r.get("n_dockq_scored")) and fnum(r.get("n_dockq_scored")) > 0)
    complexes_with_irmsd = sum(1 for r in pair_rows if fnum(r.get("n_irmsd_scored")) and fnum(r.get("n_irmsd_scored")) > 0)

    matched_without_dockq = [r for r in matched if fnum(r.get("dockq")) is None]
    matched_without_irmsd = [r for r in matched if fnum(r.get("irmsd_min")) is None]
    dockq_error_files = []
    irmsd_error_files = []
    for r in matched:
        de = (r.get("dockq_error") or "").strip()
        ie = (r.get("irmsd_error") or "").strip()
        if de:
            dockq_error_files.append({
                "filename": r.get("filename") or "",
                "complex_raw": r.get("complex_raw") or "",
                "target_token": r.get("target_token") or r.get("pdb_pair_key") or "",
                "error": de,
            })
        if ie:
            irmsd_error_files.append({
                "filename": r.get("filename") or "",
                "complex_raw": r.get("complex_raw") or "",
                "target_token": r.get("target_token") or r.get("pdb_pair_key") or "",
                "error": ie,
            })

    return {
        "name": name,
        "complexes_in_csv": count_complexes_in_csv(csv_path) if csv_path else None,
        "per_rows": per_rows,
        "pair_rows": pair_rows,
        "matched_rows": matched,
        "total_predictions": len(per_rows),
        "matched_predictions": len(matched),
        "unmatched_predictions": sum(1 for r in per_rows if (r.get("pair_match_status") or "") == "no_t_rigid_match"),
        "total_models": total_models,
        "matched_models": matched_models,
        "unmatched_models": unmatched_models,
        "dockq_scored": len(dockq_vals),
        "irmsd_scored": len(irmsd_vals),
        "dockq_errors": sum(1 for r in matched if (r.get("dockq_error") or "").strip()),
        "irmsd_errors": sum(1 for r in matched if (r.get("irmsd_error") or "").strip()),
        "dockq_values": dockq_vals,
        "irmsd_values": irmsd_vals,
        "scored_both": scored_both,
        "capri_counts": capri_counts,
        "dockq_error_counts": dockq_error_counts,
        "irmsd_error_counts": irmsd_error_counts,
        "pair_rows_count": len(pair_rows),
        "complexes_with_models": complexes_with_models,
        "complexes_with_dockq": complexes_with_dockq,
        "complexes_with_irmsd": complexes_with_irmsd,
        "matched_without_dockq": len(matched_without_dockq),
        "matched_without_irmsd": len(matched_without_irmsd),
        "dockq_error_files": dockq_error_files,
        "irmsd_error_files": irmsd_error_files,
        "top_dockq_pairs": top_dockq_pairs,
        "top_irmsd_pairs": top_irmsd_pairs,
    }


def make_plots(data_by_set, plots_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plots_dir.mkdir(parents=True, exist_ok=True)

    sets = ["rigid", "medium", "difficult"]
    colors = {"rigid": "#4C78A8", "medium": "#F28E2B", "difficult": "#76B7B2"}

    for s in sets:
        d = data_by_set[s]
        c = colors[s]

        # 1) per-set status bar (found/scored/errors)
        labels = ["Matched", "DockQ Scored", "iRMSD Scored", "DockQ Errors", "iRMSD Errors"]
        vals = [d["matched_predictions"], d["dockq_scored"], d["irmsd_scored"], d["dockq_errors"], d["irmsd_errors"]]
        plt.figure(figsize=(8, 4.8))
        plt.bar(labels, vals, color=[c, "#59A14F", "#E15759", "#B07AA1", "#9C755F"])
        plt.ylabel("Count")
        plt.title("%s Set: Model/Score Status" % s.title())
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(str(plots_dir / ("%s_status_bar.png" % s)), dpi=200)
        plt.close()

        # 2) per-set DockQ histogram
        dockq_vals = d["dockq_values"]
        if dockq_vals:
            bins = [i / 20.0 for i in range(21)]
            plt.figure(figsize=(8, 4.8))
            plt.hist(dockq_vals, bins=bins, color=c, edgecolor="white")
            plt.xlabel("DockQ")
            plt.ylabel("Count")
            plt.title("%s Set: DockQ Distribution" % s.title())
            plt.tight_layout()
            plt.savefig(str(plots_dir / ("%s_dockq_histogram.png" % s)), dpi=200)
            plt.close()

        # 3) per-set iRMSD histogram (capped)
        irmsd_vals = [min(v, 30.0) for v in d["irmsd_values"]]
        if irmsd_vals:
            bins = [i for i in range(0, 31, 2)]
            plt.figure(figsize=(8, 4.8))
            plt.hist(irmsd_vals, bins=bins, color=c, edgecolor="white")
            plt.xlabel("iRMSD (A), values >30 clipped to 30")
            plt.ylabel("Count")
            plt.title("%s Set: iRMSD Distribution" % s.title())
            plt.tight_layout()
            plt.savefig(str(plots_dir / ("%s_irmsd_histogram.png" % s)), dpi=200)
            plt.close()

        # 4) per-set DockQ vs iRMSD scatter
        pts = d["scored_both"]
        if pts:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            plt.figure(figsize=(8, 5.2))
            plt.scatter(xs, ys, alpha=0.85, s=34, color=c, edgecolors="none")
            plt.xlabel("DockQ")
            plt.ylabel("iRMSD (A)")
            plt.title("%s Set: DockQ vs iRMSD" % s.title())
            plt.xlim(-0.02, 1.02)
            plt.ylim(bottom=0)
            plt.tight_layout()
            plt.savefig(str(plots_dir / ("%s_dockq_vs_irmsd.png" % s)), dpi=200)
            plt.close()


def write_error_csv(data_by_set, out_dir):
    out_path = out_dir / "error_files.csv"
    rows = []
    for s in ["rigid", "medium", "difficult"]:
        for e in data_by_set[s]["dockq_error_files"]:
            rows.append({
                "set": s,
                "metric": "dockq",
                "filename": e["filename"],
                "complex_raw": e["complex_raw"],
                "target_token": e["target_token"],
                "error_reason": e["error"],
            })
        for e in data_by_set[s]["irmsd_error_files"]:
            rows.append({
                "set": s,
                "metric": "irmsd",
                "filename": e["filename"],
                "complex_raw": e["complex_raw"],
                "target_token": e["target_token"],
                "error_reason": e["error"],
            })
    fieldnames = ["set", "metric", "filename", "complex_raw", "target_token", "error_reason"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return out_path


def fmt(x, digits=3):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return ("%%.%df" % digits) % x
    return str(x)


def infer_insights(data_by_set):
    insights = []
    sets = ["rigid", "medium", "difficult"]
    for s in sets:
        d = data_by_set[s]
        match_rate = (100.0 * d["matched_models"] / d["total_models"]) if d["total_models"] else 0.0
        dockq_rate = (100.0 * d["dockq_scored"] / d["matched_predictions"]) if d["matched_predictions"] else 0.0
        irmsd_rate = (100.0 * d["irmsd_scored"] / d["matched_predictions"]) if d["matched_predictions"] else 0.0
        insights.append(
            "%s: %s/%s unique models matched the benchmark set (%.1f%%); DockQ scored %.1f%% of matched rows and iRMSD scored %.1f%%."
            % (s.title(), d["matched_models"], d["total_models"], match_rate, dockq_rate, irmsd_rate)
        )

    # Best overall DockQ and iRMSD across sets
    all_dockq = []
    all_irmsd = []
    for s in sets:
        for r in data_by_set[s]["matched_rows"]:
            dv = fnum(r.get("dockq"))
            iv = fnum(r.get("irmsd_min"))
            if dv is not None:
                all_dockq.append((dv, s, r))
            if iv is not None:
                all_irmsd.append((iv, s, r))
    if all_dockq:
        best_d = max(all_dockq, key=lambda x: x[0])
        r = best_d[2]
        insights.append(
            "Best DockQ across all analyzed rows was %.3f in %s (%s; pair %s, model chains %s/%s, native pair %s:%s)."
            % (
                best_d[0], best_d[1].title(), r["filename"], row_pair_id(r),
                r.get("model_receptor_chains") or "?", r.get("model_ligand_chains") or "?",
                r.get("ref_receptor_chain") or "?", r.get("ref_ligand_chain") or "?"
            )
        )
    if all_irmsd:
        best_i = min(all_irmsd, key=lambda x: x[0])
        r = best_i[2]
        insights.append(
            "Best (lowest) iRMSD across all analyzed rows was %.3f A in %s (%s; pair %s)."
            % (best_i[0], best_i[1].title(), r["filename"], row_pair_id(r))
        )

    return insights


def find_case_row(rows, filename):
    for r in rows:
        if (r.get("filename") or "") == filename:
            return r
    return None


def row_pair_id(r):
    return (r.get("pdb_pair_key") or r.get("target_token") or "").strip() or "NA"


def write_markdown(data_by_set, report_dir, plots_dir):
    sets = ["rigid", "medium", "difficult"]
    md_path = report_dir / "cross_benchmark_report.md"
    lines = []
    lines.append("# PRISM Benchmark Analysis Report (Rigid / Medium / Difficult)")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append("This report summarizes PRISM prediction analysis against `T_Rigid.csv`, `T_medium.csv`, and `T_difficult.csv` using `DockQ` and `iRMSD`.")
    lines.append("")
    lines.append("## Method Logic (Why Scores Are Not Read From PRISM Output Files)")
    lines.append("")
    lines.append("- The target metrics in this analysis are **benchmark comparison metrics** (`DockQ`, `iRMSD`) that require a comparison between:")
    lines.append("- a PRISM predicted complex model (`.pdb`) and")
    lines.append("- the benchmark native bound complex PDB with a specific chain mapping from the benchmark CSV.")
    lines.append("- PRISM output files (e.g. Rosetta score files / intermediate files) do not directly provide these benchmark-comparison metrics with the required benchmark-native chain mapping.")
    lines.append("- Therefore, the pipeline recomputes `DockQ` and `iRMSD` from the model PDB + downloaded native bound PDB + CSV-defined chain roles.")
    lines.append("- The chain fields reported in `per_prediction.csv` are taken from the model filename (or inferred fallback if omitted), as requested.")
    lines.append("")
    lines.append("## Output Locations")
    lines.append("")
    for s in sets:
        lines.append("- `%s`" % ("benchmark/prism_processed_results/prism_%s_analysis_all_jobs" % s))
    lines.append("- Plots: `%s`" % str(plots_dir))
    lines.append("")
    lines.append("## Summary Table")
    lines.append("")
    lines.append("| Set | Complexes in CSV | Total Models | Matched Models | Unmatched Models | DockQ Scored Rows | iRMSD Scored Rows | DockQ Errors | iRMSD Errors |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for s in sets:
        d = data_by_set[s]
        lines.append(
            "| %s | %s | %d | %d | %d | %d | %d | %d | %d |" % (
                s.title(), str(d["complexes_in_csv"]) if d["complexes_in_csv"] is not None else "NA",
                d["total_models"], d["matched_models"], d["unmatched_models"],
                d["dockq_scored"], d["irmsd_scored"], d["dockq_errors"], d["irmsd_errors"]
            )
        )
    lines.append("")
    lines.append("Abbreviation guide:")
    lines.append("- `DockQ` = DockQ score.")
    lines.append("- `iRMSD` = interface RMSD.")
    lines.append("- `Scored Rows` = matched model-complex rows where numeric score was produced.")
    lines.append("")
    lines.append("## Predicted Model Membership by Set")
    lines.append("")
    total_models = data_by_set["rigid"]["total_models"]
    lines.append("| Set | Models Belonging to Set | Percent of All Predicted Models |")
    lines.append("|---|---:|---:|")
    for s in sets:
        d = data_by_set[s]
        pct = (100.0 * d["matched_models"] / total_models) if total_models else 0.0
        lines.append("| %s | %d | %.2f%% |" % (s.title(), d["matched_models"], pct))
    lines.append("")
    lines.append("Abbreviation guide:")
    lines.append("- `Set` = benchmark subset (`Rigid`, `Medium`, `Difficult`).")
    lines.append("- `Models Belonging to Set` = unique model files whose PDB ID pair matches at least one row in that set.")
    lines.append("")
    lines.append("## Complex-Level Coverage")
    lines.append("")
    lines.append("| Set | Complexes With Found Models | Complexes With >=1 DockQ | Complexes With >=1 iRMSD | Matched Models Missing DockQ | Matched Models Missing iRMSD |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for s in sets:
        d = data_by_set[s]
        lines.append(
            "| %s | %d | %d | %d | %d | %d |" % (
                s.title(),
                d["complexes_with_models"],
                d["complexes_with_dockq"],
                d["complexes_with_irmsd"],
                d["matched_without_dockq"],
                d["matched_without_irmsd"],
            )
        )
    lines.append("")
    lines.append("Abbreviation guide:")
    lines.append("- `Cx` = complexes.")
    lines.append("- `DQ` = DockQ.")
    lines.append("- `IR` = iRMSD.")
    lines.append("- `Match no DQ/IR` = matched rows lacking numeric DockQ/iRMSD output.")
    lines.append("")
    lines.append("## Code Logic")
    lines.append("")
    lines.append("1. Match model files to benchmark rows using only `(PDB ID 1, PDB ID 2)` partner names from model filename, ignoring the first template/target token.")
    lines.append("2. Assign model partner chains to `PDB ID 1` and `PDB ID 2` roles, including swapped partner-name order handling.")
    lines.append("3. Expand each complex row to all native receptor-ligand chain pairs (e.g., `1AHW_AB:C` -> `AC` and `BC`).")
    lines.append("4. For each native pair, run DockQ and iRMSD against the downloaded native bound complex PDB.")
    lines.append("5. Run reciprocal iRMSD order as well, keep per-model best (`min`) across directions and native pairs.")
    lines.append("6. Aggregate by `(PDB ID 1, PDB ID 2, complex)` group: mean/variance/best across predictions/templates.")
    lines.append("")
    lines.append("## Key Insights")
    lines.append("")
    for x in infer_insights(data_by_set):
        lines.append("- %s" % x)
    lines.append("")

    lines.append("## Test Case (1AHW)")
    lines.append("")
    rigid_rows = data_by_set["rigid"]["matched_rows"]
    case_a = find_case_row(rigid_rows, "1ahwBC_1fgnH_0_1tfhA_0.rosetta_0001_0001.pdb")
    case_b = find_case_row(rigid_rows, "1ahwBC_1jptH_0_1tfhB_0.rosetta_0001_0001.pdb")
    if case_a:
        lines.append("- Benchmark row: `Complex=1AHW_AB:C`, `PDB ID 1=1FGN_LH`, `PDB ID 2=1TFH_A`.")
        lines.append("- Model A: `%s`" % case_a["filename"])
        lines.append("- Model chains used: receptor=`%s`, ligand=`%s`." % (case_a.get("model_receptor_chains"), case_a.get("model_ligand_chains")))
        lines.append("- Native pairs evaluated: `%s`." % (case_a.get("native_pairs_evaluated") or ""))
        lines.append("- DockQ per pair: `%s`." % (case_a.get("dockq_pair_scores") or ""))
        lines.append("- iRMSD forward per pair: `%s`." % (case_a.get("irmsd_pair_scores_forward") or ""))
        lines.append("- iRMSD reverse per pair: `%s`." % (case_a.get("irmsd_pair_scores_reverse") or ""))
        lines.append("- Best selected: DockQ=`%s` (%s), iRMSD=`%s`." % (case_a.get("dockq"), case_a.get("dockq_capri"), case_a.get("irmsd_min")))
        lines.append("")
        lines.append("Core logic used in this case:")
        lines.append("")
        lines.append("```python")
        lines.append("# complex row: 1AHW_AB:C -> native pairs: AC, BC")
        lines.append("native_pairs = [('A','C'), ('B','C')]")
        lines.append("# model 1ahwBC_1fgnH_0_1tfhA_0 -> model chains for PDB ID1/2 roles")
        lines.append("model_r, model_l = 'H', 'A'")
        lines.append("dockq_vals = []")
        lines.append("irmsd_vals = []")
        lines.append("for ref_r, ref_l in native_pairs:")
        lines.append("    dockq = DockQ(model, native_1ahw, mapping=f'{model_r}{model_l}:{ref_r}{ref_l}')")
        lines.append("    ir_fwd = iRMSD(model, model_r, model_l, native_1ahw, ref_r, ref_l)")
        lines.append("    ir_rev = iRMSD(model, model_l, model_r, native_1ahw, ref_l, ref_r)")
        lines.append("    dockq_vals.append(dockq)")
        lines.append("    irmsd_vals.extend([ir_fwd, ir_rev])")
        lines.append("best_dockq = max(dockq_vals)")
        lines.append("best_irmsd = min(irmsd_vals)")
        lines.append("```")
    if case_b:
        lines.append("- Model B (alt partner naming): `%s` -> model chains `%s/%s`, native pairs `%s`." % (
            case_b["filename"], case_b.get("model_receptor_chains"), case_b.get("model_ligand_chains"), case_b.get("native_pairs_evaluated")
        ))
    lines.append("")

    lines.append("## Plots")
    lines.append("")
    lines.append("Generated plot files:")
    for s in sets:
        lines.append("- `%s`" % str(plots_dir / ("%s_status_bar.png" % s)))
        lines.append("- `%s`" % str(plots_dir / ("%s_dockq_histogram.png" % s)))
        lines.append("- `%s`" % str(plots_dir / ("%s_irmsd_histogram.png" % s)))
        lines.append("- `%s`" % str(plots_dir / ("%s_dockq_vs_irmsd.png" % s)))
    lines.append("")
    lines.append("Plot confirmation notes:")
    lines.append("- Per-set status bars use `Matched`, `DockQ Scored`, `iRMSD Scored`, `DockQ Errors`, `iRMSD Errors` from each set summary.")
    lines.append("- Histograms and scatter use only rows with numeric metric values in `per_prediction.csv`.")
    lines.append("")

    for s in sets:
        d = data_by_set[s]
        lines.append("## %s Set Details" % s.title())
        lines.append("")
        lines.append("- Pair-summary rows: `%d`" % d["pair_rows_count"])
        if d["dockq_values"]:
            lines.append("- DockQ scored rows: `%d` (mean = `%s`, max = `%s`)" % (
                len(d["dockq_values"]), fmt(mean(d["dockq_values"])), fmt(max(d["dockq_values"]))
            ))
        else:
            lines.append("- DockQ scored rows: `0`")
        if d["irmsd_values"]:
            lines.append("- iRMSD scored rows: `%d` (mean = `%s`, min = `%s`)" % (
                len(d["irmsd_values"]), fmt(mean(d["irmsd_values"])), fmt(min(d["irmsd_values"]))
            ))
        else:
            lines.append("- iRMSD scored rows: `0`")

        if d["capri_counts"]:
            lines.append("- DockQ CAPRI classes: `%s`" % ", ".join("%s=%d" % (k, v) for k, v in sorted(d["capri_counts"].items())))
        if d["dockq_error_counts"]:
            lines.append("- DockQ failure reasons: `%s`" % ", ".join("%s=%d" % (k, v) for k, v in sorted(d["dockq_error_counts"].items())))
        if d["irmsd_error_counts"]:
            lines.append("- iRMSD failure reasons: `%s`" % ", ".join("%s=%d" % (k, v) for k, v in sorted(d["irmsd_error_counts"].items())))
        lines.append("")

        lines.append("Top 5 pairs by best DockQ:")
        if d["top_dockq_pairs"]:
            for val, r in d["top_dockq_pairs"]:
                lines.append("- `%s` (`%s`): best DockQ = `%s`, n = `%s`" % (r.get("pdb_pair_key") or r.get("target_token") or "NA", r["complex_raw"], fmt(val), r.get("n_predictions_total", "")))
        else:
            lines.append("- None")
        lines.append("")

        lines.append("Top 5 pairs by best (lowest) iRMSD:")
        if d["top_irmsd_pairs"]:
            for val, r in d["top_irmsd_pairs"]:
                lines.append("- `%s` (`%s`): best iRMSD = `%s`, n = `%s`" % (r.get("pdb_pair_key") or r.get("target_token") or "NA", r["complex_raw"], fmt(val), r.get("n_predictions_total", "")))
        else:
            lines.append("- None")
        lines.append("")

    lines.append("## Failure-Reason Notes")
    lines.append("")
    lines.append("- `DockQ` failures are mostly due to: `No interface found in native for selected mapping`.")
    lines.append("- This means the requested chain mapping exists syntactically, but DockQ cannot detect an interface in the native structure for that chain pair under its criteria.")
    lines.append("- `iRMSD` failures are fewer and usually stem from chain/alignment/runtime issues for specific rows.")
    lines.append("- Full error lists are exported to: `benchmark/prism_processed_results/benchmark_reports/error_files.csv`.")
    sample_err = None
    for s in sets:
        d = data_by_set[s]
        if d["dockq_error_files"]:
            e = d["dockq_error_files"][0]
            sample_err = (s, "dockq", e["filename"], e["complex_raw"], e["error"])
            break
        if d["irmsd_error_files"]:
            e = d["irmsd_error_files"][0]
            sample_err = (s, "irmsd", e["filename"], e["complex_raw"], e["error"])
            break
    if sample_err:
        lines.append("- Example error: set=`%s`, metric=`%s`, file=`%s`, complex=`%s`, reason=`%s`." % (
            sample_err[0], sample_err[1], sample_err[2], sample_err[3], sample_err[4].replace("|", "/")
        ))
    lines.append("")

    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(lines) + "\n")
    return md_path


def write_pdf(md_path, pdf_path, plots_dir, data_by_set):
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle

    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=letter,
        leftMargin=0.8 * inch,
        rightMargin=0.8 * inch,
        topMargin=0.8 * inch,
        bottomMargin=0.7 * inch,
    )
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="BodyStd",
        parent=styles["BodyText"],
        fontName="Helvetica",
        fontSize=10,
        leading=14,
        spaceAfter=6,
    ))
    styles.add(ParagraphStyle(
        name="Hdr1Std",
        parent=styles["Heading1"],
        fontName="Helvetica-Bold",
        fontSize=16,
        leading=20,
        textColor=colors.HexColor("#1F1F1F"),
        spaceAfter=10,
        spaceBefore=4,
    ))
    styles.add(ParagraphStyle(
        name="Hdr2Std",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=12,
        leading=16,
        textColor=colors.HexColor("#222222"),
        spaceAfter=8,
        spaceBefore=6,
    ))
    styles.add(ParagraphStyle(
        name="MonoSmall",
        parent=styles["BodyText"],
        fontName="Courier",
        fontSize=8.5,
        leading=11,
        spaceAfter=4,
    ))

    story = []
    story.append(Paragraph("PRISM Benchmark Analysis Report (Rigid / Medium / Difficult)", styles["Hdr1Std"]))
    story.append(Paragraph("Consolidated summary of PRISM predictions evaluated against T_Rigid, T_medium, and T_difficult benchmark CSVs using DockQ and iRMSD.", styles["BodyStd"]))
    story.append(Spacer(1, 6))

    # Bordered summary tables
    sets = ["rigid", "medium", "difficult"]
    story.append(Paragraph("Summary Table", styles["Hdr2Std"]))
    tdata = [["Set", "Cx CSV", "Models", "Match", "Unmatch", "DQ ok", "IR ok", "DQ err", "IR err"]]
    for s in sets:
        d = data_by_set[s]
        tdata.append([
            s.title(),
            str(d["complexes_in_csv"]) if d["complexes_in_csv"] is not None else "NA",
            str(d["total_models"]),
            str(d["matched_models"]),
            str(d["unmatched_models"]),
            str(d["dockq_scored"]),
            str(d["irmsd_scored"]),
            str(d["dockq_errors"]),
            str(d["irmsd_errors"]),
        ])
    tbl = Table(tdata, repeatRows=1, colWidths=[52, 44, 42, 42, 48, 42, 42, 45, 45])
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EAEAEA")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
        ("ALIGN", (1, 1), (-1, -1), "RIGHT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]))
    story.append(tbl)
    story.append(Spacer(1, 10))

    story.append(Paragraph("Complex-Level Coverage", styles["Hdr2Std"]))
    cdata = [["Set", "Cx w/Models", "Cx w/DQ", "Cx w/IR", "Match no DQ", "Match no IR"]]
    for s in sets:
        d = data_by_set[s]
        cdata.append([
            s.title(),
            str(d["complexes_with_models"]),
            str(d["complexes_with_dockq"]),
            str(d["complexes_with_irmsd"]),
            str(d["matched_without_dockq"]),
            str(d["matched_without_irmsd"]),
        ])
    ctbl = Table(cdata, repeatRows=1, colWidths=[54, 72, 62, 62, 80, 80])
    ctbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EAEAEA")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
        ("ALIGN", (1, 1), (-1, -1), "RIGHT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]))
    story.append(ctbl)
    story.append(Spacer(1, 10))

    story.append(Paragraph("Predicted Model Membership by Set", styles["Hdr2Std"]))
    total_models = data_by_set["rigid"]["total_models"]
    mdata = [["Set", "Models Belonging to Set", "Percent of All Predicted Models"]]
    for s in sets:
        d = data_by_set[s]
        pct = (100.0 * d["matched_models"] / total_models) if total_models else 0.0
        mdata.append([s.title(), str(d["matched_models"]), "%.2f%%" % pct])
    mtbl = Table(mdata, repeatRows=1, colWidths=[80, 120, 210])
    mtbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EAEAEA")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
        ("ALIGN", (1, 1), (-1, -1), "RIGHT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]))
    story.append(mtbl)
    story.append(Spacer(1, 10))

    # Insert key sections parsed from markdown (simple rendering)
    current_is_code = False
    for raw in md_path.read_text().splitlines():
        line = raw.strip()
        if line.startswith("```"):
            current_is_code = not current_is_code
            continue
        if not line:
            story.append(Spacer(1, 4))
            continue
        if line.startswith("# "):
            continue  # title already rendered
        if line.startswith("## "):
            story.append(Paragraph(line[3:], styles["Hdr2Std"]))
            continue
        if line.startswith("- "):
            story.append(Paragraph(u"• " + line[2:], styles["BodyStd"]))
            continue
        if line.startswith("|"):
            story.append(Paragraph(line.replace("&", "&amp;"), styles["MonoSmall"]))
            continue
        style = styles["MonoSmall"] if current_is_code else styles["BodyStd"]
        story.append(Paragraph(line.replace("&", "&amp;"), style))

    # Add plots as figures section with headers
    story.append(PageBreak())
    story.append(Paragraph("Figures", styles["Hdr1Std"]))
    for s in ["rigid", "medium", "difficult"]:
        story.append(Paragraph("%s Set Figures" % s.title(), styles["Hdr2Std"]))
        figure_specs = [
            ("%s Status" % s.title(), "%s_status_bar.png" % s),
            ("%s DockQ Distribution" % s.title(), "%s_dockq_histogram.png" % s),
            ("%s iRMSD Distribution" % s.title(), "%s_irmsd_histogram.png" % s),
            ("%s DockQ vs iRMSD" % s.title(), "%s_dockq_vs_irmsd.png" % s),
        ]
        for title, fname in figure_specs:
            img_path = plots_dir / fname
            if not img_path.exists():
                continue
            story.append(Paragraph(title, styles["BodyStd"]))
            img = Image(str(img_path))
            img._restrictSize(6.7 * inch, 4.8 * inch)
            story.append(img)
            story.append(Spacer(1, 6))

    doc.build(story)


def main():
    parser = argparse.ArgumentParser(description="Generate consolidated report and plots for PRISM benchmark analyses")
    parser.add_argument(
        "--results-root",
        default="benchmark/prism_processed_results",
        help="Root directory containing prism_<set>_analysis_all_jobs folders",
    )
    parser.add_argument(
        "--out-dir",
        default="benchmark/prism_processed_results/benchmark_reports",
        help="Output directory for consolidated report and plots",
    )
    parser.add_argument(
        "--data-dir",
        default="benchmark/data",
        help="Directory containing benchmark CSV files (T_Rigid.csv, T_medium.csv, T_difficult.csv)",
    )
    args = parser.parse_args()

    results_root = Path(args.results_root)
    out_dir = Path(args.out_dir)
    plots_dir = out_dir / "plots"
    data_dir = Path(args.data_dir)

    data_by_set = {
        "rigid": summarize_set("rigid", results_root, data_dir / "T_Rigid.csv"),
        "medium": summarize_set("medium", results_root, data_dir / "T_medium.csv"),
        "difficult": summarize_set("difficult", results_root, data_dir / "T_difficult.csv"),
    }

    make_plots(data_by_set, plots_dir)
    md_path = write_markdown(data_by_set, out_dir, plots_dir)
    err_csv = write_error_csv(data_by_set, out_dir)
    pdf_path = out_dir / "cross_benchmark_report.pdf"
    write_pdf(md_path, pdf_path, plots_dir, data_by_set)

    print("[INFO] Wrote", md_path)
    print("[INFO] Wrote", pdf_path)
    print("[INFO] Wrote", err_csv)
    print("[INFO] Plots in", plots_dir)


if __name__ == "__main__":
    raise SystemExit(main())
