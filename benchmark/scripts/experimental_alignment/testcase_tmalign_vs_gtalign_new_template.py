#!/usr/bin/env python3
"""
Prepare and run a reproducible PRISM alignment backend test case on new_template.

What it does:
- Uses query surfaces from processed/surface_extraction/*.asa.pdb
- Uses template interfaces from new_template/templates/interfaces
- Runs TM-align in PRISM-like per-pair mode
- Runs GTalign in batch mode over the same query/reference set
- Writes JSONL outputs and a compact summary report

This is a benchmark/test helper and does not modify PRISM pipeline code.
"""

import argparse
import csv
import json
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class PairKey:
    query: str
    template: str
    chain: str


def run(cmd: List[str], timeout: Optional[int] = None, cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=timeout,
    )


def best_effort_version(binary: Path) -> str:
    probes = [
        [str(binary), "--version"],
        [str(binary), "-version"],
        [str(binary), "-v"],
        [str(binary), "-h"],
    ]
    for cmd in probes:
        try:
            p = run(cmd, timeout=10)
        except Exception:
            continue
        text = (p.stdout or "") + "\n" + (p.stderr or "")
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        if lines:
            return lines[0][:200]
    return "unknown"


def load_queries(inputs_csv: Path, surface_dir: Path, explicit_queries: Optional[List[str]]) -> List[str]:
    if explicit_queries:
        q = sorted(set(explicit_queries))
    else:
        q = []
        with open(inputs_csv, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                for col in ("Receptor", "Ligand"):
                    val = (row.get(col) or "").strip()
                    if len(val) == 5:
                        q.append(val)
        q = sorted(set(q))

    valid = []
    for protein in q:
        p = surface_dir / f"{protein}.asa.pdb"
        if p.exists():
            valid.append(protein)
    return valid


def load_templates(template_list: Path, interfaces_dir: Path, max_templates: int, explicit_templates: Optional[List[str]]) -> List[str]:
    if explicit_templates:
        templates = explicit_templates
    else:
        with open(template_list, "r") as f:
            templates = [ln.strip() for ln in f if ln.strip()]

    valid = []
    for t in templates:
        if len(t) < 6:
            continue
        chains = list(t[4:])
        ok = True
        for c in chains:
            if not (interfaces_dir / f"{t}_{c}_int.pdb").exists():
                ok = False
                break
        if ok:
            valid.append(t)
        if len(valid) >= max_templates:
            break
    return valid


def parse_tmalign_output(text: str) -> Tuple[int, float]:
    match_count = 0
    tm_scores: List[float] = []
    for line in text.splitlines():
        if line.startswith("Aligned length"):
            m = re.search(r"Aligned length\s*=\s*(\d+)", line)
            if m:
                match_count = int(m.group(1))
        elif line.startswith("TM-score"):
            m = re.search(r"TM-score\s*=\s*([0-9.]+)", line)
            if m:
                tm_scores.append(float(m.group(1)))
    return match_count, (max(tm_scores) if tm_scores else 0.0)


def run_tmalign_pairs(
    tmalign_bin: Path,
    queries: List[str],
    templates: List[str],
    surface_dir: Path,
    interfaces_dir: Path,
) -> Tuple[Dict[PairKey, Dict], float, int]:
    out: Dict[PairKey, Dict] = {}
    pair_count = 0
    t0 = time.perf_counter()
    for q in queries:
        qpath = surface_dir / f"{q}.asa.pdb"
        for t in templates:
            for c in t[4:]:
                rpath = interfaces_dir / f"{t}_{c}_int.pdb"
                pair = PairKey(q, t, c)
                cp = run([str(tmalign_bin), str(qpath), str(rpath)], timeout=180)
                pair_count += 1
                if cp.returncode != 0:
                    out[pair] = {
                        "ok": False,
                        "returncode": cp.returncode,
                        "match_count": 0,
                        "tm_score": 0.0,
                        "stderr": (cp.stderr or "")[:500],
                        "stdout": (cp.stdout or "")[:500],
                    }
                    continue
                match_count, tm_score = parse_tmalign_output(cp.stdout)
                out[pair] = {
                    "ok": True,
                    "returncode": 0,
                    "match_count": match_count,
                    "tm_score": tm_score,
                }
    dt = time.perf_counter() - t0
    return out, dt, pair_count


def _extract_floats(line: str) -> List[float]:
    return [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", line)]


def _extract_gtalign_alignment_seq(line: str, label: str) -> Optional[str]:
    pattern = r"^\s*" + re.escape(label) + r":\s+\d+\s+([A-Za-z\-]+)\s+\d+\s*$"
    m = re.match(pattern, line)
    if m:
        return m.group(1)
    toks = line.strip().split()
    if len(toks) >= 4:
        return toks[-2]
    return None


def _extract_gtalign_query_path(lines: List[str]) -> Optional[str]:
    for i, line in enumerate(lines):
        if line.startswith(" Query ("):
            j = i + 1
            while j < len(lines) and not lines[j].startswith(" Searched:"):
                candidate = lines[j].strip()
                if candidate and not candidate.startswith("Chn:"):
                    return candidate.split(" Chn:", 1)[0].strip()
                j += 1
    return None


def _parse_gtalign_hit_block(lines: List[str], start_idx: int) -> Tuple[Optional[Dict], int]:
    i = start_idx
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines) or not lines[i].lstrip().startswith(">"):
        return None, i + 1

    ref_path = lines[i].strip()[1:].split(" Chn:", 1)[0].strip()
    i += 1

    hit = {
        "ref_path": ref_path,
        "tm_ref": None,
        "tm_query": None,
        "aligned_length": None,
        "rotation_mat": None,
        "translation": None,
        "query_aln": "",
        "ref_aln": "",
    }
    q_chunks: List[str] = []
    r_chunks: List[str] = []

    while i < len(lines):
        line = lines[i]
        if re.match(r"^\s*\d+\.\s*$", line) or line.startswith("Query length:"):
            break

        if "TM-score (Refn./Query)" in line:
            m = re.search(r"TM-score \(Refn\./Query\)\s*=\s*([0-9.]+)\s*/\s*([0-9.]+)", line)
            if m:
                hit["tm_ref"] = float(m.group(1))
                hit["tm_query"] = float(m.group(2))

        if "Matched =" in line:
            m = re.search(r"Matched\s*=\s*(\d+)/", line)
            if m:
                hit["aligned_length"] = int(m.group(1))

        if line.lstrip().startswith("Query:") and i + 2 < len(lines) and lines[i + 2].lstrip().startswith("Refn.:"):
            q = _extract_gtalign_alignment_seq(line, "Query")
            r = _extract_gtalign_alignment_seq(lines[i + 2], "Refn.")
            if q is not None and r is not None:
                q_chunks.append(q)
                r_chunks.append(r)
                i += 3
                continue

        if line.strip().startswith("Rotation [3,3] and translation [3,1] for Query:"):
            rot: List[List[float]] = []
            trans: List[float] = []
            ok = True
            for j in range(1, 4):
                if i + j >= len(lines):
                    ok = False
                    break
                vals = _extract_floats(lines[i + j])
                if len(vals) < 4:
                    ok = False
                    break
                rot.append(vals[:3])
                trans.append(vals[3])
            if ok:
                hit["rotation_mat"] = rot
                hit["translation"] = trans
            i += 4
            continue

        i += 1

    hit["query_aln"] = "".join(q_chunks)
    hit["ref_aln"] = "".join(r_chunks)
    if hit["aligned_length"] is None and hit["query_aln"] and hit["ref_aln"]:
        hit["aligned_length"] = sum(1 for a, b in zip(hit["query_aln"], hit["ref_aln"]) if a != "-" and b != "-")

    return hit, i


def parse_gtalign_output_text(text: str) -> Tuple[str, List[Dict]]:
    lines = text.splitlines()
    query_path = _extract_gtalign_query_path(lines)
    if query_path is None:
        raise ValueError("Could not parse GTalign query path")

    hits: List[Dict] = []
    i = 0
    while i < len(lines):
        if re.match(r"^\s*\d+\.\s*$", lines[i]):
            hit, i = _parse_gtalign_hit_block(lines, i + 1)
            if hit:
                hits.append(hit)
            continue
        i += 1
    return query_path, hits


def run_gtalign_batch(
    gtalign_bin: Path,
    queries: List[str],
    templates: List[str],
    surface_dir: Path,
    interfaces_dir: Path,
    out_dir: Path,
    dev_min_length: int,
    pre_score: float,
    speed: int,
    refinement: int,
) -> Tuple[Dict[PairKey, Dict], float, int]:
    out_dir.mkdir(parents=True, exist_ok=True)
    qdir = out_dir / "queries"
    rdir = out_dir / "refs"
    gout = out_dir / "gtalign_out"
    qdir.mkdir(parents=True, exist_ok=True)
    rdir.mkdir(parents=True, exist_ok=True)
    gout.mkdir(parents=True, exist_ok=True)

    selected_query_paths: Dict[str, Path] = {}
    selected_ref_paths: Dict[Tuple[str, str], Path] = {}

    for q in queries:
        p = (surface_dir / f"{q}.asa.pdb").resolve()
        selected_query_paths[q] = p
        dst = qdir / p.name
        if not dst.exists():
            dst.symlink_to(p)

    for t in templates:
        for c in t[4:]:
            p = (interfaces_dir / f"{t}_{c}_int.pdb").resolve()
            selected_ref_paths[(t, c)] = p
            dst = rdir / p.name
            if not dst.exists():
                dst.symlink_to(p)

    cmd = [
        str(gtalign_bin),
        f"--qrs={qdir}",
        f"--rfs={rdir}",
        "-o",
        str(gout),
        "-s",
        "0",
        f"--nhits={len(selected_ref_paths)}",
        f"--nalns={len(selected_ref_paths)}",
        f"--dev-min-length={int(dev_min_length)}",
        f"--pre-score={float(pre_score)}",
        f"--speed={int(speed)}",
        f"--refinement={int(refinement)}",
    ]

    t0 = time.perf_counter()
    cp = run(cmd, timeout=3600)
    dt = time.perf_counter() - t0
    if cp.returncode != 0:
        raise RuntimeError(f"GTalign failed ({cp.returncode}): {(cp.stderr or cp.stdout)[:2000]}")

    out: Dict[PairKey, Dict] = {}
    query_basename_to_id = {f"{q}.asa.pdb": q for q in queries}
    ref_basename_to_key = {f"{t}_{c}_int.pdb": (t, c) for (t, c) in selected_ref_paths.keys()}

    for out_file in sorted(gout.glob("*.out")):
        query_path, hits = parse_gtalign_output_text(out_file.read_text(errors="replace"))
        qbase = Path(query_path).name
        if qbase not in query_basename_to_id:
            continue
        qid = query_basename_to_id[qbase]

        for hit in hits:
            rbase = Path(hit["ref_path"]).name
            if rbase not in ref_basename_to_key:
                continue
            t, c = ref_basename_to_key[rbase]
            tm_candidates = [v for v in (hit.get("tm_ref"), hit.get("tm_query")) if isinstance(v, (float, int))]
            out[PairKey(qid, t, c)] = {
                "ok": True,
                "match_count": int(hit.get("aligned_length") or 0),
                "tm_score": float(max(tm_candidates) if tm_candidates else 0.0),
            }

    expected_pairs = len(queries) * sum(len(t[4:]) for t in templates)
    return out, dt, expected_pairs


def write_jsonl(path: Path, rows: List[Dict]) -> None:
    with open(path, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")


def summarize(
    tm: Dict[PairKey, Dict],
    gt: Dict[PairKey, Dict],
    tm_seconds: float,
    gt_seconds: float,
    expected_pairs: int,
) -> Dict:
    keys = sorted(set(tm.keys()) | set(gt.keys()), key=lambda x: (x.query, x.template, x.chain))

    both = 0
    only_tm = 0
    only_gt = 0
    tm_positive = 0
    gt_positive = 0
    tm_diffs: List[float] = []
    match_diffs: List[int] = []

    rows = []
    for k in keys:
        t = tm.get(k)
        g = gt.get(k)
        if t and t.get("tm_score", 0.0) > 0:
            tm_positive += 1
        if g and g.get("tm_score", 0.0) > 0:
            gt_positive += 1

        if t and g:
            both += 1
            tm_diff = abs(float(t.get("tm_score", 0.0)) - float(g.get("tm_score", 0.0)))
            match_diff = abs(int(t.get("match_count", 0)) - int(g.get("match_count", 0)))
            tm_diffs.append(tm_diff)
            match_diffs.append(match_diff)
        elif t and not g:
            only_tm += 1
        elif g and not t:
            only_gt += 1

        rows.append(
            {
                "query": k.query,
                "template": k.template,
                "chain": k.chain,
                "tmalign_tm_score": (t or {}).get("tm_score", 0.0),
                "gtalign_tm_score": (g or {}).get("tm_score", 0.0),
                "tmalign_match_count": (t or {}).get("match_count", 0),
                "gtalign_match_count": (g or {}).get("match_count", 0),
                "tm_score_abs_diff": abs(float((t or {}).get("tm_score", 0.0)) - float((g or {}).get("tm_score", 0.0))),
                "match_count_abs_diff": abs(int((t or {}).get("match_count", 0)) - int((g or {}).get("match_count", 0))),
                "present_in_tmalign": bool(t),
                "present_in_gtalign": bool(g),
            }
        )

    tm_diffs_sorted = sorted(tm_diffs)
    match_diffs_sorted = sorted(match_diffs)

    def pct(vals: List[float], p: float) -> float:
        if not vals:
            return 0.0
        idx = int(round((len(vals) - 1) * p))
        return float(vals[max(0, min(len(vals) - 1, idx))])

    return {
        "expected_pairs": expected_pairs,
        "tmalign_pairs": len(tm),
        "gtalign_pairs": len(gt),
        "both_pairs": both,
        "only_tmalign_pairs": only_tm,
        "only_gtalign_pairs": only_gt,
        "tmalign_positive_tm_score_pairs": tm_positive,
        "gtalign_positive_tm_score_pairs": gt_positive,
        "tmalign_runtime_seconds": tm_seconds,
        "gtalign_runtime_seconds": gt_seconds,
        "tmalign_pairs_per_second": (len(tm) / tm_seconds) if tm_seconds > 0 else 0.0,
        "gtalign_effective_pairs_per_second": (expected_pairs / gt_seconds) if gt_seconds > 0 else 0.0,
        "tm_score_abs_diff_mean": (sum(tm_diffs) / len(tm_diffs)) if tm_diffs else 0.0,
        "tm_score_abs_diff_p50": pct(tm_diffs_sorted, 0.5),
        "tm_score_abs_diff_p90": pct(tm_diffs_sorted, 0.9),
        "match_count_abs_diff_mean": (sum(match_diffs) / len(match_diffs)) if match_diffs else 0.0,
        "match_count_abs_diff_p50": pct(match_diffs_sorted, 0.5),
        "match_count_abs_diff_p90": pct(match_diffs_sorted, 0.9),
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="PRISM test case: TM-align vs GTalign on new_template")
    parser.add_argument("--repo_root", default=".", help="Path to PRISM repository root")
    parser.add_argument("--inputs_csv", default="inputs.csv")
    parser.add_argument("--surface_dir", default="processed/surface_extraction")
    parser.add_argument("--template_root", default="new_template/templates")
    parser.add_argument("--template_list", default=None, help="Defaults to <template_root>/calculated_templates.txt")
    parser.add_argument("--max_templates", type=int, default=25, help="How many templates to include")
    parser.add_argument("--queries", nargs="*", default=None, help="Optional explicit query IDs")
    parser.add_argument("--templates", nargs="*", default=None, help="Optional explicit template IDs")

    parser.add_argument("--tmalign", default="external_tools/TMalign", help="Path to TMalign binary")
    parser.add_argument("--gtalign", required=True, help="Path to GTalign binary")
    parser.add_argument("--gtalign_dev_min_length", type=int, default=3)
    parser.add_argument("--gtalign_pre_score", type=float, default=0.0)
    parser.add_argument("--gtalign_speed", type=int, default=0)
    parser.add_argument("--gtalign_refinement", type=int, default=3)

    parser.add_argument(
        "--output_dir",
        default="benchmark/prism_processed_results/testcase_tmalign_vs_gtalign_new_template",
        help="Directory for test-case outputs",
    )
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    inputs_csv = (root / args.inputs_csv).resolve()
    surface_dir = (root / args.surface_dir).resolve()
    template_root = (root / args.template_root).resolve()
    template_list = (root / args.template_list).resolve() if args.template_list else (template_root / "calculated_templates.txt").resolve()
    interfaces_dir = (template_root / "interfaces").resolve()

    tmalign_bin = (root / args.tmalign).resolve() if not Path(args.tmalign).is_absolute() else Path(args.tmalign)
    gtalign_bin = (root / args.gtalign).resolve() if not Path(args.gtalign).is_absolute() else Path(args.gtalign)

    if not tmalign_bin.exists():
        raise FileNotFoundError(f"TMalign binary not found: {tmalign_bin}")
    if not gtalign_bin.exists():
        raise FileNotFoundError(f"GTalign binary not found: {gtalign_bin}")
    if not inputs_csv.exists():
        raise FileNotFoundError(f"inputs.csv not found: {inputs_csv}")
    if not surface_dir.exists():
        raise FileNotFoundError(f"surface dir not found: {surface_dir}")
    if not interfaces_dir.exists():
        raise FileNotFoundError(f"interfaces dir not found: {interfaces_dir}")
    if not template_list.exists():
        raise FileNotFoundError(f"template list not found: {template_list}")

    queries = load_queries(inputs_csv, surface_dir, args.queries)
    templates = load_templates(template_list, interfaces_dir, args.max_templates, args.templates)

    if not queries:
        raise RuntimeError("No valid queries found (.asa.pdb missing for selected IDs).")
    if not templates:
        raise RuntimeError("No valid templates found with complete chain interface files.")

    output_dir = (root / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Queries ({len(queries)}): {queries}")
    print(f"Templates ({len(templates)}): {templates[:5]}{' ...' if len(templates) > 5 else ''}")
    print(f"Total expected pairs: {len(queries) * sum(len(t[4:]) for t in templates)}")

    print("Running TM-align per-pair...")
    tm_res, tm_seconds, tm_pairs = run_tmalign_pairs(
        tmalign_bin,
        queries,
        templates,
        surface_dir,
        interfaces_dir,
    )
    print(f"TM-align pairs={tm_pairs}, runtime={tm_seconds:.3f}s")

    print("Running GTalign batch...")
    gt_res, gt_seconds, expected_pairs = run_gtalign_batch(
        gtalign_bin,
        queries,
        templates,
        surface_dir,
        interfaces_dir,
        out_dir=output_dir / "tmp",
        dev_min_length=args.gtalign_dev_min_length,
        pre_score=args.gtalign_pre_score,
        speed=args.gtalign_speed,
        refinement=args.gtalign_refinement,
    )
    print(f"GTalign parsed_pairs={len(gt_res)}, runtime={gt_seconds:.3f}s")

    summary = summarize(tm_res, gt_res, tm_seconds, gt_seconds, expected_pairs)
    summary["metadata"] = {
        "repo_root": str(root),
        "template_root": str(template_root),
        "template_list": str(template_list),
        "surface_dir": str(surface_dir),
        "queries": queries,
        "templates": templates,
        "tmalign_bin": str(tmalign_bin),
        "gtalign_bin": str(gtalign_bin),
        "tmalign_version": best_effort_version(tmalign_bin),
        "gtalign_version": best_effort_version(gtalign_bin),
        "gtalign_params": {
            "dev_min_length": args.gtalign_dev_min_length,
            "pre_score": args.gtalign_pre_score,
            "speed": args.gtalign_speed,
            "refinement": args.gtalign_refinement,
        },
    }

    tm_rows = []
    for k, v in sorted(tm_res.items(), key=lambda kv: (kv[0].query, kv[0].template, kv[0].chain)):
        tm_rows.append({"query": k.query, "template": k.template, "chain": k.chain, **v})
    gt_rows = []
    for k, v in sorted(gt_res.items(), key=lambda kv: (kv[0].query, kv[0].template, kv[0].chain)):
        gt_rows.append({"query": k.query, "template": k.template, "chain": k.chain, **v})

    write_jsonl(output_dir / "tmalign_results.jsonl", tm_rows)
    write_jsonl(output_dir / "gtalign_results.jsonl", gt_rows)
    write_jsonl(output_dir / "pairwise_comparison.jsonl", summary["rows"])

    summary_no_rows = dict(summary)
    summary_no_rows.pop("rows", None)
    with open(output_dir / "summary.json", "w") as f:
        json.dump(summary_no_rows, f, indent=2)

    with open(output_dir / "README_test_case.txt", "w") as f:
        f.write("PRISM test case: TM-align vs GTalign on new_template\n")
        f.write(f"queries={len(queries)}, templates={len(templates)}, expected_pairs={expected_pairs}\n")
        f.write(f"tmalign_runtime_seconds={tm_seconds:.6f}\n")
        f.write(f"gtalign_runtime_seconds={gt_seconds:.6f}\n")
        f.write(f"tm_score_abs_diff_mean={summary_no_rows['tm_score_abs_diff_mean']:.6f}\n")
        f.write(f"tm_score_abs_diff_p90={summary_no_rows['tm_score_abs_diff_p90']:.6f}\n")

    print("Done. Outputs:")
    print(f"- {output_dir / 'summary.json'}")
    print(f"- {output_dir / 'pairwise_comparison.jsonl'}")
    print(f"- {output_dir / 'tmalign_results.jsonl'}")
    print(f"- {output_dir / 'gtalign_results.jsonl'}")


if __name__ == "__main__":
    main()
