#!/usr/bin/env python3
"""Run a matched, diagnostic-only TMalign calibration against a MultiProt panel."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import sys
from pathlib import Path


def _read_panel(panel_root: Path) -> tuple[list[str], list[str]]:
    with (panel_root / "inputs.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    queries = sorted({row["Receptor"] for row in rows} | {row["Ligand"] for row in rows})
    templates = [
        line.strip()
        for line in (panel_root / "templates" / "calculated_templates.txt").read_text().splitlines()
        if line.strip()
    ]
    if not queries or not templates:
        raise RuntimeError("panel is missing query or template records")
    return queries, templates


def _prepare_workdir(panel_root: Path, output_root: Path) -> None:
    (output_root / "processed").mkdir(parents=True, exist_ok=True)
    for source, destination in (
        (panel_root / "processed" / "surface_extraction", output_root / "processed" / "surface_extraction"),
        (panel_root / "templates", output_root / "templates"),
    ):
        if destination.exists() or destination.is_symlink():
            if destination.is_symlink() and destination.resolve() == source.resolve():
                continue
            raise RuntimeError(f"refusing to replace existing path: {destination}")
        destination.symlink_to(source, target_is_directory=True)


def _load_json(root: Path, query: str, template: str, chain: str) -> dict:
    path = root / "processed" / "alignment" / f"{query}_{template}_{chain}.json"
    if not path.is_file():
        return {"status": "missing", "match_count": 0, "tm_score": 0.0}
    return json.loads(path.read_text())


def _side_gate(payload: dict, template: str, chain: str, template_size: dict[str, int]) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    if payload.get("status") != "success":
        reasons.append(str(payload.get("status", "missing")))
    if int(payload.get("match_count", 0)) < 15:
        reasons.append("match_count_below_minimum")
    if float(payload.get("tm_score", 0.0)) < 0.5:
        reasons.append("tm_score_below_threshold")
    size = float(template_size.get(f"{template}_{chain}", 0))
    count = int(payload.get("match_count", 0))
    if size > 150:
        if (count / size) * 100.0 <= 30.0:
            reasons.append("match_percentage_below_threshold")
    elif size > 0 and (count / size) * 100.0 <= 50.0:
        reasons.append("match_percentage_below_threshold")
    return not reasons, reasons


def _pearson(values_x: list[float], values_y: list[float]) -> float | None:
    if len(values_x) < 2:
        return None
    mean_x = sum(values_x) / len(values_x)
    mean_y = sum(values_y) / len(values_y)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in zip(values_x, values_y))
    denom_x = math.sqrt(sum((x - mean_x) ** 2 for x in values_x))
    denom_y = math.sqrt(sum((y - mean_y) ** 2 for y in values_y))
    if denom_x == 0.0 or denom_y == 0.0:
        return None
    return numerator / (denom_x * denom_y)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[2]
    panel_root = args.panel_root.resolve()
    output_root = args.output_root.resolve()
    if output_root == panel_root or panel_root in output_root.parents:
        raise RuntimeError("output root must be isolated from the retained panel")
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"output root is not empty: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    _prepare_workdir(panel_root, output_root)

    # The alignment module uses project-relative processed paths, so import it
    # only after entering the isolated output root.
    os.chdir(output_root)
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from src.alignment import align as align_tmalign

    queries, templates = _read_panel(panel_root)
    os.environ["PRISM_TMALIGN"] = str(project_root / "external_tools" / "TMalign")
    os.environ["PRISM_TMALIGN_WORKERS"] = str(args.workers)
    os.environ["PRISM_STAGE_STATUS_PATH"] = str(output_root / "stage_status.jsonl")
    align_tmalign(queries, templates)

    template_size: dict[str, int] = {}
    for template in templates:
        data = json.loads((output_root / "templates" / "interfaces_lists" / f"{template}.json").read_text())
        for chain in template[4:]:
            template_size[f"{template}_{chain}"] = len(data[chain])

    rows: list[dict[str, object]] = []
    paired_mp_scores: list[float] = []
    paired_tm_scores: list[float] = []
    for query in queries:
        for template in templates:
            chains = template[4:]
            for chain in chains:
                mp = _load_json(panel_root, query, template, chain)
                tm = _load_json(output_root, query, template, chain)
                mp_ok, mp_reasons = _side_gate(mp, template, chain, template_size)
                tm_ok, tm_reasons = _side_gate(tm, template, chain, template_size)
                if mp.get("status") == "success" and tm.get("status") == "success":
                    paired_mp_scores.append(float(mp.get("tm_score", 0.0)))
                    paired_tm_scores.append(float(tm.get("tm_score", 0.0)))
                rows.append(
                    {
                        "query": query,
                        "template": template,
                        "chain": chain,
                        "multiprot_status": mp.get("status", "missing"),
                        "tmalign_status": tm.get("status", "missing"),
                        "multiprot_match_count": mp.get("match_count", 0),
                        "tmalign_match_count": tm.get("match_count", 0),
                        "multiprot_tm_score": mp.get("tm_score", 0.0),
                        "tmalign_tm_score": tm.get("tm_score", 0.0),
                        "multiprot_gate": mp_ok,
                        "tmalign_gate": tm_ok,
                        "multiprot_gate_reasons": ";".join(mp_reasons),
                        "tmalign_gate_reasons": ";".join(tm_reasons),
                    }
                )

    pair_rows: list[dict[str, object]] = []
    for template in templates:
        c1, c2 = template[4], template[5]
        for orientation, left_chain, right_chain in (("o1", c1, c2), ("o2", c2, c1)):
            left_mp = _load_json(panel_root, queries[0], template, left_chain)
            right_mp = _load_json(panel_root, queries[1], template, right_chain)
            left_tm = _load_json(output_root, queries[0], template, left_chain)
            right_tm = _load_json(output_root, queries[1], template, right_chain)
            left_mp_ok, _ = _side_gate(left_mp, template, left_chain, template_size)
            right_mp_ok, _ = _side_gate(right_mp, template, right_chain, template_size)
            left_tm_ok, _ = _side_gate(left_tm, template, left_chain, template_size)
            right_tm_ok, _ = _side_gate(right_tm, template, right_chain, template_size)
            pair_rows.append(
                {
                    "template": template,
                    "orientation": orientation,
                    "multiprot_both_success": left_mp.get("status") == "success" and right_mp.get("status") == "success",
                    "tmalign_both_success": left_tm.get("status") == "success" and right_tm.get("status") == "success",
                    "multiprot_pair_gate": left_mp_ok and right_mp_ok,
                    "tmalign_pair_gate": left_tm_ok and right_tm_ok,
                }
            )

    summary = {
        "panel_root": str(panel_root),
        "output_root": str(output_root),
        "queries": queries,
        "template_count": len(templates),
        "record_count": len(rows),
        "multiprot_success_count": sum(row["multiprot_status"] == "success" for row in rows),
        "tmalign_success_count": sum(row["tmalign_status"] == "success" for row in rows),
        "paired_success_count": len(paired_mp_scores),
        "multiprot_pair_both_success_count": sum(row["multiprot_both_success"] for row in pair_rows),
        "tmalign_pair_both_success_count": sum(row["tmalign_both_success"] for row in pair_rows),
        "multiprot_pair_gate_count": sum(row["multiprot_pair_gate"] for row in pair_rows),
        "tmalign_pair_gate_count": sum(row["tmalign_pair_gate"] for row in pair_rows),
        "score_pearson_multiprot_vs_tmalign": _pearson(paired_mp_scores, paired_tm_scores),
        "score_median_abs_difference": (
            sorted(abs(x - y) for x, y in zip(paired_mp_scores, paired_tm_scores))[len(paired_mp_scores) // 2]
            if paired_mp_scores
            else None
        ),
        "interpretation": (
            "Diagnostic calibration only: this compares the retained MultiProt records "
            "with TMalign on identical query/template/chain inputs. It does not change "
            "the production threshold or establish structural-quality equivalence."
        ),
    }
    fields = list(rows[0])
    with (output_root / "matched_alignment_ledger.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with (output_root / "matched_pair_ledger.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)
    (output_root / "matched_alignment_calibration.json").write_text(
        json.dumps({"summary": summary, "pairs": pair_rows}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
