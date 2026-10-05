#!/usr/bin/env python3
"""Diagnose the MultiProt alignment-to-transformation threshold gate."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--query-left", required=True)
    parser.add_argument("--query-right", required=True)
    parser.add_argument("--filter-mode", default="geometry_only_experimental")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    os.environ["PRISM_FILTER_MODE"] = args.filter_mode
    project_root = Path(__file__).resolve().parents[2]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    from src import transformation as tr

    template_root = args.root / "templates"
    alignment_root = args.root / "processed" / "alignment"
    templates = [
        line.strip()
        for line in (template_root / "calculated_templates.txt").read_text().splitlines()
        if line.strip()
    ]
    for template in templates:
        data = json.loads((template_root / "interfaces_lists" / f"{template}.json").read_text())
        for chain in template[4:]:
            tr.template_size[f"{template}_{chain}"] = len(data[chain])

    def load(query: str, template: str, chain: str) -> dict:
        path = alignment_root / f"{query}_{template}_{chain}.json"
        return json.loads(path.read_text())

    def side(template: str, chain: str, payload: dict) -> tuple[bool, list[str]]:
        reasons: list[str] = []
        if payload.get("status") != "success":
            reasons.append(str(payload.get("status", "missing")))
        if int(payload.get("match_count", 0)) < tr.MINIMUM_RESIDUE_MATCH_COUNT:
            reasons.append("match_count_below_minimum")
        if not tr.alignment_score_passes(payload):
            reasons.append(
                "multiprot_score_gate_failed"
                if str(payload.get("aligner", "")).strip().lower() == "multiprot"
                else "tm_score_below_threshold"
            )
        size = float(tr.template_size.get(f"{template}_{chain}", 0))
        if size > tr.TEMPLATE_RESIDUE_COUNT:
            if (int(payload.get("match_count", 0)) / size) * 100.0 <= (
                tr.MINIMUM_RESIDUE_MATCH_PERCENTAGE - tr.DIFF_PERCENTAGE
            ):
                reasons.append("match_percentage_below_threshold")
        elif size > 0 and (int(payload.get("match_count", 0)) / size) * 100.0 <= tr.MINIMUM_RESIDUE_MATCH_PERCENTAGE:
            reasons.append("match_percentage_below_threshold")
        return not reasons, reasons

    rows: list[dict[str, object]] = []
    for template in templates:
        c1, c2 = template[4], template[5]
        for orientation, left_chain, right_chain in (("o1", c1, c2), ("o2", c2, c1)):
            left = load(args.query_left, template, left_chain)
            right = load(args.query_right, template, right_chain)
            left_ok, left_reasons = side(template, left_chain, left)
            right_ok, right_reasons = side(template, right_chain, right)
            rows.append(
                {
                    "template": template,
                    "orientation": orientation,
                    "left_chain": left_chain,
                    "right_chain": right_chain,
                    "left_status": left.get("status", "missing"),
                    "right_status": right.get("status", "missing"),
                    "left_match_count": left.get("match_count", 0),
                    "right_match_count": right.get("match_count", 0),
                    "left_tm_score": left.get("tm_score", 0.0),
                    "right_tm_score": right.get("tm_score", 0.0),
                    "left_threshold_eligible": left_ok,
                    "right_threshold_eligible": right_ok,
                    "left_failure_reasons": ";".join(left_reasons),
                    "right_failure_reasons": ";".join(right_reasons),
                    "pair_threshold_eligible": left_ok and right_ok,
                }
            )

    both_success = sum(
        row["left_status"] == "success" and row["right_status"] == "success" for row in rows
    )
    summary = {
        "root": str(args.root.resolve()),
        "query_left": args.query_left,
        "query_right": args.query_right,
        "filter_mode": args.filter_mode,
        "template_count": len(templates),
        "orientation_count": len(rows),
        "both_success_count": both_success,
        "both_threshold_eligible_count": sum(row["pair_threshold_eligible"] for row in rows),
        "successful_side_count": sum(row["left_status"] == "success" for row in rows)
        + sum(row["right_status"] == "success" for row in rows),
        "successful_side_score_gate_count": sum(
            row["left_threshold_eligible"] for row in rows
        )
        + sum(row["right_threshold_eligible"] for row in rows),
        "interpretation": (
            "MultiProt records use their native match-count/coverage score "
            "contract; TMalign/GTalign records use the shared TM-score gate. "
            "The remaining pair eligibility is determined by these gates and "
            "the downstream transformation/clash checks."
        ),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with (args.output_dir / "multiprot_gate_ledger.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "multiprot_gate_diagnosis.json").write_text(
        json.dumps({"summary": summary, "orientations": rows}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (args.output_dir / "multiprot_gate_diagnosis.md").write_text(
        "# MultiProt alignment-to-transformation gate\n\n"
        + json.dumps(summary, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
