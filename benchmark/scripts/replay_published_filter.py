#!/usr/bin/env python3
"""Replay published hotspot/contact filtering on retained alignment JSONs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.pdb_download import normalize_target_id
from src.template_filtering import evaluate_protocol_candidate, load_filter_assets


FIELDS = ("row_index", "receptor", "ligand", "template_id", "orientation", "status", "reason", "hotspot_matches", "complementary_contacts", "hotspots_sha256", "contacts_sha256")


def _template_ids(path: Path) -> list[str]:
    return sorted({line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.lstrip().startswith("#")})


def _alignment_index(root: Path, template_ids: list[str]) -> dict[tuple[str, str, str], dict]:
    index: dict[tuple[str, str, str], dict] = {}
    for path in sorted(root.rglob("*.json")):
        stem = path.stem
        for template_id in template_ids:
            for chain in template_id[4:]:
                suffix = f"_{template_id}_{chain}"
                if stem.endswith(suffix):
                    query = stem[: -len(suffix)]
                    try:
                        index[(normalize_target_id(query), template_id, chain)] = json.loads(path.read_text(encoding="utf-8"))
                    except (OSError, json.JSONDecodeError):
                        index[(normalize_target_id(query), template_id, chain)] = {}
                    break
            else:
                continue
            break
    return index


def replay(inputs_csv: Path, template_list: Path, alignment_root: Path, asset_root: Path, output: Path) -> list[dict[str, str]]:
    templates = _template_ids(template_list)
    alignments = _alignment_index(alignment_root, templates)
    with inputs_csv.open(newline="", encoding="utf-8") as handle:
        inputs = list(csv.DictReader(handle))
    rows: list[dict[str, str]] = []
    for row_index, row in enumerate(inputs, 1):
        receptor, ligand = normalize_target_id(row["Receptor"]), normalize_target_id(row["Ligand"])
        for template_id in templates:
            chains = template_id[4:]
            try:
                assets = load_filter_assets(template_id, asset_root)
            except (FileNotFoundError, ValueError) as exc:
                for orientation in ("o1", "o2"):
                    rows.append({"row_index": str(row_index), "receptor": receptor, "ligand": ligand, "template_id": template_id, "orientation": orientation, "status": "failed", "reason": f"protocol_assets:{exc}", "hotspot_matches": "0", "complementary_contacts": "0", "hotspots_sha256": "", "contacts_sha256": ""})
                continue
            for orientation, left_query, right_query, left_chain, right_chain in (("o1", receptor, ligand, chains[0], chains[1]), ("o2", receptor, ligand, chains[1], chains[0])):
                left = alignments.get((left_query, template_id, left_chain))
                right = alignments.get((right_query, template_id, right_chain))
                if not left or not right:
                    decision = None
                    status, reason = "failed", "alignment_missing"
                else:
                    decision = evaluate_protocol_candidate(left.get("match_dict", {}), right.get("match_dict", {}), assets["hotspots"], assets["hotspots"], assets["contacts"])
                    status, reason = ("passed", decision.reason) if decision.passed else ("failed", decision.reason)
                rows.append({
                    "row_index": str(row_index), "receptor": receptor, "ligand": ligand, "template_id": template_id, "orientation": orientation,
                    "status": status, "reason": reason,
                    "hotspot_matches": str(decision.hotspot_matches if decision else 0), "complementary_contacts": str(decision.complementary_contacts if decision else 0),
                    "hotspots_sha256": assets["hotspots_sha256"], "contacts_sha256": assets["contacts_sha256"],
                })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--template-list", type=Path, required=True)
    parser.add_argument("--alignment-root", type=Path, required=True)
    parser.add_argument("--asset-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = replay(args.inputs, args.template_list, args.alignment_root, args.asset_root, args.output)
    print(json.dumps({"rows": len(rows), "passed": sum(row["status"] == "passed" for row in rows)}, sort_keys=True))
