#!/usr/bin/env python3
"""Compare modern JSON filter assets with derived legacy protocol assets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


FIELDS = (
    "template_id",
    "current_hotspots_status",
    "current_contacts_status",
    "legacy_hotspot_count",
    "current_hotspot_count",
    "legacy_contact_count",
    "current_contact_count",
    "hotspots_exact",
    "contacts_numeric_exact",
    "parity_status",
    "current_hotspots_sha256",
    "current_contacts_sha256",
    "legacy_hotspots_sha256",
    "legacy_contacts_sha256",
    "difference",
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else ""


def _hotspots(value: object) -> set[tuple[str, str, int]]:
    records: set[tuple[str, str, int]] = set()
    if isinstance(value, dict):
        for chain, entries in value.items():
            for entry in entries if isinstance(entries, list) else []:
                if isinstance(entry, (list, tuple)) and len(entry) >= 2:
                    records.add((str(chain), str(entry[1]).upper(), int(entry[0])))
    elif isinstance(value, list):
        for entry in value:
            if isinstance(entry, (list, tuple)) and len(entry) >= 3:
                records.add((str(entry[0]), str(entry[1]).upper(), int(entry[2])))
    return records


def _legacy_contact(value: object) -> set[tuple[str, int, str, int]]:
    records: set[tuple[str, int, str, int]] = set()
    if not isinstance(value, list):
        return records
    for pair in value:
        if not isinstance(pair, (list, tuple)) or len(pair) < 2:
            continue
        endpoints = []
        for token in pair[:2]:
            parts = str(token).split(".")
            if len(parts) < 3:
                continue
            endpoints.append((parts[0], int(parts[2])))
        if len(endpoints) == 2:
            records.add((endpoints[0][0], endpoints[0][1], endpoints[1][0], endpoints[1][1]))
    return records


def _current_contact(value: object, chains: str) -> set[tuple[str, int, str, int]]:
    records: set[tuple[str, int, str, int]] = set()
    if not isinstance(value, list) or len(chains) < 2:
        return records
    left, right = chains[0], chains[1]
    for pair in value:
        if isinstance(pair, (list, tuple)) and len(pair) >= 2:
            records.add((left, int(pair[0]), right, int(pair[1])))
    return records


def compare_assets(template_list: Path, modern_root: Path, legacy_root: Path, output: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    ids = sorted({line.strip() for line in template_list.read_text(encoding="utf-8").splitlines() if line.strip() and not line.lstrip().startswith("#")})
    for template_id in ids:
        modern_hotspot = modern_root / "hotspots" / f"{template_id}.json"
        modern_contact = modern_root / "contacts" / f"{template_id}.json"
        legacy_hotspot = legacy_root / "hotspots" / f"{template_id}.json"
        legacy_contact = legacy_root / "contacts" / f"{template_id}.json"
        difference: list[str] = []
        try:
            current_h = _hotspots(json.loads(modern_hotspot.read_text(encoding="utf-8"))) if modern_hotspot.is_file() else set()
            current_h_status = "present" if modern_hotspot.is_file() else "missing"
        except (OSError, ValueError, TypeError):
            current_h, current_h_status = set(), "invalid"
        try:
            current_c = _current_contact(json.loads(modern_contact.read_text(encoding="utf-8")), template_id[4:]) if modern_contact.is_file() else set()
            current_c_status = "present" if modern_contact.is_file() else "missing"
        except (OSError, ValueError, TypeError):
            current_c, current_c_status = set(), "invalid"
        legacy_h = _hotspots(json.loads(legacy_hotspot.read_text(encoding="utf-8"))) if legacy_hotspot.is_file() else set()
        legacy_c = _legacy_contact(json.loads(legacy_contact.read_text(encoding="utf-8"))) if legacy_contact.is_file() else set()
        hotspots_exact = current_h == legacy_h and current_h_status == "present"
        contacts_exact = current_c == legacy_c and current_c_status == "present"
        if current_h_status != "present":
            difference.append(f"hotspots_{current_h_status}")
        elif not hotspots_exact:
            difference.append("hotspots_semantic_disagreement")
        if current_c_status != "present":
            difference.append(f"contacts_{current_c_status}")
        elif not contacts_exact:
            difference.append("contacts_numeric_disagreement")
        status = "exact" if not difference else ("current_assets_missing" if any(item.endswith("missing") for item in difference) else "disagreement")
        rows.append({
            "template_id": template_id,
            "current_hotspots_status": current_h_status,
            "current_contacts_status": current_c_status,
            "legacy_hotspot_count": str(len(legacy_h)),
            "current_hotspot_count": str(len(current_h)),
            "legacy_contact_count": str(len(legacy_c)),
            "current_contact_count": str(len(current_c)),
            "hotspots_exact": str(int(hotspots_exact)),
            "contacts_numeric_exact": str(int(contacts_exact)),
            "parity_status": status,
            "current_hotspots_sha256": _sha(modern_hotspot),
            "current_contacts_sha256": _sha(modern_contact),
            "legacy_hotspots_sha256": _sha(legacy_hotspot),
            "legacy_contacts_sha256": _sha(legacy_contact),
            "difference": ";".join(difference),
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--template-list", type=Path, required=True)
    parser.add_argument("--modern-root", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = compare_assets(args.template_list, args.modern_root, args.legacy_root, args.output)
    counts = {status: sum(row["parity_status"] == status for row in rows) for status in ("exact", "disagreement", "current_assets_missing")}
    print(json.dumps({"template_count": len(rows), **counts}, sort_keys=True))
