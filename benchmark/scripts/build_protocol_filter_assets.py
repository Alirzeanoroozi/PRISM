#!/usr/bin/env python3
"""Convert retained legacy hotspot/contact text into derived JSON assets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def parse_legacy_hotspots(path: Path) -> list[tuple[str, str, str]]:
    records = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 2 and not line.lstrip().startswith("#"):
            residue = parts[0].split(".")
            if len(residue) >= 3:
                records.append((residue[0], residue[1], residue[2]))
    return records


def parse_legacy_contacts(path: Path) -> list[tuple[str, str]]:
    records = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 2 and not line.lstrip().startswith("#"):
            records.append((parts[0], parts[1]))
    return records


def build_assets(template_list: Path, legacy_root: Path, output_root: Path) -> dict[str, int]:
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"refusing nonempty derived asset root: {output_root}")
    (output_root / "hotspots").mkdir(parents=True, exist_ok=True)
    (output_root / "contacts").mkdir(parents=True, exist_ok=True)
    template_ids = sorted({line.split()[0] for line in template_list.read_text(encoding="utf-8").splitlines() if line.strip() and not line.lstrip().startswith("#")})
    hotspot_count = contact_count = 0
    manifest_rows = []
    for template_id in template_ids:
        hotspot = legacy_root / "hotspot" / f"hotspot{template_id}"
        contact = legacy_root / "contact" / f"{template_id}.txt"
        if not hotspot.is_file() or not contact.is_file():
            raise RuntimeError(f"missing legacy asset for {template_id}")
        hotspot_path = output_root / "hotspots" / f"{template_id}.json"
        contact_path = output_root / "contacts" / f"{template_id}.json"
        hotspot_path.write_text(json.dumps(parse_legacy_hotspots(hotspot), sort_keys=True) + "\n", encoding="utf-8")
        contact_path.write_text(json.dumps(parse_legacy_contacts(contact), sort_keys=True) + "\n", encoding="utf-8")
        for asset_type, path in (("hotspots", hotspot_path), ("contacts", contact_path)):
            manifest_rows.append({"template_id": template_id, "asset_type": asset_type, "path": str(path.relative_to(output_root)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "size": path.stat().st_size})
        hotspot_count += 1
        contact_count += 1
    with (output_root / "asset_manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("template_id", "asset_type", "path", "sha256", "size"), delimiter="\t")
        writer.writeheader(); writer.writerows(manifest_rows)
    return {"template_count": len(template_ids), "hotspot_count": hotspot_count, "contact_count": contact_count, "asset_manifest": str(output_root / "asset_manifest.tsv")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--template-list", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_assets(args.template_list, args.legacy_root, args.output_root), sort_keys=True))
