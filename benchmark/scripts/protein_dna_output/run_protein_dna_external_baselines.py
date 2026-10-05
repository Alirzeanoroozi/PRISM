#!/usr/bin/env python3
"""
Probe and stage state-of-the-art protein-DNA baseline tools alongside the PRISM benchmark.

This script does not modify PRISM core behavior. It records which external baselines are
available on the current system, which are web-only, and which benchmark pairs would be
eligible for future execution once the relevant binaries or services are wired in.
"""

import argparse
import csv
import json
import shutil
from pathlib import Path


def load_manifest(manifest_path):
    with open(manifest_path, newline="") as handle:
        return list(csv.DictReader(handle))


def load_tool_registry(registry_path):
    with open(registry_path, "r") as handle:
        return json.load(handle)


def probe_tool(tool):
    available = None
    for candidate in tool.get("binary_candidates", []):
        path = shutil.which(candidate)
        if path:
            available = path
            break
    status = "available" if available else "skipped"
    reason = "" if available else ("web_server_only" if tool.get("integration_kind") == "web_server" else "tool_not_installed")
    return {
        "tool_id": tool["tool_id"],
        "tool_name": tool["name"],
        "integration_kind": tool["integration_kind"],
        "supports": "|".join(tool.get("supports", [])),
        "available": bool(available),
        "binary_path": available or "",
        "status": status,
        "reason": reason,
        "official_url": tool.get("official_url", ""),
        "source_url": tool.get("source_url", ""),
        "notes": tool.get("notes", ""),
    }


def build_plan_rows(manifest_rows, tool_rows):
    plan_rows = []
    for manifest_row in manifest_rows:
        for tool in tool_rows:
            plan_rows.append(
                {
                    "pair_id": manifest_row["pair_id"],
                    "case_id": manifest_row["case_id"],
                    "template_id": manifest_row["template_id"],
                    "tool_id": tool["tool_id"],
                    "tool_name": tool["tool_name"],
                    "integration_kind": tool["integration_kind"],
                    "available": tool["available"],
                    "status": "ready" if tool["available"] else "skipped",
                    "reason": tool["reason"],
                    "official_url": tool["official_url"],
                    "source_url": tool["source_url"],
                    "notes": tool["notes"],
                }
            )
    return plan_rows


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_probe(manifest_path, registry_path, output_root, tools=None):
    manifest_rows = load_manifest(manifest_path)
    registry = load_tool_registry(registry_path)
    if tools:
        registry = [tool for tool in registry if tool["tool_id"] in tools]
    tool_rows = [probe_tool(tool) for tool in registry]
    plan_rows = build_plan_rows(manifest_rows, tool_rows)

    output_root = Path(output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    write_csv(
        output_root / "tool_capability_report.csv",
        tool_rows,
        fieldnames=[
            "tool_id",
            "tool_name",
            "integration_kind",
            "supports",
            "available",
            "binary_path",
            "status",
            "reason",
            "official_url",
            "source_url",
            "notes",
        ],
    )
    write_csv(
        output_root / "tool_baseline_plan.csv",
        plan_rows,
        fieldnames=[
            "pair_id",
            "case_id",
            "template_id",
            "tool_id",
            "tool_name",
            "integration_kind",
            "available",
            "status",
            "reason",
            "official_url",
            "source_url",
            "notes",
        ],
    )
    summary = {
        "manifest": str(Path(manifest_path).resolve()),
        "registry": str(Path(registry_path).resolve()),
        "tool_count": len(tool_rows),
        "available_tools": [row["tool_id"] for row in tool_rows if row["available"]],
        "skipped_tools": [row["tool_id"] for row in tool_rows if not row["available"]],
        "pair_count": len(manifest_rows),
        "plan_rows": len(plan_rows),
        "output_root": str(output_root),
    }
    with open(output_root / "run_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main():
    parser = argparse.ArgumentParser(description="Probe state-of-the-art protein-DNA baseline tools for Dockground benchmarking")
    parser.add_argument("--manifest", default="benchmark/data/protein_dna_dockground_subset_manifest.csv")
    parser.add_argument("--registry", default="benchmark/data/protein_dna_external_tools.json")
    parser.add_argument("--output-root", default="benchmark/prism_processed/results/protein_dna_external_baselines")
    parser.add_argument(
        "--tools",
        default="lightdock,haddock3,hdock,pydockdna",
        help="Comma-separated tool ids to include",
    )
    args = parser.parse_args()

    selected_tools = [tool.strip() for tool in args.tools.split(",") if tool.strip()]
    summary = run_probe(args.manifest, args.registry, args.output_root, selected_tools)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
