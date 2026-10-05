#!/usr/bin/env python3
"""Create the frozen provenance and template preflight artifacts for a run.

This command is intentionally non-destructive.  It writes only to the output
directory supplied by the caller and can be used as the gate before alignment,
refinement, or scoring jobs are submitted.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

# Support both ``python -m benchmark.scripts.run_investigation_preflight`` and
# direct execution by path from the repository root.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.scripts.investigation_provenance import (
    build_provenance_manifest,
    preflight_template_assets,
    write_provenance_manifest,
    write_template_preflight_tsv,
)


def _split_values(values: list[str]) -> list[str]:
    result = []
    for value in values:
        result.extend(item.strip() for item in value.split(",") if item.strip())
    return result


def build_run_preflight(args: argparse.Namespace) -> dict:
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    repo_root = Path(args.repo_root).resolve()

    files = _split_values(args.file)
    arm_specs = []
    for spec in args.arm:
        parts = spec.split("=", 2)
        if len(parts) != 3 or not all(parts):
            raise ValueError("--arm must use NAME=MANIFEST=ASSET_ROOT")
        arm_specs.append(tuple(parts))
    manifest_paths = [args.template_manifest] if args.template_manifest else []
    manifest_paths.extend(spec[1] for spec in arm_specs)
    files.extend(manifest_paths)
    executables = _split_values(args.executable)
    packages = _split_values(args.package)
    environment = _split_values(args.env)
    seeds = dict(item.split("=", 1) for item in args.seed if "=" in item)

    manifest = build_provenance_manifest(
        repo_root,
        files=files,
        executables=executables,
        environment=environment or None,
        packages=packages,
        command=args.command or None,
        seeds=seeds or None,
        config=args.config,
    )
    provenance_path = write_provenance_manifest(output_dir / "provenance.json", manifest)

    template_report = None
    arm_reports = {}
    if arm_specs:
        for name, manifest_path, asset_root in arm_specs:
            report = preflight_template_assets(manifest_path, asset_root, template_column=args.template_column)
            arm_reports[name] = report
            write_template_preflight_tsv(output_dir / f"template_assets_{name}.tsv", report)
        (output_dir / "template_preflight_arms.json").write_text(
            json.dumps(arm_reports, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
    elif args.template_manifest:
        if not args.asset_root:
            raise ValueError("--asset-root is required with --template-manifest")
        template_report = preflight_template_assets(
            args.template_manifest,
            args.asset_root,
            template_column=args.template_column,
        )
        write_template_preflight_tsv(output_dir / "template_assets.tsv", template_report)
        (output_dir / "template_preflight.json").write_text(
            json.dumps(template_report, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )

    summary = {
        "provenance": str(provenance_path),
        "template_assets": str(output_dir / "template_assets.tsv") if template_report is not None else None,
        "template_preflight": template_report,
        "template_preflight_arms": arm_reports,
    }
    (output_dir / "preflight_summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--template-manifest")
    parser.add_argument("--template-column")
    parser.add_argument("--asset-root")
    parser.add_argument(
        "--arm",
        action="append",
        default=[],
        help="Named template arm in the form NAME=MANIFEST=ASSET_ROOT; repeat for current/legacy comparison",
    )
    parser.add_argument("--file", action="append", default=[])
    parser.add_argument("--executable", action="append", default=[])
    parser.add_argument("--package", action="append", default=[])
    parser.add_argument("--env", action="append", default=[])
    parser.add_argument("--seed", action="append", default=[])
    parser.add_argument("--config")
    parser.add_argument("--command", action="append", default=[])
    parser.add_argument(
        "--require-all-templates",
        action="store_true",
        help="Return exit code 2 when any listed template is invalid or unresolved",
    )
    args = parser.parse_args(argv)

    try:
        summary = build_run_preflight(args)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))

    reports = []
    if summary.get("template_preflight") is not None:
        reports.append(("default", summary["template_preflight"]))
    reports.extend(summary.get("template_preflight_arms", {}).items())
    for name, report in reports:
        print(
            f"template preflight ({name}): "
            f"listed={report['listed']} unique={report['unique']} "
            f"valid={report['valid']} fully_resolvable={report['fully_resolvable']} "
            f"missing={report['missing']}"
        )
    if args.require_all_templates and any(report["missing"] for _, report in reports):
        print("template preflight failed: unresolved templates remain", file=sys.stderr)
        return 2
    print(f"wrote preflight artifacts to {Path(args.output_dir).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
