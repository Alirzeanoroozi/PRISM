#!/usr/bin/env python3
"""Validate the repository runtime contract without running scientific tools."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(root: Path = ROOT) -> list[str]:
    errors: list[str] = []
    environment = root / "environment.yaml"
    if not environment.is_file():
        errors.append("environment.yaml is missing")
    manifest_path = root / "runtime_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"runtime_manifest.json is invalid: {exc}"]
    required = ("schema_version", "python_scoring_environment", "external_tools", "hpc")
    errors.extend(f"runtime manifest missing {key}" for key in required if key not in manifest)
    tools = manifest.get("external_tools", {})
    tmalign = tools.get("tmalign", {})
    expected_tmalign = root / "external_tools/TMalign"
    if expected_tmalign.is_file() and tmalign.get("sha256") != sha256(expected_tmalign):
        errors.append("TMalign hash does not match runtime_manifest.json")
    fiberdock = tools.get("fiberdock", {})
    expected_fiber = root / "working_version/multiprot/external_tools/fiberdock/FiberDock"
    if expected_fiber.is_file() and fiberdock.get("fiberdock_sha256") != sha256(expected_fiber):
        errors.append("FiberDock hash does not match runtime_manifest.json")
    if fiberdock.get("status") != "energy-only-confirmed; full-refinement-blocked-pending-reduce2-runtime":
        errors.append("FiberDock status must remain fail-closed until full refinement is proven")
    return errors


def main() -> int:
    errors = validate()
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("runtime manifest valid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
