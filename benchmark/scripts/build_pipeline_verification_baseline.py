#!/usr/bin/env python3
"""Freeze retained PRISM-run evidence into a claim and artifact baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections.abc import Iterable
from pathlib import Path


def _read_exit_status(run_root: Path) -> dict[str, object]:
    path = run_root / "status" / "exit.json"
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _scheduler_cancelled(run_root: Path, scheduler_logs: Iterable[Path] = ()) -> bool:
    paths = [*run_root.glob("slurm-*.err"), *scheduler_logs]
    for path in paths:
        if "CANCELLED" in path.read_text(encoding="utf-8", errors="replace").upper():
            return True
    return False


def classify_retained_run(run_root: Path, scheduler_logs: Iterable[Path] = ()) -> dict[str, str]:
    """Classify retained execution evidence without inferring model quality."""
    exit_status = _read_exit_status(run_root)
    scientific_status = str(exit_status.get("scientific_status", "missing"))
    if _scheduler_cancelled(run_root, scheduler_logs):
        reason = "scheduler_cancelled_status_contradiction"
        classification = "unsupported" if scientific_status == "completed" else "unresolved"
    elif scientific_status == "completed":
        reason = "retained_status_completed"
        classification = "supported"
    else:
        reason = f"retained_status_{scientific_status}"
        classification = "unresolved"
    return {
        "run_root": str(run_root),
        "classification": classification,
        "reason": reason,
        "scientific_status": scientific_status,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({field for row in rows for field in row})
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _artifact_rows(run_root: Path, external_paths: Iterable[Path] = ()) -> list[dict[str, str]]:
    rows = [
        {
            "run_root": str(run_root),
            "relative_path": str(path.relative_to(run_root)),
            "sha256": _sha256(path),
            "bytes": str(path.stat().st_size),
        }
        for path in sorted(run_root.rglob("*"))
        if path.is_file()
    ]
    for path in external_paths:
        if not path.is_file():
            raise FileNotFoundError(f"declared evidence path does not exist: {path}")
        try:
            relative_path = str(path.relative_to(run_root))
        except ValueError:
            relative_path = f"external:{path.name}"
        if any(row["relative_path"] == relative_path for row in rows):
            continue
        rows.append(
            {
                "run_root": str(run_root),
                "relative_path": relative_path,
                "sha256": _sha256(path),
                "bytes": str(path.stat().st_size),
            }
        )
    return rows


def build_baseline(config_path: Path, output_dir: Path) -> None:
    """Write deterministic retained-run claim and artifact manifests."""
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to mix baseline output into non-empty directory: {output_dir}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    claims: list[dict[str, str]] = []
    artifacts: list[dict[str, str]] = []
    for item in config.get("runs", []):
        run_root = Path(item["run_root"])
        if not run_root.is_dir():
            raise FileNotFoundError(f"retained run root does not exist: {run_root}")
        scheduler_logs = [Path(value) for value in item.get("scheduler_logs", [])]
        claim = classify_retained_run(run_root, scheduler_logs)
        claim.update(
            {
                "claim_id": str(item["claim_id"]),
                "claim": str(item.get("claim", "")),
                "scope": str(item.get("scope", "execution")),
                "evidence_paths": str(run_root),
            }
        )
        claims.append(claim)
        artifacts.extend(_artifact_rows(run_root, scheduler_logs))
    _write_tsv(output_dir / "claims.tsv", claims)
    _write_tsv(output_dir / "artifact_manifest.tsv", artifacts)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build_baseline(args.config, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
