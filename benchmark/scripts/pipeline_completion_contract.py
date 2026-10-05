#!/usr/bin/env python3
"""Classify PRISM pipeline completion from explicit stage evidence."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


CURRENT_STAGES = ("input", "alignment", "transformation", "refinement")
TERMINAL_EVENTS = {
    "completed",
    "skipped",
    "failed",
    "timeout",
    "partial",
    "unavailable",
    "unknown",
}


def write_stage_event(path: Path, stage: str, event: str, return_code: int | None = None, detail: str = "") -> None:
    """Append one durable stage event when stage-status recording is enabled."""
    record: dict[str, object] = {
        "stage": stage,
        "event": event,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "detail": detail,
    }
    if return_code is not None:
        record["return_code"] = return_code
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def _completed_stages(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    completed: set[str] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if record.get("event") in TERMINAL_EVENTS:
            completed.add(str(record.get("stage", "")))
    return completed


def _paired_transformations(run_root: Path) -> int:
    directory = run_root / "processed" / "transformation"
    if not directory.is_dir():
        return 0
    left_stems = {
        path.name[: -len("_L.pdb")]
        for path in directory.glob("*_L.pdb")
        if path.stat().st_size > 0
    }
    right_stems = {
        path.name[: -len("_R.pdb")]
        for path in directory.glob("*_R.pdb")
        if path.stat().st_size > 0
    }
    return len(left_stems.intersection(right_stems))


def _refined_models(run_root: Path) -> int:
    external = run_root / "processed" / "rosetta_refinement"
    pyrosetta = run_root / "processed" / "pyrosetta_refinement" / "structures"
    return sum(
        1
        for directory in (external, pyrosetta)
        if directory.is_dir()
        for path in directory.glob("*.pdb")
        if path.stat().st_size > 0
    )


def _valid_refined_models(run_root: Path) -> int:
    """Count nonempty refined PDBs containing at least two monotonic chains."""

    directories = (
        run_root / "processed" / "rosetta_refinement",
        run_root / "processed" / "pyrosetta_refinement" / "structures",
    )
    valid = 0
    for directory in directories:
        if not directory.is_dir():
            continue
        for path in directory.glob("*.pdb"):
            if not path.is_file() or path.stat().st_size == 0:
                continue
            chains: dict[str, list[int]] = {}
            try:
                with path.open(errors="replace") as handle:
                    for line in handle:
                        if not line.startswith("ATOM  ") or len(line) < 27:
                            continue
                        chain = line[21].strip() or "_"
                        try:
                            residue = int(line[22:26])
                        except ValueError:
                            continue
                        chains.setdefault(chain, []).append(residue)
            except OSError:
                continue
            if len(chains) < 2 or any(
                any(later < earlier for earlier, later in zip(numbers, numbers[1:]))
                for numbers in chains.values()
            ):
                continue
            valid += 1
    return valid


def classify_run(run_root: Path, process_return_code: int, termination_signal: str | None) -> dict[str, object]:
    """Return a fail-closed scheduler and scientific status for a current run."""
    paired_count = _paired_transformations(run_root)
    refined_count = _refined_models(run_root)
    valid_refined_count = _valid_refined_models(run_root)
    base: dict[str, object] = {
        "paired_transformation_count": paired_count,
        "refined_model_count": refined_count,
        "valid_refined_model_count": valid_refined_count,
        "process_return_code": process_return_code,
        "termination_signal": termination_signal or "",
    }
    if termination_signal:
        return base | {
            "process_status": "cancelled",
            "scientific_status": "cancelled",
            "reason": f"termination_signal:{termination_signal}",
        }
    if process_return_code != 0:
        return base | {
            "process_status": "failed",
            "scientific_status": "failed",
            "reason": f"return_code:{process_return_code}",
        }
    status_dir = run_root / "status"
    if not (status_dir / "pipeline_returned.json").is_file():
        return base | {
            "process_status": "completed",
            "scientific_status": "incomplete",
            "reason": "missing_pipeline_returned",
        }
    missing = [stage for stage in CURRENT_STAGES if stage not in _completed_stages(status_dir / "stages.jsonl")]
    if missing:
        return base | {
            "process_status": "completed",
            "scientific_status": "incomplete",
            "reason": f"missing_terminal_stages:{','.join(missing)}",
        }
    if refined_count == 0:
        return base | {
            "process_status": "completed",
            "scientific_status": "completed_no_predictions",
            "reason": "no_refined_models",
        }
    if valid_refined_count == 0:
        return base | {
            "process_status": "completed",
            "scientific_status": "incomplete",
            "reason": "no_valid_refined_models",
        }
    return base | {
        "process_status": "completed",
        "scientific_status": "completed",
        "reason": "terminal_stages_and_refined_models_present",
    }


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--return-code", type=int, required=True)
    parser.add_argument("--termination-signal", default="")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pipeline", required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--elapsed-seconds", type=int, required=True)
    parser.add_argument("--slurm-job-id", default="")
    parser.add_argument("--array-task-id", default="")
    args = parser.parse_args(argv)
    result = classify_run(args.run_root, args.return_code, args.termination_signal or None)
    result.update(
        {
            "return_code": args.return_code,
            "elapsed_seconds": args.elapsed_seconds,
            "pipeline": args.pipeline,
            "aligner": args.aligner,
            "slurm_job_id": args.slurm_job_id,
            "array_task_id": args.array_task_id,
        }
    )
    _write_json(args.output, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
