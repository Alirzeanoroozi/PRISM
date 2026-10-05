#!/usr/bin/env python3
"""Run isolated cancellation/no-prediction/positive-pose completion canaries."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

try:
    from benchmark.scripts.pipeline_completion_contract import classify_run
except ModuleNotFoundError:  # Direct Slurm/script invocation.
    from pipeline_completion_contract import classify_run


def _events(status: Path) -> None:
    status.mkdir(parents=True, exist_ok=True)
    with (status / "stages.jsonl").open("w", encoding="utf-8") as handle:
        for stage in ("input", "alignment", "transformation", "refinement"):
            handle.write(json.dumps({"stage": stage, "event": "completed", "return_code": 0}) + "\n")
    (status / "pipeline_returned.json").write_text("{}\n", encoding="utf-8")


def _valid_pdb(path: Path) -> None:
    lines = []
    for serial, chain in ((1, "A"), (2, "B")):
        lines.append(f"ATOM  {serial:5d}  CA  ALA {chain}   1      0.000   0.000   0.000  1.00 20.00           C\n")
    path.write_text("".join(lines), encoding="ascii")


def run(output: Path) -> dict[str, object]:
    if output.exists():
        raise RuntimeError(f"refusing nonempty canary output: {output}")
    output.mkdir(parents=True)
    cancelled = output / "cancelled"
    _events(cancelled / "status")
    cancellation = classify_run(cancelled, 0, "SIGTERM")

    no_prediction = output / "no_prediction"
    _events(no_prediction / "status")
    no_prediction_result = classify_run(no_prediction, 0, None)

    positive = output / "positive_full_pose"
    _events(positive / "status")
    refined = positive / "processed" / "rosetta_refinement"
    refined.mkdir(parents=True)
    _valid_pdb(refined / "complete_rosetta_0001.pdb")
    positive_result = classify_run(positive, 0, None)

    results = {"cancelled": cancellation, "no_prediction": no_prediction_result, "positive_full_pose": positive_result}
    expected = {
        "cancelled": "cancelled",
        "no_prediction": "completed_no_predictions",
        "positive_full_pose": "completed",
    }
    for name, scientific_status in expected.items():
        if results[name]["scientific_status"] != scientific_status:
            raise RuntimeError(f"canary {name} mismatch: {results[name]}")
    (output / "canary_results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), sort_keys=True))
