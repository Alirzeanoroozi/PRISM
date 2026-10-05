#!/usr/bin/env python3
"""Summarize resumable FiberDock/Rosetta/DockQ refinement checkpoints."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    files = sorted((root / "checkpoints").glob("*.json"))
    record_status = Counter()
    stage_status: dict[str, Counter[str]] = defaultdict(Counter)
    pipeline_stage: dict[str, Counter[str]] = defaultdict(Counter)
    scored_values: dict[tuple[str, str], list[float]] = defaultdict(list)
    both_scored: dict[str, set[str]] = defaultdict(set)
    seen_scored: dict[str, set[str]] = defaultdict(set)
    parse_errors = Counter()

    for path in files:
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:  # pragma: no cover - diagnostic robustness
            parse_errors[type(exc).__name__] += 1
            continue
        key = str(record.get("key", path.stem))
        pipeline = next((p for p in ("multiprot", "tmalign", "usalign") if p in key), "other")
        record_status[str(record.get("status"))] += 1
        stages = record.get("stages") or {}
        for name, value in stages.items():
            if not isinstance(value, dict):
                continue
            status = str(value.get("status"))
            stage_status[name][status] += 1
            pipeline_stage[pipeline][f"{name}:{status}"] += 1
            if name.startswith("dockq_") and status == "scored":
                seen_scored[pipeline].add(key)
                score = value.get("dockq_global", value.get("dockq"))
                if isinstance(score, (int, float)):
                    scored_values[(pipeline, name)].append(float(score))
            if name in {"dockq_fiberdock", "dockq_rosetta"} and status == "scored":
                other = "dockq_rosetta" if name == "dockq_fiberdock" else "dockq_fiberdock"
                if isinstance(stages.get(other), dict) and stages[other].get("status") == "scored":
                    both_scored[pipeline].add(key)

    print(f"checkpoint_files={len(files)}")
    print(f"record_status={dict(sorted(record_status.items()))}")
    print(f"stage_status={{{', '.join(f'{k!r}: {dict(sorted(v.items()))!r}' for k, v in sorted(stage_status.items()))}}}")
    print(f"pipeline_stage={{{', '.join(f'{k!r}: {dict(sorted(v.items()))!r}' for k, v in sorted(pipeline_stage.items()))}}}")
    for (pipeline, stage), values in sorted(scored_values.items()):
        print(
            f"dockq_stats pipeline={pipeline} stage={stage} n={len(values)} "
            f"mean={sum(values) / len(values):.6f} min={min(values):.6f} max={max(values):.6f}"
        )
    for pipeline in sorted(seen_scored):
        print(f"dockq_any_scored pipeline={pipeline} n={len(seen_scored[pipeline])}")
        print(f"dockq_both_scored pipeline={pipeline} n={len(both_scored[pipeline])}")
    print(f"parse_errors={dict(parse_errors)}")
    results = [p for p in sorted((root / "results").glob("**/*")) if p.is_file()]
    print(f"result_files={len(results)}")
    for path in results:
        print(f"result_file={path.relative_to(root)} size={path.stat().st_size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
