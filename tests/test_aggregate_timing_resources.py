from __future__ import annotations

import csv
import math
from pathlib import Path

from benchmark.scripts.aggregate_timing_resources import aggregate


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_aggregate_timing_keeps_refinement_and_resource_fields(tmp_path: Path) -> None:
    timing = tmp_path / "timing.tsv"
    refinement = tmp_path / "refinement.tsv"
    _write(timing, [{"stage": "alignment", "wall_seconds": "10", "cpus_per_task": "4", "status": "completed"}])
    _write(refinement, [{
        "fiberdock_elapsed_seconds": "2", "fiberdock_status": "scored",
        "external_rosetta_elapsed_seconds": "3", "external_rosetta_status": "no_model",
    }])
    result = aggregate(tmp_path / "out", [("tmalign", timing)], [("tmalign", refinement)])
    assert result["row_count"] == 3
    with (tmp_path / "out" / "timing_resources.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    alignment = next(row for row in rows if row["stage"] == "alignment")
    assert math.isclose(float(alignment["cpu_core_hours"]), 10 / 3600 * 4)
    assert {row["stage"] for row in rows} == {
        "alignment", "refinement_fiberdock", "refinement_external_rosetta"
    }
