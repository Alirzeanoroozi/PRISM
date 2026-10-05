#!/usr/bin/env python3
"""Compact GTalign DockQ JSON into global and requested-interface tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def number(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def compact(scores_csv: Path, dataset_manifest: Path, output_dir: Path) -> dict[str, Any]:
    manifest = json.loads(dataset_manifest.read_text(encoding="utf-8"))
    dataset = {row["case_id"]: row for row in manifest["pairs"]}
    score_rows = list(csv.DictReader(scores_csv.open(newline="", encoding="utf-8")))
    compact_rows: list[dict[str, Any]] = []
    interface_rows: list[dict[str, Any]] = []
    parse_failures = 0
    for source in score_rows:
        case = dataset.get(source.get("case_id", ""), {})
        row: dict[str, Any] = {
            "case_id": source.get("case_id", ""),
            "split": case.get("split", source.get("benchmark_split", "")),
            "benchmark_complex": case.get("benchmark_complex", ""),
            "native_receptor_chains": case.get("native_receptor_chains", ""),
            "native_ligand_chains": case.get("native_ligand_chains", ""),
            "template": source.get("template", ""),
            "orientation": (
                source.get("orientation", "")
                if str(source.get("orientation", "")).startswith("o")
                else f"o{source.get('orientation', '')}"
            ),
            "receptor": source.get("receptor", ""),
            "ligand": source.get("ligand", ""),
            "query_left": source.get("receptor", ""),
            "query_right": source.get("ligand", ""),
            "chain_left": source.get("model_receptor_chains", ""),
            "chain_right": source.get("model_ligand_chains", ""),
            "stage": source.get("stage", "transformed"),
            "source_status": source.get("status", ""),
            "score_status": "score_failed" if source.get("status") == "score_failed" else "valid_unscored",
            "score_error": source.get("error", ""),
            "dockq_global": "",
            "dockq_best_internal_diagnostic": "",
            "dockq_cross_best": "",
            "dockq_cross_mean": "",
            "dockq_cross_fnat_mean": "",
            "dockq_cross_irmsd_mean": "",
            "dockq_cross_lrmsd_mean": "",
            "dockq_cross_interface_count": 0,
            "dockq_mapping": source.get("dockq_mapping", ""),
            "raw_dockq_json": source.get("raw_dockq_json", ""),
            "raw_dockq_json_sha256": source.get("raw_dockq_json_sha256", ""),
        }
        raw_value = source.get("raw_dockq_json", "")
        raw = Path(raw_value) if raw_value else None
        if raw is None or not raw.is_file():
            compact_rows.append(row)
            continue
        try:
            payload = json.loads(raw.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            parse_failures += 1
            row["score_status"] = "score_failed"
            row["score_error"] = f"invalid_raw_dockq_json:{exc}"
            compact_rows.append(row)
            continue
        row["raw_dockq_json_sha256"] = sha256(raw)
        global_dockq = number(payload.get("GlobalDockQ"))
        row["dockq_global"] = "" if global_dockq is None else global_dockq
        row["dockq_best_internal_diagnostic"] = number(payload.get("best_dockq")) or ""
        # Preserve a terminal source failure even when a partial/raw JSON
        # happened to be emitted.  A parsed payload is evidence for diagnosis,
        # not permission to promote the execution status to scientific
        # success.
        row["score_status"] = (
            "score_failed"
            if source.get("status") == "score_failed"
            else "scored"
            if global_dockq is not None
            else "valid_unscored"
        )
        row["dockq_mapping"] = payload.get("best_mapping_str", row["dockq_mapping"])
        best_result = payload.get("best_result") or {}
        cross_keys = [
            f"{receptor}{ligand}"
            for receptor in case.get("native_receptor_chains", "")
            for ligand in case.get("native_ligand_chains", "")
        ]
        selected = []
        for interface in cross_keys:
            component = best_result.get(interface)
            if not isinstance(component, dict):
                continue
            selected.append(component)
            interface_rows.append({
                "case_id": row["case_id"],
                "template": row["template"],
                "orientation": row["orientation"],
                "interface": interface,
                "requested_cross_interface": True,
                "GlobalDockQ": row["dockq_global"],
                "DockQ": component.get("DockQ", ""),
                "Fnat": component.get("fnat", component.get("Fnat", "")),
                "iRMSD": component.get("iRMSD", ""),
                "LRMSD": component.get("LRMSD", ""),
                "mapping": row["dockq_mapping"],
                "raw_dockq_json": str(raw),
                "raw_dockq_json_sha256": row["raw_dockq_json_sha256"],
            })
        def values(key: str) -> list[float]:
            return [value for item in selected if (value := number(item.get(key))) is not None]
        dockq_values = values("DockQ")
        fnat_values = values("fnat")
        irmsd_values = values("iRMSD")
        lrmsd_values = values("LRMSD")
        row["dockq_cross_best"] = max(dockq_values) if dockq_values else ""
        row["dockq_cross_mean"] = mean(dockq_values) if dockq_values else ""
        row["dockq_cross_fnat_mean"] = mean(fnat_values) if fnat_values else ""
        row["dockq_cross_irmsd_mean"] = mean(irmsd_values) if irmsd_values else ""
        row["dockq_cross_lrmsd_mean"] = mean(lrmsd_values) if lrmsd_values else ""
        row["dockq_cross_interface_count"] = len(selected)
        compact_rows.append(row)

    output_dir.mkdir(parents=True, exist_ok=True)

    def write(name: str, rows: list[dict[str, Any]]) -> str:
        path = output_dir / name
        fields = sorted({key for row in rows for key in row}) or ["status"]
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        return sha256(path)

    rows_hash = write("gtalign_transformed_dockq.tsv", compact_rows)
    interfaces_hash = write("gtalign_transformed_interfaces.tsv", interface_rows)
    status_counts: dict[str, int] = {}
    for row in compact_rows:
        status = row["score_status"]
        status_counts[status] = status_counts.get(status, 0) + 1
    summary = {
        "schema_version": "prism-gtalign-dockq-compact/v1",
        "input_rows": len(score_rows),
        "output_rows": len(compact_rows),
        "interface_rows": len(interface_rows),
        "parse_failures": parse_failures,
        "status_counts": dict(sorted(status_counts.items())),
        "rows_sha256": rows_hash,
        "interfaces_sha256": interfaces_hash,
        "source_scores_csv": str(scores_csv.resolve()),
        "source_dataset_manifest": str(dataset_manifest.resolve()),
    }
    (output_dir / "compact_status.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores-csv", type=Path, required=True)
    parser.add_argument("--dataset-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compact(args.scores_csv, args.dataset_manifest, args.output_dir), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
