#!/usr/bin/env python3
"""Score compacted USalign transformation candidates with the common evaluator.

This stage consumes only the common transformation outputs and candidate
ledger.  It never reruns USalign.  Each candidate gets an explicit
GlobalDockQ, requested cross-interface DockQ components, mapping, iRMSD
status, evaluator provenance, and a resumable checkpoint before disposable
combined/raw scoring files are removed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Bio.PDB import PDBParser

from benchmark.scripts.score_bijective_benchmark_models import (
    assign_group,
    dockq_command,
    dockq_version,
    is_recoverable_complete_mapping_error,
    parse_complex,
    run_dockq,
    run_irmsd,
    run_pairwise_cross_dockq,
    sha256_file,
    standardize_dockq_json,
)
from benchmark.scripts.standardized_evaluator import validate_raw_pdb_chain_contract


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def candidate_pair_paths(run_root: Path, row: dict[str, str]) -> tuple[Path, Path]:
    prefix = "{}_{}_{}_{}".format(
        row["template"], row["query_left"], row["query_right"], row["orientation"]
    )
    root = run_root / "processed" / "transformation"
    return root / f"{prefix}_L.pdb", root / f"{prefix}_R.pdb"


def source_chain_ids(path: Path) -> str:
    """Recover source partner chains from the common transformation filename."""

    fields = path.stem.split("_")
    if len(fields) >= 5:
        target = fields[1] if fields[-1] == "L" else fields[2]
        if len(target) >= 5:
            return "".join(dict.fromkeys(target[4:]))
    chains = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM") and len(line) >= 22:
                chain = line[21].strip() or "_"
                if chain not in chains:
                    chains.append(chain)
    if not chains:
        raise ValueError(f"no source chains found in {path}")
    return "".join(chains)


def partner_chain_groups(left: str, right: str) -> tuple[str, str]:
    """Match the common Rosetta chain-renaming contract without side effects."""

    if not set(left).intersection(right):
        return left, right
    available = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    used: set[str] = set()
    new_left = ""
    for _ in left:
        candidate = next(candidate for candidate in available if candidate not in used)
        used.add(candidate)
        new_left += candidate
    new_right = ""
    for _ in right:
        candidate = next(candidate for candidate in available if candidate not in used)
        used.add(candidate)
        new_right += candidate
    return new_left, new_right


def combine_transformed_pair(left: Path, right: Path, output: Path) -> tuple[str, str]:
    """Assemble transformation halves using the same chain contract as Rosetta."""

    left_source = source_chain_ids(left)
    right_source = source_chain_ids(right)
    left_group, right_group = partner_chain_groups(left_source, right_source)
    output.parent.mkdir(parents=True, exist_ok=True)
    with left.open(encoding="ascii", errors="replace") as left_handle, right.open(
        encoding="ascii", errors="replace"
    ) as right_handle, output.open("w", encoding="ascii") as combined:
        for handle, source, target in (
            (left_handle, left_source, left_group),
            (right_handle, right_source, right_group),
        ):
            for line in handle:
                if not line.startswith("ATOM") or len(line) < 22:
                    continue
                chain = line[21].strip()
                if chain not in source:
                    continue
                replacement = target[source.index(chain)]
                combined.write(line[:21] + replacement + line[22:])
            combined.write("TER\n")
        combined.write("END\n")
    return left_group, right_group


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _input_index(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    return {
        (row.get("Receptor", "").strip().lower(), row.get("Ligand", "").strip().lower()): row
        for row in rows
    }


def iter_bounded_results(tasks, worker, workers: int):
    """Yield scored candidates with at most 2*workers futures in flight."""

    iterator = iter(tasks)
    pending_limit = max(1, 2 * workers)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        pending = {}
        exhausted = False
        while pending or not exhausted:
            while not exhausted and len(pending) < pending_limit:
                try:
                    task = next(iterator)
                except StopIteration:
                    exhausted = True
                    break
                pending[executor.submit(worker, *task)] = task
            if not pending:
                break
            future = next(iter(as_completed(tuple(pending))))
            pending.pop(future)
            yield future.result()


def _base_result(index: int, row: dict[str, str]) -> dict[str, Any]:
    result: dict[str, Any] = dict(row)
    result.update(
        {
            "candidate_index": index,
            "score_status": "not_attempted",
            "score_scope": "",
            "score_error": "",
            "dockq_global": "",
            "dockq_global_status": "not_attempted",
            "dockq_cross_best": "",
            "dockq_cross_mean": "",
            "dockq_cross_components": "",
            "dockq_cross_requested_count": "",
            "dockq_cross_scoreable_count": "",
            "dockq_cross_unscoreable_count": "",
            "dockq_cross_unscoreable_interfaces": "",
            "irmsd_status": "not_attempted",
            "irmsd_grouped_forward": "",
            "irmsd_grouped_reverse": "",
            "irmsd_grouped_min": "",
            "dockq_mapping": "",
            "mapping_validation_status": "",
            "receptor_assignment": "",
            "ligand_assignment": "",
            "dockq_version": "",
            "score_elapsed_seconds": "",
            "model_sha256": "",
            "native_sha256": "",
            "raw_dockq_json_sha256": "",
            "transformed_left": "",
            "transformed_right": "",
            "combined_model": "",
        }
    )
    return result


def score_candidate(
    index: int,
    row: dict[str, str],
    inputs: dict[tuple[str, str], dict[str, str]],
    run_root: Path,
    native_root: Path,
    score_python: Path,
    score_root: Path,
    timeout: int,
    installed_dockq_version: str,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    checkpoint = checkpoint_dir / f"candidate_{index:08d}.json"
    if checkpoint.is_file():
        return json.loads(checkpoint.read_text())

    started = time.perf_counter()
    result = _base_result(index, row)

    def persist() -> dict[str, Any]:
        result["score_elapsed_seconds"] = round(time.perf_counter() - started, 6)
        _write_json(checkpoint, result)
        return result

    left, right = candidate_pair_paths(run_root, row)
    result["transformed_left"] = str(left)
    result["transformed_right"] = str(right)
    if row.get("status") != "generated":
        result.update(score_status="not_scoreable", score_error=row.get("error_reason") or row.get("status", "not_generated"))
        return persist()
    if not left.is_file() or not right.is_file():
        result.update(score_status="not_scoreable", score_error="transformed_pair_missing")
        return persist()

    input_row = inputs.get((row.get("query_left", "").strip().lower(), row.get("query_right", "").strip().lower()))
    if input_row is None:
        result.update(score_status="not_scoreable", score_error="benchmark_input_pair_missing")
        return persist()

    benchmark_set = input_row["benchmark_set"]
    pdb_id, native_receptor, native_ligand = parse_complex(input_row["complex"])
    native = native_root / f"native_bound_complexes_t_{benchmark_set}" / f"{pdb_id}.pdb"
    result.update(
        {
            "pair_id": input_row.get("pair_id", ""),
            "benchmark_set": benchmark_set,
            "complex": input_row["complex"],
            "native_pdb": str(native),
            "native_receptor_chains": native_receptor,
            "native_ligand_chains": native_ligand,
            "native_sha256": sha256_file(native) if native.is_file() else "",
        }
    )
    if not native.is_file():
        result.update(score_status="not_scoreable", score_error="native_pdb_missing")
        return persist()

    combined = score_root / "combined" / f"candidate_{index:08d}.pdb"
    raw_json = score_root / "dockq_json" / f"candidate_{index:08d}.json"
    result["combined_model"] = str(combined)
    try:
        model_receptor_hint, model_ligand_hint = combine_transformed_pair(left, right, combined)
        result["model_sha256"] = sha256_file(combined)
        contract = validate_raw_pdb_chain_contract(combined, model_receptor_hint, model_ligand_hint)
        if not contract.valid:
            raise ValueError("combined model chain contract rejected: " + "; ".join(contract.errors))
        parser = PDBParser(QUIET=True)
        model_structure = parser.get_structure("model", str(combined))
        native_structure = parser.get_structure("native", str(native))
        model_receptor, receptor_diag = assign_group(
            model_structure, native_structure, model_receptor_hint, native_receptor
        )
        model_ligand, ligand_diag = assign_group(
            model_structure, native_structure, model_ligand_hint, native_ligand
        )
        mapping = f"{model_receptor}{model_ligand}:{native_receptor}{native_ligand}"
        result.update(
            {
                "dockq_mapping": mapping,
                "mapping_validation_status": "aligned_default",
                "receptor_assignment": receptor_diag,
                "ligand_assignment": ligand_diag,
                "dockq_version": installed_dockq_version,
            }
        )
        score_scope = "complete_bijective_mapping"
        full_error = ""
        fallback_components: list[dict[str, Any]] = []
        fallback_commands: list[list[str]] = []
        try:
            dockq = run_dockq(score_python, None, combined, native, mapping, raw_json, timeout, False)
            standardized = standardize_dockq_json(raw_json)
            raw_hash = standardized.raw_json_sha256
        except RuntimeError as exc:
            if not is_recoverable_complete_mapping_error(exc):
                raise
            full_error = str(exc)
            dockq, fallback_components, fallback_commands = run_pairwise_cross_dockq(
                score_python, None, combined, native,
                model_receptor, model_ligand, native_receptor, native_ligand,
                raw_json, timeout, False,
            )
            standardized = None
            raw_hash = sha256_file(raw_json)
            score_scope = "requested_cross_interfaces_only"

        cross_keys = [f"{r}{l}" for r in native_receptor for l in native_ligand]
        best_result = dockq.get("best_result", {})
        cross = [best_result[key] for key in cross_keys if key in best_result]
        if not cross:
            raise RuntimeError(f"no requested cross interface in DockQ result: expected={cross_keys}")
        fallback_status = {
            str(item.get("interface", "")): str(item.get("status", ""))
            for item in fallback_components
        }
        unscoreable = [key for key in cross_keys if fallback_status.get(key) == "no_native_interface"]
        grouped = reciprocal = None
        irmsd_status = "scored"
        irmsd_error = ""
        try:
            grouped = run_irmsd(score_python, REPO_ROOT / "benchmark/scripts/irmsd.py", combined, model_receptor, model_ligand, native, native_receptor, native_ligand, timeout)
            reciprocal = run_irmsd(score_python, REPO_ROOT / "benchmark/scripts/irmsd.py", combined, model_ligand, model_receptor, native, native_ligand, native_receptor, timeout)
        except Exception as exc:
            irmsd_status = "failed_auxiliary"
            irmsd_error = str(exc)
        result.update(
            {
                "score_status": (
                    "scored_cross_only"
                    if score_scope == "requested_cross_interfaces_only"
                    else "scored"
                    if dockq.get("GlobalDockQ") is not None
                    else "valid_unscored"
                ),
                "score_scope": score_scope,
                "score_error": full_error,
                "dockq_global": "" if score_scope != "complete_bijective_mapping" else str(dockq.get("GlobalDockQ", "")),
                "dockq_global_status": (
                    "unavailable_cross_only"
                    if score_scope != "complete_bijective_mapping"
                    else "scored"
                    if dockq.get("GlobalDockQ") is not None
                    else "valid_unscored"
                ),
                "dockq_cross_best": str(max(float(item["DockQ"]) for item in cross)),
                "dockq_cross_mean": str(sum(float(item["DockQ"]) for item in cross) / len(cross)),
                "dockq_cross_components": json.dumps(cross, sort_keys=True),
                "dockq_cross_requested_count": str(len(cross_keys)),
                "dockq_cross_scoreable_count": str(len(cross)),
                "dockq_cross_unscoreable_count": str(len(unscoreable)),
                "dockq_cross_unscoreable_interfaces": json.dumps(unscoreable),
                "irmsd_status": irmsd_status,
                "irmsd_error": irmsd_error,
                "irmsd_grouped_forward": "" if grouped is None else str(grouped),
                "irmsd_grouped_reverse": "" if reciprocal is None else str(reciprocal),
                "irmsd_grouped_min": "" if grouped is None or reciprocal is None else str(min(grouped, reciprocal)),
                "raw_dockq_json_sha256": raw_hash,
                "dockq_best_internal": "" if score_scope != "complete_bijective_mapping" else str(dockq.get("best_dockq", "")),
                "dockq_argv": json.dumps(dockq_command(score_python, None, combined, native, mapping, raw_json, False)),
            }
        )
        if standardized is not None:
            result["dockq_interface_records"] = [dict(item) for item in standardized]
        else:
            result["dockq_interface_records"] = [
                {"interface": key, **best_result[key]} for key in sorted(best_result)
            ]
    except Exception as exc:
        result.update(
            score_status="score_failed",
            score_scope="",
            score_error=str(exc),
            dockq_version=installed_dockq_version,
        )

    # The checkpoint is the compact scientific record.  It is written before
    # the combined model/raw JSON cleanup so an interrupted task is resumable.
    persist()
    for disposable in (combined, raw_json):
        if disposable.is_file():
            disposable.unlink()
    result["disposable_scoring_files_deleted"] = True
    persist()
    return result


def write_outputs(results: list[dict[str, Any]], output: Path, interfaces_output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in results for key in row if key != "dockq_interface_records"})
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    interface_rows = []
    for row in results:
        records = row.get("dockq_interface_records", [])
        if isinstance(records, list):
            for record in records:
                interface_rows.append(
                    {
                        "candidate_index": row.get("candidate_index", ""),
                        "pair_id": row.get("pair_id", ""),
                        "template": row.get("template", ""),
                        "score_status": row.get("score_status", ""),
                        "requested_cross_interface": str(record.get("interface", "")) in {
                            f"{r}{l}" for r in str(row.get("native_receptor_chains", "")) for l in str(row.get("native_ligand_chains", ""))
                        },
                        **record,
                    }
                )
    interfaces_output.parent.mkdir(parents=True, exist_ok=True)
    interface_fields = sorted({key for row in interface_rows for key in row}) or ["candidate_index"]
    with interfaces_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=interface_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(interface_rows)


def cleanup_transformed_halves(
    results: list[dict[str, Any]],
    output: Path,
    *,
    preserve_for_downstream: bool = False,
) -> dict[str, Any]:
    paths = []
    for row in results:
        # A transformed pair is disposable only after its compact score record
        # is valid.  Keep score failures and unresolved non-scoreable rows for
        # audit/debugging instead of making the only structural evidence
        # unrecoverable.
        if row.get("status") == "generated" and row.get("score_status") in {
            "scored", "scored_cross_only", "valid_unscored",
        }:
            paths.extend([row.get("transformed_left", ""), row.get("transformed_right", "")])
    deleted = []
    retained = []
    if preserve_for_downstream:
        retained = sorted(set(paths))
    else:
        for raw_path in sorted(set(paths)):
            if not raw_path:
                continue
            path = Path(raw_path)
            if path.is_file():
                path.unlink()
                deleted.append(str(path))
            else:
                retained.append(str(path))
    manifest = {
        "schema_version": "prism-usalign-transformed-dockq-cleanup/v1",
        "candidate_count": len(results),
        "planned_transformed_half_count": len(set(paths)),
        "deleted_transformed_half_count": len(deleted),
        "missing_or_already_deleted_count": len(retained),
        "preserved_for_downstream": preserve_for_downstream,
        "retention_reason": (
            "common_refinement_pending"
            if preserve_for_downstream
            else "validated_score_and_no_structural_consumer"
        ),
        "retained_unresolved_candidate_count": sum(
            row.get("status") == "generated" and row.get("score_status") not in {
                "scored", "scored_cross_only", "valid_unscored",
            }
            for row in results
        ),
        "deleted_paths_sha256": hashlib.sha256("\n".join(deleted).encode()).hexdigest(),
    }
    _write_json(output, manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--score-python", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--retain-transformed-inputs",
        action="store_true",
        help="Retain transformed halves for a downstream refinement consumer.",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")

    run_root = args.run_root.resolve()
    status_dir = run_root / "status"
    candidates_path = status_dir / "candidate_generated.csv"
    inputs_path = run_root / "inputs.csv"
    output = status_dir / "transformed_dockq.tsv"
    interfaces_output = status_dir / "transformed_dockq_interfaces.tsv"
    cleanup_output = status_dir / "transformed_dockq_cleanup.json"
    if not candidates_path.is_file() or not inputs_path.is_file():
        _write_json(status_dir / "transformed_dockq_status.json", {"status": "skipped", "reason": "missing_compacted_inputs"})
        return 0

    stage_started = time.perf_counter()
    candidates = read_csv(candidates_path)
    inputs = _input_index(read_csv(inputs_path))
    score_root = status_dir / "transformed_dockq_work"
    checkpoint_dir = score_root / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    score_python = args.score_python.resolve()
    installed_version = dockq_version(score_python)
    tasks = (
        (
            index, row, inputs, run_root, args.native_root.resolve(), score_python,
            score_root, args.timeout, installed_version, checkpoint_dir,
        )
        for index, row in enumerate(candidates)
    )
    results = list(iter_bounded_results(tasks, score_candidate, args.workers))
    results.sort(key=lambda row: int(row["candidate_index"]))
    write_outputs(results, output, interfaces_output)
    cleanup = cleanup_transformed_halves(
        results,
        cleanup_output,
        preserve_for_downstream=args.retain_transformed_inputs,
    )
    summary = {
        "schema_version": "prism-usalign-transformed-dockq/v1",
        "status": "validated_compacted",
        "candidate_rows": len(results),
        "generated_rows": sum(row.get("status") == "generated" for row in results),
        "scored_rows": sum(row.get("score_status") in {"scored", "scored_cross_only"} for row in results),
        "valid_unscored_rows": sum(row.get("score_status") == "valid_unscored" for row in results),
        "score_failed_rows": sum(row.get("score_status") == "score_failed" for row in results),
        "not_scoreable_rows": sum(row.get("score_status") == "not_scoreable" for row in results),
        "dockq_version": installed_version,
        "score_python": str(score_python),
        "workers": args.workers,
        "timeout": args.timeout,
        "wall_seconds": round(time.perf_counter() - stage_started, 6),
        "output": str(output),
        "interfaces_output": str(interfaces_output),
        "cleanup": cleanup,
    }
    _write_json(status_dir / "transformed_dockq_status.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
