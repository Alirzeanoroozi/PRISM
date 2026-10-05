#!/usr/bin/env python3
"""Build a labeled-free candidate table from PRISM alignment artifacts.

This is intentionally a data-extraction tool, not a ranking decision. Missing
alignment or transformation files become explicit status rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.candidate_audit import alignment_features
from benchmark.scripts.stage_current_models_for_main_benchmark import select_final_external_models


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _alignment_path(processed_root: Path, filename: str) -> tuple[Path | None, str]:
    """Resolve one alignment record without silently choosing a run."""

    candidates = []
    direct = processed_root / "alignment" / filename
    if direct.is_file():
        candidates.append(direct)
    candidates.extend(sorted((processed_root / "alignment_gtalign").glob(f"*/{filename}")))
    unique = list(dict.fromkeys(path.resolve() for path in candidates))
    if not unique:
        return None, "missing_alignment_json"
    if len(unique) != 1:
        return None, "ambiguous_alignment_json"
    return unique[0], ""


def _payload_number(payload: dict, field: str):
    value = payload.get(field, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return ""


def _batch_processed_root(model: Path) -> Path:
    for parent in model.parents:
        if parent.name.startswith("batch_"):
            processed = parent / "processed"
            if processed.is_dir():
                return processed
    raise ValueError(f"model is not beneath a batch processed root: {model}")


def build_rows_from_stage_manifest(stage_manifest: Path) -> list[dict]:
    """Build ranking rows from canonical staged models and alignment JSON."""

    with stage_manifest.open(newline="", encoding="utf-8") as handle:
        stages = list(csv.DictReader(handle))
    rows: list[dict] = []
    alignment_cache: dict[Path, dict] = {}
    template_cache: dict[Path, dict] = {}
    for stage in stages:
        if stage.get("status") != "staged_symlink" or stage.get("source_gate_status") == "audit_only":
            continue
        model = Path(stage["source_model_path"])
        template_1 = stage.get("template_1", "")
        template_2 = stage.get("template_2", "")
        if len(template_1) < 5 or len(template_2) < 5 or template_1[:4].lower() != template_2[:4].lower():
            rows.append({**stage, "status": "alignment_failed", "error_reason": "invalid_template_partner_tokens"})
            continue
        template = template_1[:4] + template_1[4:].replace("_", "") + template_2[4:].replace("_", "")
        if stage.get("orientation") == "2":
            chain_left, chain_right = template[5], template[4]
        else:
            chain_left, chain_right = template[4], template[5]
        query_left = stage.get("receptor", "")
        query_right = stage.get("ligand", "")
        processed = _batch_processed_root(model)
        template_interface_path = (
            processed.parent / "templates" / "interfaces_lists" / f"{template}.json"
        )
        template_sizes: dict = {}
        coverage_status = "missing_template_interface"
        if template_interface_path.is_file():
            try:
                resolved_template_path = template_interface_path.resolve()
                if resolved_template_path not in template_cache:
                    template_cache[resolved_template_path] = json.loads(
                        resolved_template_path.read_text(encoding="utf-8")
                    )
                template_sizes = template_cache[resolved_template_path]
                coverage_status = "available"
            except (OSError, json.JSONDecodeError):
                coverage_status = "invalid_template_interface"
        left_path, left_path_error = _alignment_path(
            processed, f"{query_left}_{template}_{chain_left}.json"
        )
        right_path, right_path_error = _alignment_path(
            processed, f"{query_right}_{template}_{chain_right}.json"
        )
        left: dict = {}
        right: dict = {}
        error_reason = ""
        if left_path_error or right_path_error:
            error_reason = ";".join(
                reason
                for reason in (
                    f"left:{left_path_error}" if left_path_error else "",
                    f"right:{right_path_error}" if right_path_error else "",
                )
                if reason
            )
        else:
            try:
                if left_path not in alignment_cache:
                    alignment_cache[left_path] = json.loads(left_path.read_text(encoding="utf-8"))
                if right_path not in alignment_cache:
                    alignment_cache[right_path] = json.loads(right_path.read_text(encoding="utf-8"))
                left = alignment_cache[left_path]
                right = alignment_cache[right_path]
            except (OSError, json.JSONDecodeError):
                error_reason = "invalid_alignment_json"
        if not error_reason:
            payload_errors = []
            for side, payload in (("left", left), ("right", right)):
                if payload.get("status") != "success":
                    payload_errors.append(f"{side}:alignment_status_not_success")
                if payload.get("aligner") != "GTalign":
                    payload_errors.append(f"{side}:unexpected_aligner")
                if not payload.get("raw_output_sha256"):
                    payload_errors.append(f"{side}:missing_raw_output_sha256")
            error_reason = ";".join(payload_errors)
        actual_model_hash = sha256_file(model) if model.is_file() else ""
        expected_model_hash = stage.get("source_model_sha256", "")
        if not error_reason and not actual_model_hash:
            error_reason = "missing_source_model"
        elif not error_reason and expected_model_hash and actual_model_hash != expected_model_hash:
            error_reason = "source_model_sha256_mismatch"
        left_size = len(template_sizes.get(chain_left, []) or []) or None
        right_size = len(template_sizes.get(chain_right, []) or []) or None
        if coverage_status == "available" and (left_size is None or right_size is None):
            coverage_status = "missing_template_chain"
        left_features = alignment_features(left, left_size)
        right_features = alignment_features(right, right_size)
        rows.append(
            {
                "dataset_row_id": stage.get("dataset_row_id", ""),
                "native_complex_id": stage.get("dataset_row_id", ""),
                "complex": stage.get("complex", ""),
                "query_left": query_left,
                "query_right": query_right,
                "template": template,
                "chain_left": chain_left,
                "chain_right": chain_right,
                "orientation": f"o{stage.get('orientation', '')}",
                "source_pipeline": stage.get("refinement_backend", ""),
                "batch": stage.get("batch", ""),
                "pair_id": stage.get("pair_id", ""),
                "source_gate_status": stage.get("source_gate_status", ""),
                "status": "alignment_failed" if error_reason else "refinement_accepted",
                "error_reason": error_reason,
                "match_count_left": left_features["match_count"],
                "match_count_right": right_features["match_count"],
                "tm_score_left": left_features["tm_score"],
                "tm_score_right": right_features["tm_score"],
                "tm_score_ref_left": _payload_number(left, "tm_score_ref"),
                "tm_score_ref_right": _payload_number(right, "tm_score_ref"),
                "tm_score_query_left": _payload_number(left, "tm_score_query"),
                "tm_score_query_right": _payload_number(right, "tm_score_query"),
                "mapping_count_left": left_features["mapping_count"],
                "mapping_count_right": right_features["mapping_count"],
                "match_coverage_left": left_features["match_coverage"],
                "match_coverage_right": right_features["match_coverage"],
                "template_size_left": left_size or "",
                "template_size_right": right_size or "",
                "template_coverage_status": coverage_status,
                "template_interface_path": (
                    str(template_interface_path.resolve()) if template_interface_path.is_file() else ""
                ),
                "template_interface_sha256": (
                    sha256_file(template_interface_path) if template_interface_path.is_file() else ""
                ),
                "model_complex": str(model.resolve()),
                "source_model_sha256": expected_model_hash or actual_model_hash,
                "observed_source_model_sha256": actual_model_hash,
                "alignment_left_path": str(left_path) if left_path else "",
                "alignment_right_path": str(right_path) if right_path else "",
                "alignment_left_sha256": sha256_file(left_path) if left_path else "",
                "alignment_right_sha256": sha256_file(right_path) if right_path else "",
                "alignment_left_raw_output_sha256": left.get("raw_output_sha256", ""),
                "alignment_right_raw_output_sha256": right.get("raw_output_sha256", ""),
                "alignment_left_aligner": left.get("aligner", ""),
                "alignment_right_aligner": right.get("aligner", ""),
                "alignment_left_status": left.get("status", ""),
                "alignment_right_status": right.get("status", ""),
            }
        )
    return rows


def load_alignments(
    alignment_dir: Path, limit_files: int | None = None
) -> dict[tuple[str, str, str], dict]:
    result = {}
    paths = sorted(alignment_dir.glob("*.json"))
    if limit_files is not None:
        paths = paths[:limit_files]
    for path in paths:
        parts = path.stem.split("_")
        if len(parts) != 3:
            continue
        query, template, chain = parts
        try:
            result[(query, template, chain)] = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
    return result


def load_selected_alignments(
    alignment_dir: Path,
    query_pairs: list[tuple[str, str]],
    templates: set[str],
) -> dict[tuple[str, str, str], dict]:
    """Read only requested filenames; avoids scanning historical output trees."""
    result = {}
    queries = {query for pair in query_pairs for query in pair}
    for query in queries:
        for template in templates:
            if len(template) < 6:
                continue
            for chain in (template[4], template[5]):
                path = alignment_dir / f"{query}_{template}_{chain}.json"
                if not path.is_file():
                    continue
                try:
                    result[(query, template, chain)] = json.loads(path.read_text())
                except (OSError, json.JSONDecodeError):
                    continue
    return result


def build_rows(
    alignment_dir: Path,
    transformation_dir: Path,
    query_pairs: list[tuple[str, str]],
    limit: int | None = None,
    template_filter: set[str] | None = None,
    native_complex_id: str | None = None,
    dataset_row_id: str | None = None,
    refinement_dir: Path | None = None,
) -> list[dict]:
    if template_filter is not None:
        templates = sorted(template_filter)
        if limit is not None:
            if limit <= 0:
                raise ValueError("limit must be positive")
            templates = templates[:limit]
        alignments = load_selected_alignments(alignment_dir, query_pairs, set(templates))
    else:
        # A bounded pilot limits both parsing work and the resulting template set.
        alignments = load_alignments(
            alignment_dir,
            limit_files=(limit * 20 if limit is not None else None),
        )
    rows = []
    templates = sorted({template for _, template, _ in alignments})
    if limit is not None:
        if limit <= 0:
            raise ValueError("limit must be positive")
        templates = templates[:limit]
    for left_query, right_query in query_pairs:
        for template in templates:
            if len(template) < 6:
                continue
            chain_a, chain_b = template[4], template[5]
            for orientation, left_chain, right_chain in (
                ("o1", chain_a, chain_b),
                ("o2", chain_b, chain_a),
            ):
                left = alignments.get((left_query, template, left_chain))
                right = alignments.get((right_query, template, right_chain))
                status = "generated"
                error_reason = None
                if left is None or right is None:
                    status = "alignment_failed"
                    error_reason = "missing_alignment_json"
                left_features = alignment_features(left or {})
                right_features = alignment_features(right or {})
                prefix = f"{template}_{left_query}_{right_query}_{orientation}"
                model_left = transformation_dir / f"{prefix}_L.pdb"
                model_right = transformation_dir / f"{prefix}_R.pdb"
                refined_candidates = (
                    sorted(refinement_dir.glob(f"{prefix}*.pdb"))
                    if refinement_dir is not None
                    else []
                )
                refined_models = select_final_external_models(refined_candidates) or refined_candidates[:1]
                if status == "generated" and not (model_left.exists() and model_right.exists()):
                    status = "transformation_failed"
                    error_reason = "missing_transformed_pair"
                rows.append({
                    "query_left": left_query,
                    "query_right": right_query,
                    "native_complex_id": native_complex_id,
                    "dataset_row_id": dataset_row_id,
                    "template": template,
                    "chain_left": left_chain,
                    "chain_right": right_chain,
                    "orientation": orientation,
                    "source_pipeline": "tmalign_rosetta",
                    "status": status,
                    "error_reason": error_reason,
                    "match_count_left": left_features["match_count"],
                    "match_count_right": right_features["match_count"],
                    "tm_score_left": left_features["tm_score"],
                    "tm_score_right": right_features["tm_score"],
                    "model_left": str(model_left) if model_left.exists() else None,
                    "model_right": str(model_right) if model_right.exists() else None,
                    "model_complex": str(refined_models[0]) if refined_models else None,
                    "source_model_sha256": sha256_file(refined_models[0]) if refined_models else None,
                })
    return rows


def read_pairs(path: Path) -> list[tuple[str, str]]:
    with path.open(newline="") as handle:
        return [(row["Receptor"].strip(), row["Ligand"].strip()) for row in csv.DictReader(handle)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-dir", type=Path, default=Path("processed/alignment"))
    parser.add_argument("--transformation-dir", type=Path, default=Path("processed/transformation"))
    parser.add_argument("--inputs", type=Path, default=Path("inputs.csv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--stage-manifest",
        type=Path,
        help="canonical stage CSV; derive rows from staged models and their batch alignment JSON",
    )
    parser.add_argument("--limit-templates", type=int, default=None)
    parser.add_argument(
        "--query-pair",
        action="append",
        default=[],
        metavar="RECEPTOR,LIGAND",
        help="override inputs.csv; may be repeated",
    )
    parser.add_argument(
        "--template",
        action="append",
        default=[],
        help="restrict to an explicit template; may be repeated",
    )
    parser.add_argument(
        "--native-complex-id",
        default=None,
        help="benchmark complex group identifier to preserve for later labels",
    )
    parser.add_argument(
        "--dataset-row-id",
        default=None,
        help="durable benchmark row identity preserved for canonical score joins",
    )
    parser.add_argument("--refinement-dir", type=Path, default=None)
    args = parser.parse_args()
    if args.stage_manifest:
        rows = build_rows_from_stage_manifest(args.stage_manifest)
    else:
        pairs = (
            [tuple(value.split(",", 1)) for value in args.query_pair]
            if args.query_pair
            else read_pairs(args.inputs)
        )
        if any(len(pair) != 2 or not pair[0] or not pair[1] for pair in pairs):
            parser.error("--query-pair values must be RECEPTOR,LIGAND")
        rows = build_rows(
            args.alignment_dir,
            args.transformation_dir,
            pairs,
            limit=args.limit_templates,
            template_filter=set(args.template) if args.template else None,
            native_complex_id=args.native_complex_id,
            dataset_row_id=args.dataset_row_id,
            refinement_dir=args.refinement_dir,
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["status"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} candidate rows to {args.output}")


if __name__ == "__main__":
    main()
