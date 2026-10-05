#!/usr/bin/env python3
"""Map completed current/legacy models to benchmark rows and score them."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path
import sys
import re

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def normalize_target_id(target):
    value = str(target).strip().replace(" ", "")
    if len(value) < 4:
        raise ValueError(f"invalid target identifier: {target!r}")
    suffix = re.sub(r"\([^)]*\)", "", value[4:].replace("_", ""))
    return value[:4].lower() + "".join(sorted(set(c for c in suffix if c.isalnum())))


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def resolve_n_cpu(n_cpu: int | None = None) -> int:
    value = n_cpu if n_cpu is not None else os.environ.get("DOCKQ_N_CPU", "1")
    try:
        value = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("DockQ n_cpu must be a positive integer") from exc
    if value < 1:
        raise ValueError("DockQ n_cpu must be a positive integer")
    return value


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def model_chain_groups(path: Path, left_hint: str, right_hint: str):
    chains = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip() or "_"
                if chain not in chains:
                    chains.append(chain)
    left_count, right_count = len(left_hint), len(right_hint)
    if left_hint and right_hint and set(left_hint).issubset(chains) and set(right_hint).issubset(chains) and not set(left_hint).intersection(right_hint):
        return left_hint, right_hint, "source_chain_ids_present_in_model"
    if left_count and right_count and len(chains) >= left_count + right_count:
        return "".join(chains[:left_count]), "".join(chains[left_count:left_count + right_count]), "inferred_from_target_chain_counts"
    # Legacy FiberDock and single-chain current outputs normally have one
    # chain per partner; retain an explicit ambiguity label otherwise.
    if len(chains) == 2:
        return chains[0], chains[1], "inferred_from_two_model_chains"
    return "", "", "chain_mapping_ambiguous"


def benchmark_groups(complex_id: str):
    value = complex_id.strip().split("_", 1)
    if len(value) != 2 or ":" not in value[1]:
        return "", ""
    receptor, ligand = value[1].split(":", 1)
    return receptor.strip(), ligand.strip()


def discover_models(current_roots: list[Path], legacy_root: Path):
    models = []
    for root in current_roots:
        # The external refiner has emitted both the historical direct output
        # (`rosetta_refinement/*_0001.pdb`) and the newer copied structure
        # layout (`rosetta_refinement/structures/*_0001_0001.pdb`).  Discover
        # both, but never score the unnumbered convenience copy.
        candidates = list(root.glob("batch_*/processed/rosetta_refinement/*_0001.pdb"))
        candidates.extend(root.glob("batch_*/processed/rosetta_refinement/structures/*_0001_0001.pdb"))
        seen: set[Path] = set()
        for path in candidates:
            if path in seen or path.name.endswith("_rosetta.pdb"):
                continue
            seen.add(path)
            models.append(("tmalign_rosetta", path))
    for path in legacy_root.glob("batch_*/jobs/batch_*/fiberdock/*.ref.pdb"):
        models.append(("multiprot_fiberdock", path))
    return sorted(models, key=lambda item: str(item[1]))


def find_pair(path: Path, batch_rows: dict[str, list[dict[str, str]]]):
    match = re.search(r"batch_(\d{4})", str(path))
    if not match:
        return None
    rows = batch_rows.get(f"batch_{match.group(1)}", [])
    name = path.name.lower()
    hits = []
    for row in rows:
        left = normalize_target_id(row["Receptor"]).lower()
        right = normalize_target_id(row["Ligand"]).lower()
        left_match = re.search(rf"(?<![a-z0-9]){re.escape(left)}(?![a-z0-9])", name)
        right_match = re.search(rf"(?<![a-z0-9]){re.escape(right)}(?![a-z0-9])", name)
        if left_match and right_match and left_match.start() < right_match.start():
            hits.append(row)
    return hits[0] if len(hits) == 1 else None


def build_manifest(manifest: Path, current_roots: list[Path], legacy_root: Path, output: Path):
    rows = read_csv(manifest)
    by_batch = {}
    for row in rows:
        batch_number = (int(row["source_row"]) - 1) // 10 + 1 if row["benchmark_set"] == "rigid" else None
        # Use the prepared batch files for exact provenance rather than
        # reconstructing batch numbering across benchmark sets.
    batch_root = manifest.parent
    for batch_path in sorted(batch_root.glob("batch_*/inputs.csv")):
        by_batch[batch_path.parent.name] = read_csv(batch_path)

    result = []
    for pipeline, model in discover_models(current_roots, legacy_root):
        row = find_pair(model, by_batch)
        if row is None:
            result.append({"pipeline": pipeline, "model_path": str(model), "status": "unmatched_model"})
            continue
        native = Path("benchmark/prism_processed/results") / f"native_bound_complexes_t_{row['benchmark_set']}" / f"{row['complex'][:4].lower()}.pdb"
        native_receptor, native_ligand = benchmark_groups(row["complex"])
        left_hint = normalize_target_id(row["Receptor"])[4:]
        right_hint = normalize_target_id(row["Ligand"])[4:]
        model_receptor, model_ligand, mapping_status = model_chain_groups(model, left_hint, right_hint)
        result.append({
            "pipeline": pipeline,
            "pair_id": row["pair_id"],
            "benchmark_set": row["benchmark_set"],
            "complex": row["complex"],
            "model_path": str(model),
            "native_path": str(native),
            "model_receptor": model_receptor,
            "model_ligand": model_ligand,
            "native_receptor": native_receptor,
            "native_ligand": native_ligand,
            "mapping_status": mapping_status,
            "status": "ready" if native.exists() and model_receptor and model_ligand else "not_scoreable",
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    write_rows(output, result or [{"status": ""}])
    print(f"wrote {len(result)} model rows to {output}")


def score_manifest(
    model_manifest: Path,
    output: Path,
    score_python: Path,
    repo_root: Path,
    dockq_json_dir: Path | None = None,
    *,
    dockq_no_align: bool = True,
    n_cpu: int | None = None,
):
    dockq_cpu = resolve_n_cpu(n_cpu)
    rows = read_csv(model_manifest)
    scored = []
    for row in rows:
        if row.get("status") != "ready":
            scored.append(row | {
                "irmsd": "", "dockq": "", "dockq_global": "", "dockq_interface_count": "",
                "dockq_raw_json_path": "", "dockq_raw_json_sha256": "",
                "dockq_json_status": "not_attempted", "dockq_n_cpu": dockq_cpu,
                "score_error": row.get("status", "not_scoreable"),
            })
            continue
        out_csv = output.parent / (Path(row["model_path"]).stem + ".score.csv")
        cmd = [
            str(score_python), str(repo_root / "benchmark/scripts/score_single_prism_pair.py"),
            row["model_path"], row["native_path"],
            "--model-receptor", row["model_receptor"], "--model-ligand", row["model_ligand"],
            "--native-receptor", row["native_receptor"], "--native-ligand", row["native_ligand"],
            "--score-python", str(score_python), "--out-csv", str(out_csv),
            "--n-cpu", str(dockq_cpu),
        ]
        if dockq_no_align:
            cmd.append("--dockq-no-align")
        if dockq_json_dir is not None:
            cmd.extend(["--dockq-json-dir", str(dockq_json_dir)])
        row = dict(row)
        try:
            proc = subprocess.run(cmd, cwd=str(repo_root), text=True, capture_output=True, timeout=300)
        except subprocess.TimeoutExpired:
            row.update({
                "irmsd": "", "dockq": "", "dockq_global": "", "dockq_interface_count": "",
                "dockq_raw_json_path": "", "dockq_raw_json_sha256": "",
                "dockq_json_status": "timeout", "dockq_n_cpu": dockq_cpu,
                "score_error": "scoring timed out after 300s",
            })
            scored.append(row)
            continue
        if proc.returncode == 0 and out_csv.exists():
            scored_row = read_csv(out_csv)[0]
            for field in (
                "irmsd", "dockq", "dockq_global", "dockq_interface_count",
                "dockq_raw_json_path", "dockq_raw_json_sha256", "dockq_mapping",
                "mapping_validation_status", "mapping_validation_errors",
                "dockq_json_status", "dockq_n_cpu",
            ):
                row[field] = scored_row.get(field, "")
            score_error = scored_row.get("error_irmsd", "") or scored_row.get("error_dockq", "")
            # Preserve the full diagnostic in the scorer's raw sidecar while
            # keeping the tabular benchmark report bounded and readable.
            row["score_error"] = score_error[:2000]
        else:
            row.update({
                "irmsd": "",
                "dockq": "",
                "dockq_json_status": "failed",
                "dockq_n_cpu": dockq_cpu,
                "score_error": (proc.stderr or proc.stdout)[-2000:],
            })
        scored.append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_rows(output, scored or [{"status": ""}])
    print(f"wrote {len(scored)} scored model rows to {output}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--current-root", type=Path, action="append", default=[])
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--score-output", type=Path)
    parser.add_argument("--score-python", type=Path)
    parser.add_argument("--dockq-json-dir", type=Path)
    parser.add_argument(
        "--n-cpu",
        type=int,
        default=None,
        help="Explicit CPU count passed to DockQ; defaults to DOCKQ_N_CPU or 1",
    )
    parser.add_argument(
        "--dockq-align",
        action="store_true",
        help="allow DockQ's alignment step; default is the strict no-align contract",
    )
    args = parser.parse_args()
    build_manifest(args.manifest, args.current_root, args.legacy_root, args.output)
    if args.score_output:
        if not args.score_python:
            parser.error("--score-python is required with --score-output")
        score_manifest(
            args.output,
            args.score_output,
            args.score_python,
            REPO_ROOT,
            args.dockq_json_dir,
            dockq_no_align=not args.dockq_align,
            n_cpu=args.n_cpu,
        )


if __name__ == "__main__":
    main()
