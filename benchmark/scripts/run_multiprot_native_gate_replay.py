#!/usr/bin/env python3
"""Diagnostic MultiProt native-gate transform/refiner replay.

This experiment deliberately bypasses only the uncalibrated MultiProt TM-score
proxy. It reuses the current transformation writer and refiner entry points in
an isolated workdir; it is not a production threshold change.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import traceback
from pathlib import Path


PAIRS = (
    ("1ahwAF", "A", "F"),
    ("1ahwBC", "B", "C"),
)
QUERY_LEFT = "1fgnHL"
QUERY_RIGHT = "1tfhA"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _structure_summary(path: Path) -> dict:
    from Bio.PDB import PDBParser

    if not path.is_file():
        return {"path": str(path), "exists": False, "valid": False}
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    chains = []
    atoms = 0
    residues = 0
    for model in structure:
        for chain in model:
            chains.append(chain.id)
            for residue in chain:
                if residue.id[0] == " ":
                    residues += 1
                atoms += len(residue.child_list)
    return {
        "path": str(path),
        "exists": True,
        "valid": bool(chains and atoms),
        "chains": "".join(chains),
        "residues": residues,
        "atoms": atoms,
    }


def _stage(output_root: Path, panel_root: Path, project_root: Path) -> None:
    (output_root / "processed").mkdir(parents=True, exist_ok=True)
    (output_root / "external_tools").mkdir(parents=True, exist_ok=True)
    for source, destination in (
        (panel_root / "processed" / "pdbs", output_root / "processed" / "pdbs"),
        (panel_root / "templates", output_root / "templates"),
    ):
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() or destination.is_symlink():
            raise RuntimeError(f"refusing to replace existing staged path: {destination}")
        destination.symlink_to(source, target_is_directory=True)
    shutil.copy2(panel_root / "inputs.csv", output_root / "inputs.csv")
    shutil.copytree(
        project_root / "external_tools" / "fiberdock",
        output_root / "external_tools" / "fiberdock",
        symlinks=True,
    )
    alignment_dir = output_root / "processed" / "alignment"
    alignment_dir.mkdir(parents=True, exist_ok=True)
    source_alignment = panel_root / "processed" / "alignment"
    for template, left_chain, right_chain in PAIRS:
        for query, chain in (
            (QUERY_LEFT, left_chain),
            (QUERY_RIGHT, right_chain),
        ):
            source = source_alignment / f"{query}_{template}_{chain}.json"
            if not source.is_file():
                raise FileNotFoundError(source)
            shutil.copy2(source, alignment_dir / source.name)


def _run_refiner(name: str, passed_pairs: list[tuple[str, str]]) -> dict:
    result = {"status": "not_started", "results": [], "error": ""}
    try:
        if name == "fiberdock":
            from src.fiberdock_refinement import refine_pairs
        else:
            from src.rosetta_refinement import refiner as refine_pairs
        returned = refine_pairs(passed_pairs)
        result["status"] = "completed"
        result["results"] = returned if returned is not None else []
    except Exception as exc:  # preserve diagnostic failure state
        result["status"] = "failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()
    return result


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_multiprot_native_gate_replay.py PANEL_ROOT OUTPUT_ROOT")
    panel_root = Path(sys.argv[1]).resolve()
    output_root = Path(sys.argv[2]).resolve()
    project_root = Path(__file__).resolve().parents[2]
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"output root is not empty: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    _stage(output_root, panel_root, project_root)
    os.chdir(output_root)
    os.environ["PRISM_FILTER_MODE"] = "geometry_only_experimental"
    os.environ["PRISM_FIBERDOCK_DIR"] = str(output_root / "external_tools" / "fiberdock")
    os.environ["PRISM_STAGE_STATUS_PATH"] = str(output_root / "stage_status.jsonl")
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    from src import transformation as tr

    passed_pairs: list[tuple[str, str]] = []
    candidates: list[dict] = []
    for template, left_chain, right_chain in PAIRS:
        left = _load(output_root / "processed" / "alignment" / f"{QUERY_LEFT}_{template}_{left_chain}.json")
        right = _load(output_root / "processed" / "alignment" / f"{QUERY_RIGHT}_{template}_{right_chain}.json")
        status = tr.create_transformed_pair(
            template, QUERY_LEFT, QUERY_RIGHT, left, right, passed_pairs, "o1"
        )
        candidates.append(
            {
                "template": template,
                "orientation": "o1",
                "gate": "status + match_count + match_percentage; MultiProt tm_score ignored diagnostically",
                "left_match_count": left.get("match_count", 0),
                "right_match_count": right.get("match_count", 0),
                "left_proxy_tm_score": left.get("tm_score", 0.0),
                "right_proxy_tm_score": right.get("tm_score", 0.0),
                "left_rmsd": left.get("rmsd"),
                "right_rmsd": right.get("rmsd"),
                "transform_status": status,
                "left_transformed": str(output_root / "processed" / "transformation" / f"{template}_{QUERY_LEFT}_{QUERY_RIGHT}_o1_L.pdb"),
                "right_transformed": str(output_root / "processed" / "transformation" / f"{template}_{QUERY_LEFT}_{QUERY_RIGHT}_o1_R.pdb"),
            }
        )

    transform_summaries = []
    for candidate in candidates:
        candidate["left_integrity"] = _structure_summary(Path(candidate["left_transformed"]))
        candidate["right_integrity"] = _structure_summary(Path(candidate["right_transformed"]))
        transform_summaries.append(candidate)

    refiner_results = {
        "fiberdock": _run_refiner("fiberdock", passed_pairs),
        "external_rosetta": _run_refiner("external_rosetta", passed_pairs),
    }
    structure_roots = {
        "fiberdock": output_root / "processed" / "fiberdock_refinement",
        "external_rosetta": output_root / "processed" / "rosetta_refinement",
    }
    for name, root in structure_roots.items():
        paths = sorted(root.rglob("*.pdb")) if root.exists() else []
        refiner_results[name]["pdb_count"] = len(paths)
        refiner_results[name]["structures"] = [_structure_summary(path) for path in paths]

    summary = {
        "experiment": "multiprot_native_gate_replay",
        "panel_root": str(panel_root),
        "output_root": str(output_root),
        "queries": [QUERY_LEFT, QUERY_RIGHT],
        "candidate_count": len(candidates),
        "passed_pair_count": len(passed_pairs),
        "candidates": transform_summaries,
        "refiners": refiner_results,
        "interpretation": (
            "Diagnostic only. The MultiProt tm_score proxy was ignored for this "
            "replay so the two paired-success orientations could reach the "
            "existing transform/refiner contracts. This is not a production "
            "threshold change or a quality claim."
        ),
    }
    (output_root / "native_gate_replay.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
