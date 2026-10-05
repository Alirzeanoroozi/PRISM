#!/usr/bin/env python3
"""Matched MultiProt/TMalign transform and refiner replay.

The two arms use the same query pair, template, source PDBs, transform code,
clash thresholds, and refiner entry points. MultiProt has an additional
diagnostic candidate (1ahwBC) whose proxy score fails; the common 1ahwAF arm is
the causal comparison target.
"""

from __future__ import annotations

import importlib
import json
import os
import shutil
import sys
import traceback
from pathlib import Path

QUERY_LEFT = "1fgnHL"
QUERY_RIGHT = "1tfhA"
COMMON_CANDIDATE = ("1ahwAF", "A", "F")
MULTIPROT_EXTRA = ("1ahwBC", "B", "C")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def pdb_summary(path: Path) -> dict:
    from Bio.PDB import PDBParser

    if not path.is_file():
        return {"path": str(path), "exists": False, "valid": False}
    structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
    chains = []
    atoms = 0
    residues = 0
    for model in structure:
        for chain in model:
            chains.append(chain.id)
            for residue in chain:
                atoms += len(residue.child_list)
                residues += residue.id[0] == " "
    return {
        "path": str(path),
        "exists": True,
        "valid": bool(chains and atoms),
        "chains": "".join(chains),
        "residues": residues,
        "atoms": atoms,
    }


def stage_arm(arm_root: Path, panel_root: Path, project_root: Path) -> None:
    (arm_root / "processed").mkdir(parents=True, exist_ok=True)
    (arm_root / "external_tools").mkdir(parents=True, exist_ok=True)
    (arm_root / "processed" / "pdbs").symlink_to(
        panel_root / "processed" / "pdbs", target_is_directory=True
    )
    (arm_root / "templates").symlink_to(panel_root / "templates", target_is_directory=True)
    shutil.copytree(
        project_root / "external_tools" / "fiberdock",
        arm_root / "external_tools" / "fiberdock",
        symlinks=True,
    )
    (arm_root / "processed" / "alignment").mkdir(parents=True, exist_ok=True)
    (arm_root / "processed" / "transformation").mkdir(parents=True, exist_ok=True)


def copy_alignment(source_root: Path, arm_root: Path, template: str, chain_left: str, chain_right: str) -> None:
    for query, chain in ((QUERY_LEFT, chain_left), (QUERY_RIGHT, chain_right)):
        source = source_root / "processed" / "alignment" / f"{query}_{template}_{chain}.json"
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, arm_root / "processed" / "alignment" / source.name)


def run_transformations(arm_root: Path, alignment_source: Path, candidates: tuple[tuple[str, str, str], ...], gate: str) -> tuple[list[tuple[str, str]], list[dict]]:
    os.chdir(arm_root)
    os.environ["PRISM_FILTER_MODE"] = "geometry_only_experimental"
    os.environ["PRISM_FIBERDOCK_DIR"] = str(arm_root / "external_tools" / "fiberdock")
    if str(Path(__file__).resolve().parents[2]) not in sys.path:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from src import transformation as tr

    tr.passed_pairs.clear()
    tr.template_size.clear()
    passed: list[tuple[str, str]] = []
    records = []
    for template, left_chain, right_chain in candidates:
        copy_alignment(alignment_source, arm_root, template, left_chain, right_chain)
        left = load_json(arm_root / "processed" / "alignment" / f"{QUERY_LEFT}_{template}_{left_chain}.json")
        right = load_json(arm_root / "processed" / "alignment" / f"{QUERY_RIGHT}_{template}_{right_chain}.json")
        status = tr.create_transformed_pair(template, QUERY_LEFT, QUERY_RIGHT, left, right, passed, "o1")
        left_path = arm_root / "processed" / "transformation" / f"{template}_{QUERY_LEFT}_{QUERY_RIGHT}_o1_L.pdb"
        right_path = arm_root / "processed" / "transformation" / f"{template}_{QUERY_LEFT}_{QUERY_RIGHT}_o1_R.pdb"
        records.append({
            "template": template,
            "orientation": "o1",
            "gate": gate,
            "left_match_count": left.get("match_count", 0),
            "right_match_count": right.get("match_count", 0),
            "left_tm_score": left.get("tm_score", 0.0),
            "right_tm_score": right.get("tm_score", 0.0),
            "left_rmsd": left.get("rmsd"),
            "right_rmsd": right.get("rmsd"),
            "transform_status": status,
            "left_integrity": pdb_summary(left_path),
            "right_integrity": pdb_summary(right_path),
        })
    return passed, records


def run_refiners(arm_root: Path, passed: list[tuple[str, str]]) -> dict:
    os.chdir(arm_root)
    results: dict[str, dict] = {}
    os.environ["PRISM_FIBERDOCK_DIR"] = str(arm_root / "external_tools" / "fiberdock")
    try:
        import src.fiberdock_refinement as fd
        fd = importlib.reload(fd)
        results["fiberdock"] = {"status": "completed", "results": fd.refine_pairs(passed)}
    except Exception as exc:
        results["fiberdock"] = {"status": "failed", "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()}
    try:
        import src.rosetta_refinement as rr
        rr = importlib.reload(rr)
        rr.refiner(passed)
        results["external_rosetta"] = {"status": "completed", "results": []}
    except Exception as exc:
        results["external_rosetta"] = {"status": "failed", "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()}
    for name, relative in (
        ("fiberdock", Path("processed/fiberdock_refinement")),
        ("external_rosetta", Path("processed/rosetta_refinement")),
    ):
        root = arm_root / relative
        paths = sorted(root.rglob("*.pdb")) if root.exists() else []
        results[name]["pdb_count"] = len(paths)
        results[name]["structures"] = [pdb_summary(path) for path in paths]
    return results


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: run_matched_aligner_refiner_replay.py PANEL_ROOT TMALIGN_ROOT OUTPUT_ROOT")
    panel_root, tmalign_root, output_root = (Path(value).resolve() for value in sys.argv[1:])
    project_root = Path(__file__).resolve().parents[2]
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"output root is not empty: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    configs = {
        "multiprot": {
            "alignment_source": panel_root,
            "candidates": (COMMON_CANDIDATE, MULTIPROT_EXTRA),
            "gate": "diagnostic status + match_count + match_percentage; MultiProt tm_score ignored",
        },
        "tmalign": {
            "alignment_source": tmalign_root,
            "candidates": (COMMON_CANDIDATE,),
            "gate": "current TMalign side/pair gate; selected 1ahwAF/o1 from exact calibration",
        },
    }
    summary = {
        "experiment": "matched_aligner_refiner_replay",
        "panel_root": str(panel_root),
        "tmalign_root": str(tmalign_root),
        "output_root": str(output_root),
        "queries": [QUERY_LEFT, QUERY_RIGHT],
        "arms": {},
        "interpretation": (
            "Diagnostic only. The common 1ahwAF/o1 record uses identical query, "
            "template, source PDBs, transformation code, clash thresholds, and "
            "refiner entry points. MultiProt additionally replays 1ahwBC/o1 "
            "because its proxy score fails despite paired alignment success."
        ),
    }
    for name, config in configs.items():
        arm_root = output_root / name
        stage_arm(arm_root, panel_root, project_root)
        passed, records = run_transformations(
            arm_root, config["alignment_source"], config["candidates"], config["gate"]
        )
        refiners = run_refiners(arm_root, passed)
        summary["arms"][name] = {
            "candidate_count": len(records),
            "passed_pair_count": len(passed),
            "candidates": records,
            "refiners": refiners,
        }
    output_root.joinpath("matched_aligner_refiner_replay.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
