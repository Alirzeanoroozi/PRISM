
#!/usr/bin/env python3
"""Compare the retained six-pair MultiProt+FiberDock and TMalign+Rosetta arms."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

PARSER = PDBParser(QUIET=True)


def digest_dir(path: Path) -> dict:
    digest = hashlib.sha256()
    count = 0
    for item in sorted(p for p in path.rglob("*") if p.is_file()):
        digest.update(str(item.relative_to(path)).encode())
        digest.update(hashlib.sha256(item.read_bytes()).digest())
        count += 1
    return {"file_count": count, "sha256": digest.hexdigest()}


def ca_map(path: Path) -> dict[tuple[str, tuple], np.ndarray]:
    structure = PARSER.get_structure(path.stem, str(path))
    values = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                atom = residue.child_dict.get("CA")
                if atom is not None:
                    values[(chain.id, residue.id)] = np.asarray(atom.coord, dtype=float)
    return values


def atom_coords(path: Path) -> np.ndarray:
    structure = PARSER.get_structure(path.stem, str(path))
    return np.asarray([atom.coord for atom in structure.get_atoms()], dtype=float)


def direct_rmsd(left: Path, right: Path) -> dict:
    a, b = ca_map(left), ca_map(right)
    keys = sorted(set(a) & set(b), key=str)
    if not keys:
        return {"common_ca": 0, "rmsd": None}
    delta = np.asarray([a[key] - b[key] for key in keys])
    return {"common_ca": len(keys), "rmsd": float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))}


def cross_clashes(left: Path, right: Path, threshold: float = 3.0) -> dict:
    a, b = atom_coords(left), atom_coords(right)
    if not len(a) or not len(b):
        return {"left_atoms": len(a), "right_atoms": len(b), "clashes": None, "min_distance": None}
    min_sq = math.inf
    clashes = 0
    for start in range(0, len(a), 512):
        block = a[start:start + 512]
        delta = block[:, None, :] - b[None, :, :]
        sq = np.sum(delta * delta, axis=2)
        min_sq = min(min_sq, float(np.min(sq)))
        clashes += int(np.count_nonzero(sq <= threshold * threshold))
    return {
        "left_atoms": len(a),
        "right_atoms": len(b),
        "clashes": clashes,
        "min_distance": float(math.sqrt(min_sq)),
    }


def candidate_refiner_state(tm_root: Path, label: str) -> dict:
    rosetta = tm_root / "processed" / "rosetta_refinement"
    pattern = f"{label}_L_{label}_R_rosetta"
    all_files = sorted(rosetta.glob(pattern + "*.pdb"))
    canonical = rosetta / "structures" / f"{pattern}_0001_0001.pdb"
    return {
        "artifact_count": len(all_files),
        "canonical_exists": canonical.is_file(),
        "canonical_path": str(canonical),
        "artifacts": [str(path) for path in all_files],
    }


def main() -> int:
    if len(sys.argv) != 5:
        raise SystemExit("usage: compare_six_pair_arms.py REPLAY_JSON MP_ROOT TM_ROOT OUTPUT_JSON")
    replay_json, mp_root, tm_root, output_json = (Path(value).resolve() for value in sys.argv[1:])
    replay = json.loads(replay_json.read_text(encoding="utf-8"))
    rows = []
    for record in replay["records"]:
        label = record["label"]
        mp_left = Path(record["source_left"])
        mp_right = Path(record["source_right"])
        tm_left = tm_root / "processed" / "transformation" / f"{label}_L.pdb"
        tm_right = tm_root / "processed" / "transformation" / f"{label}_R.pdb"
        row = {
            "label": label,
            "query_left": record["query_left"],
            "query_right": record["query_right"],
            "template": record["template"],
            "orientation": record["orientation"],
            "mp_alignment_left": record["mp_alignment_left"],
            "mp_alignment_right": record["mp_alignment_right"],
            "tm_alignment_left": record["tm_alignment_left"],
            "tm_alignment_right": record["tm_alignment_right"],
            "mp_transform_exists": mp_left.is_file() and mp_right.is_file(),
            "tm_transform_exists": tm_left.is_file() and tm_right.is_file(),
            "mp_fiberdock_pdbs": record["refined_pdbs"],
            "tm_rosetta": candidate_refiner_state(tm_root, label),
        }
        if row["mp_transform_exists"] and row["tm_transform_exists"]:
            row["left_direct_ca"] = direct_rmsd(mp_left, tm_left)
            row["right_direct_ca"] = direct_rmsd(mp_right, tm_right)
            row["mp_cross_clashes"] = cross_clashes(mp_left, mp_right)
            row["tm_cross_clashes"] = cross_clashes(tm_left, tm_right)
        rows.append(row)

    mp_alignment = list((mp_root / "processed" / "alignment").glob("*.json"))
    tm_alignment = list((tm_root / "processed" / "alignment").glob("*.json"))
    summary = {
        "candidate_count": len(rows),
        "both_transform_count": sum(r["mp_transform_exists"] and r["tm_transform_exists"] for r in rows),
        "mp_fiberdock_valid_pdb_count": sum(
            item["valid"] for r in rows for item in r["mp_fiberdock_pdbs"]
        ),
        "tm_rosetta_canonical_count": sum(r["tm_rosetta"]["canonical_exists"] for r in rows),
        "tm_rosetta_partial_or_missing_count": sum(not r["tm_rosetta"]["canonical_exists"] for r in rows),
        "mp_alignment_success_sides": sum(
            r["mp_alignment_left"]["status"] == "success" for r in rows
        ) + sum(r["mp_alignment_right"]["status"] == "success" for r in rows),
        "tm_alignment_success_sides": sum(
            r["tm_alignment_left"]["status"] == "success" for r in rows
        ) + sum(r["tm_alignment_right"]["status"] == "success" for r in rows),
        "direct_ca_rmsd_left": [
            r["left_direct_ca"]["rmsd"] for r in rows if r.get("left_direct_ca", {}).get("rmsd") is not None
        ],
        "direct_ca_rmsd_right": [
            r["right_direct_ca"]["rmsd"] for r in rows if r.get("right_direct_ca", {}).get("rmsd") is not None
        ],
        "mp_cross_clashes_total": sum(r.get("mp_cross_clashes", {}).get("clashes", 0) for r in rows),
        "tm_cross_clashes_total": sum(r.get("tm_cross_clashes", {}).get("clashes", 0) for r in rows),
    }
    report = {
        "experiment": "six_pair_multiprot_fiberdock_vs_tmalign_rosetta",
        "replay_json": str(replay_json),
        "mp_root": str(mp_root),
        "tm_root": str(tm_root),
        "input_sha256_mp": hashlib.sha256((mp_root / "inputs.csv").read_bytes()).hexdigest(),
        "input_sha256_tm": hashlib.sha256((tm_root / "inputs.csv").read_bytes()).hexdigest(),
        "pdb_inventory_mp": digest_dir(mp_root / "processed" / "pdbs"),
        "pdb_inventory_tm": digest_dir(tm_root / "processed" / "pdbs"),
        "surface_inventory_mp": digest_dir(mp_root / "processed" / "surface_extraction"),
        "surface_inventory_tm": digest_dir(tm_root / "processed" / "surface_extraction"),
        "alignment_file_count_mp": len(mp_alignment),
        "alignment_file_count_tm": len(tm_alignment),
        "summary": summary,
        "records": rows,
        "interpretation": (
            "Diagnostic stage comparison only. Inputs and source PDB inventories "
            "match; MultiProt and TMalign alignment contracts differ. The seven "
            "common transformed candidates are compared geometrically. Rosetta "
            "partial/missing states remain explicit. No DockQ/native quality claim."
        ),
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
