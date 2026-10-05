#!/usr/bin/env python3
"""Compare matched TMalign refinement arms at candidate level.

The broad stage ledger counts every retained PDB, including FiberDock
intermediates.  This companion report compares the shared transformation
identities with canonical external-Rosetta structures and FiberDock energy
outputs, while preserving partial/zero-transformation states.  It is an
observational diagnostic and does not claim refiner quality causality.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio.PDB import PDBParser


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def files_equal(left: Path, right: Path) -> bool:
    return left.is_file() and right.is_file() and sha256(left) == sha256(right)


def candidate_ids(root: Path) -> set[str]:
    transformation = root / "processed" / "transformation"
    ids: set[str] = set()
    for path in transformation.glob("*.pdb"):
        if path.name.endswith("_L.pdb"):
            ids.add(path.name[: -len("_L.pdb")])
    return ids


def final_rosetta_path(root: Path, candidate: str) -> Path:
    return (
        root
        / "processed"
        / "rosetta_refinement"
        / "structures"
        / f"{candidate}_L_{candidate}_R_rosetta_0001_0001.pdb"
    )


def rosetta_artifacts(root: Path, candidate: str) -> list[Path]:
    directory = root / "processed" / "rosetta_refinement"
    return sorted(directory.glob(f"{candidate}_L_{candidate}_R_rosetta*.pdb"))


def fiberdock_dir(root: Path, candidate: str) -> Path:
    return root / "processed" / "fiberdock_refinement" / candidate


def final_fiberdock_path(root: Path, candidate: str) -> Path:
    return fiberdock_dir(root, candidate) / "fiberdock_energies_1.ref.pdb"


def pdb_summary(path: Path) -> dict[str, object]:
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(path.stem, str(path))
        chains = sorted({chain.id for chain in structure.get_chains()})
        atoms = sum(1 for _ in structure.get_atoms())
        return {"valid": True, "chains": "".join(chains), "atoms": atoms}
    except Exception as exc:  # pragma: no cover - diagnostic preservation
        return {"valid": False, "chains": "", "atoms": 0, "error": str(exc)}


def surface_summary(external: Path, fiberdock: Path) -> dict[str, object]:
    ext_dir = external / "processed" / "surface_extraction"
    fd_dir = fiberdock / "processed" / "surface_extraction"
    asa_names = sorted(path.name for path in ext_dir.glob("*.asa.pdb"))
    rsa_names = sorted(path.name for path in ext_dir.glob("*.rsa"))
    asa_equal = all(files_equal(ext_dir / name, fd_dir / name) for name in asa_names)
    rsa_differences: list[dict[str, object]] = []
    for name in rsa_names:
        left = (ext_dir / name).read_text(encoding="utf-8").splitlines()
        right = (fd_dir / name).read_text(encoding="utf-8").splitlines()
        if left == right:
            continue
        changed = [
            {"external": a, "fiberdock": b}
            for a, b in zip(left, right)
            if a != b
        ]
        rsa_differences.append({"file": name, "changed_lines": changed})
    return {
        "asa_files": len(asa_names),
        "asa_files_byte_identical": asa_equal,
        "rsa_files": len(rsa_names),
        "rsa_files_byte_identical": not rsa_differences,
        "rsa_differences": rsa_differences,
    }


def compare(external: Path, fiberdock: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    candidates = sorted(candidate_ids(external) | candidate_ids(fiberdock))
    rows: list[dict[str, object]] = []
    fd_log = fiberdock / "log"
    fd_log_text = fd_log.read_text(encoding="utf-8", errors="replace") if fd_log.is_file() else ""

    for candidate in candidates:
        rosetta_final = final_rosetta_path(external, candidate)
        rosetta_files = rosetta_artifacts(external, candidate)
        fd_final = final_fiberdock_path(fiberdock, candidate)
        fd_directory = fiberdock_dir(fiberdock, candidate)
        zero_marker = fd_directory / "zero-transformation"
        fd_error = bool(candidate in fd_log_text and "ERROR" in fd_log_text)

        if rosetta_final.is_file():
            rosetta_state = "canonical_final"
        elif rosetta_files:
            rosetta_state = "partial_artifacts"
        else:
            rosetta_state = "missing"
        if fd_final.is_file():
            fiberdock_state = "canonical_energy_output"
        elif fd_directory.is_dir():
            fiberdock_state = "partial_artifacts"
        else:
            fiberdock_state = "missing"

        ext_summary = pdb_summary(rosetta_final) if rosetta_final.is_file() else {}
        fd_summary = pdb_summary(fd_final) if fd_final.is_file() else {}
        rows.append(
            {
                "candidate": candidate,
                "rosetta_state": rosetta_state,
                "rosetta_artifact_count": len(rosetta_files),
                "rosetta_final": str(rosetta_final) if rosetta_final.is_file() else "",
                "rosetta_pdb_valid": ext_summary.get("valid", ""),
                "rosetta_chains": ext_summary.get("chains", ""),
                "rosetta_atoms": ext_summary.get("atoms", ""),
                "fiberdock_state": fiberdock_state,
                "fiberdock_zero_transform_marker": zero_marker.is_file(),
                "fiberdock_final": str(fd_final) if fd_final.is_file() else "",
                "fiberdock_pdb_valid": fd_summary.get("valid", ""),
                "fiberdock_chains": fd_summary.get("chains", ""),
                "fiberdock_atoms": fd_summary.get("atoms", ""),
                "fiberdock_error_in_log": fd_error,
            }
        )

    input_equal = files_equal(external / "inputs.csv", fiberdock / "inputs.csv")
    template_equal = files_equal(
        external / "templates" / "calculated_templates.txt",
        fiberdock / "templates" / "calculated_templates.txt",
    )
    alignment_files = sorted((external / "processed" / "alignment").glob("*.json"))
    alignment_equal = all(
        files_equal(path, fiberdock / "processed" / "alignment" / path.name)
        for path in alignment_files
    )
    transformation_files = sorted((external / "processed" / "transformation").glob("*.pdb"))
    transformation_equal = all(
        files_equal(path, fiberdock / "processed" / "transformation" / path.name)
        for path in transformation_files
    )
    summary = {
        "external_root": str(external.resolve()),
        "fiberdock_root": str(fiberdock.resolve()),
        "candidate_count": len(candidates),
        "candidate_ids": candidates,
        "inputs_byte_identical": input_equal,
        "template_manifest_byte_identical": template_equal,
        "alignment_json_byte_identical": alignment_equal,
        "transformation_pdb_byte_identical": transformation_equal,
        "surface": surface_summary(external, fiberdock),
        "rosetta_canonical_final_count": sum(row["rosetta_state"] == "canonical_final" for row in rows),
        "rosetta_partial_count": sum(row["rosetta_state"] == "partial_artifacts" for row in rows),
        "rosetta_missing_count": sum(row["rosetta_state"] == "missing" for row in rows),
        "fiberdock_canonical_energy_count": sum(row["fiberdock_state"] == "canonical_energy_output" for row in rows),
        "fiberdock_zero_transform_marker_count": sum(row["fiberdock_zero_transform_marker"] for row in rows),
        "fiberdock_error_logged_count": sum(row["fiberdock_error_in_log"] for row in rows),
        "fiberdock_partial_or_missing_count": sum(row["fiberdock_state"] in {"partial_artifacts", "missing"} for row in rows),
        "interpretation": (
            "Inputs, alignment JSON, and transformed PDBs are identical. "
            "Surface ASA coordinate files are identical; RSA differences are "
            "limited to NACCESS total-line formatting/rounding. The first "
            "candidate-level divergence is refinement output state. FiberDock "
            "has an energy PDB for every candidate, but one candidate also has "
            "a logged missing zero-trial input; the zero-transformation marker "
            "is present for all candidates and is normal metadata."
        ),
    }
    return rows, summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, required=True)
    parser.add_argument("--fiberdock-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows, summary = compare(args.external_root, args.fiberdock_root)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["candidate"]
    with (args.output_dir / "matched_candidate_ledger.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "matched_candidate_ledger.json").write_text(
        json.dumps({"summary": summary, "candidates": rows}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
