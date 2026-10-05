
#!/usr/bin/env python3
"""Diagnose the FiberDock parameter/output filename contract.

This is intentionally isolated from src/fiberdock_refinement.py. It reads
completed FiberDock work directories, derives the declared energy path from
fd_params.txt, compares it with the current implementation's fd_params.ref
guess, parses the declared .ref solution, and validates refined PDB artifacts.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def declared_value(params_path: Path, key: str) -> str | None:
    prefix = key + " "
    for line in params_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            return line[len(prefix):].strip()
    return None


def parse_energy(ref_path: Path) -> dict | None:
    if not ref_path.is_file():
        return None
    for line in ref_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) >= 3 and parts[0].isdigit():
            try:
                global_energy = float(parts[1])
            except ValueError:
                global_energy = parts[1]
            return {
                "solution": int(parts[0]),
                "global_energy": global_energy,
                "raw_line": line,
            }
    return None


def pdb_integrity(path: Path) -> dict:
    result = {"path": str(path), "exists": path.is_file(), "valid": False}
    if not path.is_file():
        return result
    try:
        from Bio.PDB import PDBParser
        structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
        atoms = 0
        residues = 0
        chains = []
        for model in structure:
            for chain in model:
                chains.append(chain.id)
                for residue in chain:
                    residues += residue.id[0] == " "
                    atoms += len(residue.child_list)
        result.update({
            "valid": bool(atoms and chains),
            "atoms": atoms,
            "residues": residues,
            "chains": "".join(chains),
        })
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def inspect_workdir(workdir: Path) -> dict:
    params = workdir / "fd_params.txt"
    declared = declared_value(params, "energiesOutFileName")
    declared_ref = Path(declared + ".ref") if declared else None
    current_guess = workdir / "fd_params.ref"
    pdbs = sorted(workdir.glob("fiberdock_energies*.ref.pdb"))
    return {
        "workdir": str(workdir),
        "params_file": str(params),
        "declared_energy_stem": declared,
        "declared_ref": str(declared_ref) if declared_ref else None,
        "declared_ref_exists": bool(declared_ref and declared_ref.is_file()),
        "current_parser_guess": str(current_guess),
        "current_parser_guess_exists": current_guess.is_file(),
        "declared_solution": parse_energy(declared_ref) if declared_ref else None,
        "refined_pdbs": [pdb_integrity(path) for path in pdbs],
        "contract_status": (
            "declared_output_found_current_guess_wrong"
            if declared_ref and declared_ref.is_file() and not current_guess.is_file()
            else "declared_output_missing"
            if not declared_ref or not declared_ref.is_file()
            else "current_guess_matches_declared_output"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path, help="Replay root or one FiberDock work directory")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if (root / "fd_params.txt").is_file():
        workdirs = [root]
    else:
        workdirs = sorted(path.parent for path in root.rglob("fd_params.txt"))
    report = {
        "experiment": "fiberdock_output_contract",
        "root": str(root),
        "workdir_count": len(workdirs),
        "records": [inspect_workdir(path) for path in workdirs],
    }
    report["summary"] = {
        "declared_outputs_found": sum(item["declared_ref_exists"] for item in report["records"]),
        "current_guess_found": sum(item["current_parser_guess_exists"] for item in report["records"]),
        "corrected_energy_values": [
            item["declared_solution"]["global_energy"]
            for item in report["records"]
            if item["declared_solution"] is not None
        ],
        "valid_refined_pdbs": sum(
            pdb["valid"] for item in report["records"] for pdb in item["refined_pdbs"]
        ),
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
