import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import freesasa
import rust_sasa_python

CURRENT_DIR = Path(__file__).resolve().parent
DIFFMASIF_DIR = CURRENT_DIR.parent / "diffmasif_replacement"
if str(DIFFMASIF_DIR) not in sys.path:
    sys.path.insert(0, str(DIFFMASIF_DIR))

from compare_with_naccess import parse_rsa_file, summarize_comparison
from export_residue_accessibility import STANDARD_ASA


def resolve_structure(structure_kind: str, structure_id: str) -> Tuple[Path, List[str]]:
    if structure_kind == "target":
        if len(structure_id) < 5:
            raise ValueError("Target ids must include a chain, for example `1FGNH`.")
        return Path(f"processed/pdbs/{structure_id[:4].lower()}.pdb"), [structure_id[4]]
    if len(structure_id) < 6:
        raise ValueError("Template ids must include two chains, for example `2AI9AB`.")
    return Path(f"templates/pdbs/{structure_id[:4].lower()}.pdb"), [structure_id[4], structure_id[5]]


def extract_selected_chains(source_pdb: Path, chains: Iterable[str], output_pdb: Path) -> None:
    allowed = set(chains)
    kept = []
    with source_pdb.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            chain_id = line[21].strip()
            if chain_id not in allowed:
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            kept.append(line.rstrip("\n"))
    if not kept:
        raise ValueError(f"No ATOM records kept for chains {sorted(allowed)} in {source_pdb}")
    output_pdb.write_text("\n".join(kept) + "\n")


def relative_asa_key(resname: str, residue_number: str, chain_id: str) -> str:
    return f"{resname}_{residue_number}_{chain_id}"


def run_naccess(clean_pdb: Path, allowed_chains: Iterable[str]) -> Dict[str, float]:
    naccess_exec = Path("external_tools/naccess/naccess").resolve()
    if not naccess_exec.exists():
        raise FileNotFoundError(f"Naccess executable not found: {naccess_exec}")

    with tempfile.TemporaryDirectory(prefix="asa_naccess_", dir=".") as tmpdir:
        tmpdir_path = Path(tmpdir)
        local_pdb = tmpdir_path / clean_pdb.name
        shutil.copy2(clean_pdb, local_pdb)
        result = subprocess.run(
            [str(naccess_exec), local_pdb.name],
            cwd=str(tmpdir_path),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            check=False,
        )
        rsa_candidates = [
            tmpdir_path / f"{local_pdb.stem}.rsa",
            tmpdir_path / f"{local_pdb.stem.lower()}.rsa",
            tmpdir_path / f"{local_pdb.stem[:4].lower()}.rsa",
        ]
        rsa_path = next((candidate for candidate in rsa_candidates if candidate.exists()), None)
        if result.returncode != 0 or rsa_path is None:
            raise RuntimeError(
                f"Naccess failed for {clean_pdb.name} (exit code {result.returncode}).\n{result.stdout}"
            )
        return parse_rsa_file(str(rsa_path), allowed_chains=allowed_chains)


def run_freesasa(clean_pdb: Path, allowed_chains: Iterable[str]) -> Dict[str, float]:
    allowed = set(allowed_chains)
    structure = freesasa.Structure(str(clean_pdb))
    result = freesasa.calc(structure)
    residue_areas = result.residueAreas()
    relative: Dict[str, float] = {}
    for chain_id, residues in residue_areas.items():
        if chain_id not in allowed:
            continue
        for residue_number, area in residues.items():
            resname = area.residueType.upper()
            standard = STANDARD_ASA.get(resname)
            if not standard:
                continue
            relative[relative_asa_key(resname, residue_number, chain_id)] = area.total * 100.0 / standard
    return relative


def run_rustsasa(clean_pdb: Path, allowed_chains: Iterable[str]) -> Dict[str, float]:
    allowed = set(allowed_chains)
    residues = rust_sasa_python.calculate_residue_sasa(str(clean_pdb))
    relative: Dict[str, float] = {}
    for residue in residues:
        if residue.chain_id not in allowed:
            continue
        resname = residue.residue_name.upper()
        standard = STANDARD_ASA.get(resname)
        if not standard:
            continue
        relative[
            relative_asa_key(resname, str(residue.residue_number), residue.chain_id)
        ] = float(residue.sasa) * 100.0 / standard
    return relative


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure-kind", required=True, choices=["target", "template"])
    parser.add_argument("--structure-id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    source_pdb, chains = resolve_structure(args.structure_kind, args.structure_id)
    if not source_pdb.exists():
        raise FileNotFoundError(f"Structure PDB not found: {source_pdb}")

    with tempfile.TemporaryDirectory(prefix="asa_clean_pdb_", dir=".") as tmpdir:
        clean_pdb = Path(tmpdir) / f"{args.structure_id}.pdb"
        extract_selected_chains(source_pdb, chains, clean_pdb)

        naccess = run_naccess(clean_pdb, allowed_chains=chains)
        freesasa_relative = run_freesasa(clean_pdb, allowed_chains=chains)
        rustsasa_relative = run_rustsasa(clean_pdb, allowed_chains=chains)

    payload = {
        "structure_kind": args.structure_kind,
        "structure_id": args.structure_id,
        "source_pdb": str(source_pdb),
        "chains": chains,
        "naccess_relative_asa": naccess,
        "freesasa_relative_asa": freesasa_relative,
        "rustsasa_relative_asa": rustsasa_relative,
        "summary_vs_naccess": {
            "freesasa": summarize_comparison(
                naccess,
                {**{k: {"relative_accessibility_proxy": v} for k, v in freesasa_relative.items()}, "_meta": {}},
            ),
            "rustsasa": summarize_comparison(
                naccess,
                {**{k: {"relative_accessibility_proxy": v} for k, v in rustsasa_relative.items()}, "_meta": {}},
            ),
        },
    }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(output_path)


if __name__ == "__main__":
    main()
