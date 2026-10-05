#!/usr/bin/env python3
"""Stage immutable four-character benchmark PDB inputs before Slurm compute jobs."""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
from pathlib import Path
import sys
import urllib.request

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.pdb_download import normalize_target_id
from Bio.PDB import MMCIFParser, PDBIO


def required_pdb_ids(manifest: Path) -> list[str]:
    ids = set()
    with manifest.open(newline="") as handle:
        for row in csv.DictReader(handle):
            ids.add(normalize_target_id(row["Receptor"])[:4])
            ids.add(normalize_target_id(row["Ligand"])[:4])
    return sorted(ids)


def stage_pdbs(manifest: Path, output_dir: Path, source_roots: list[Path]) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    missing = []
    for pdb_id in required_pdb_ids(manifest):
        destination = output_dir / f"{pdb_id}.pdb"
        if destination.exists() and destination.stat().st_size > 0:
            continue
        candidates = [root / f"{pdb_id}.pdb" for root in source_roots]
        source = next((path for path in candidates if path.exists() and path.stat().st_size > 0), None)
        if source is not None:
            shutil.copy2(source, destination)
            continue
        archive = output_dir / f"{pdb_id}.ent.gz"
        url = f"https://files.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{pdb_id}.ent.gz"
        try:
            urllib.request.urlretrieve(url, archive)
            with gzip.open(archive, "rb") as source_handle, destination.open("wb") as destination_handle:
                shutil.copyfileobj(source_handle, destination_handle)
            archive.unlink()
        except Exception:
            if archive.exists():
                archive.unlink()
            if destination.exists():
                destination.unlink()
            cif_path = output_dir / f"{pdb_id}.cif"
            try:
                urllib.request.urlretrieve(
                    f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", cif_path
                )
                structure = MMCIFParser(QUIET=True).get_structure(pdb_id, cif_path)
                io = PDBIO()
                io.set_structure(structure)
                io.save(destination)
                cif_path.unlink()
            except Exception:
                if cif_path.exists():
                    cif_path.unlink()
                if destination.exists():
                    destination.unlink()
                missing.append(pdb_id)
    return missing


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, action="append", default=[])
    args = parser.parse_args()
    sources = args.source_root or [REPO_ROOT / "benchmark/jobs"]
    # benchmark/jobs is nested by joblist; search it once so callers can use
    # the default without manually enumerating every job directory.
    expanded = []
    for root in sources:
        expanded.extend({path.parent for path in root.rglob("*.pdb")} if root.exists() else [])
    missing = stage_pdbs(args.manifest, args.output_dir, expanded or sources)
    print(f"staged={len(required_pdb_ids(args.manifest)) - len(missing)} missing={len(missing)}")
    if missing:
        print("missing_pdb_ids=" + ",".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
