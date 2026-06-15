import os
import gzip
import shutil
import urllib.request
import pandas as pd

TARGET_DIR = "processed/pdbs"
os.makedirs(TARGET_DIR, exist_ok=True)

PDB_URLS = (
    "https://files.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{pdb_id}.ent.gz",
    "https://files.rcsb.org/download/{pdb_id_upper}.pdb.gz",
)


def split_target_id(target: str):
    """Split a target id like `3i6eEF` into (`3i6e`, [`E`,`F`]).

    Accepts:
      - 4-char PDB id only (`3i6e`) -> all chains
      - 4-char PDB id + 1+ chain letters (`3i6eE`, `3i6eEF`, `3i6eABC`)
    """
    target = target.strip()
    if len(target) < 4:
        raise ValueError(f"Invalid target id {target!r}; need at least 4 chars (PDB id)")
    pdb_id = target[:4].lower()
    chains = list(target[4:])
    return pdb_id, chains


def download_pdb_file(pdb_id, pdb_dir):
    pdb_id = pdb_id.lower()
    final_pdb = f"{pdb_dir}/{pdb_id}.pdb"
    gz_file = f"{pdb_dir}/{pdb_id}.ent.gz"

    if os.path.exists(final_pdb):
        return True

    last_exc = None
    for url_template in PDB_URLS:
        url = url_template.format(pdb_id=pdb_id, pdb_id_upper=pdb_id.upper())
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                with open(gz_file, "wb") as fh:
                    fh.write(response.read())
            with gzip.open(gz_file, "rb") as f_in, open(final_pdb, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_file)
            return True
        except Exception as exc:
            last_exc = exc
            if os.path.exists(gz_file):
                try:
                    os.remove(gz_file)
                except OSError:
                    pass
    print(f"PDB download failed for {pdb_id}: {last_exc}")
    return False


def pdb_downloader(args):
    receptor_targets = []
    ligand_targets = []

    df = pd.read_csv(args.inputs_csv)
    df.columns = [c.strip() for c in df.columns]
    if "Receptor" not in df.columns or "Ligand" not in df.columns:
        raise ValueError(
            f"Input CSV {args.inputs_csv} must have 'Receptor' and 'Ligand' columns (found {list(df.columns)})"
        )

    for receptor, ligand in zip(df["Receptor"], df["Ligand"]):
        try:
            receptor_pdb_id, _ = split_target_id(receptor)
            ligand_pdb_id, _ = split_target_id(ligand)
        except ValueError as exc:
            print(f"Skipping invalid row ({receptor!r}, {ligand!r}): {exc}")
            continue

        if not download_pdb_file(receptor_pdb_id, TARGET_DIR):
            print(f"Failed to download PDB {receptor_pdb_id}; skipping pair")
            continue
        if not download_pdb_file(ligand_pdb_id, TARGET_DIR):
            print(f"Failed to download PDB {ligand_pdb_id}; skipping pair")
            continue

        receptor_targets.append(receptor.strip())
        ligand_targets.append(ligand.strip())

    assert len(receptor_targets) == len(ligand_targets)
    return receptor_targets, ligand_targets
