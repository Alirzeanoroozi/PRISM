import os
import gzip
import shutil
import re
import urllib.request
import pandas as pd

TARGET_DIR = "processed/pdbs"
# Local change: allow test runs to override the pair list without editing inputs.csv.
INPUTS_CSV = os.environ.get("PRISM_INPUTS_CSV", "inputs.csv")
os.makedirs(TARGET_DIR, exist_ok=True)


def split_target_id(target):
    """Split a target id like ``3i6eEF`` into ``('3i6e', ['E', 'F'])``.

    Accepts:
      - 4-char PDB id only (``3i6e``) -> all chains (via probing the file)
      - 4-char PDB id + 1+ chain letters (``3i6eE``, ``3i6eEF``)
    """
    target = str(target).strip()
    if len(target) < 4:
        raise ValueError(f"Invalid target id {target!r}; need at least 4 chars (PDB id)")
    pdb_id = target[:4].lower()
    suffix = target[4:].strip()
    if not suffix:
        return pdb_id, []
    chains = [c for c in suffix if c.isalnum()]
    return pdb_id, chains


def normalize_target_id(target):
    """Return the legacy-compatible ``pdbid + sorted unique chains`` token."""
    value = str(target).strip().replace(" ", "")
    if len(value) < 4:
        raise ValueError(f"target identifier must contain a four-character PDB ID: {target!r}")
    pdb_id = value[:4].lower()
    suffix = re.sub(r"\([^)]*\)", "", value[4:].replace("_", ""))
    chains = "".join(sorted(set(c for c in suffix if c.isalnum())))
    return pdb_id + chains


def target_chain_ids(target):
    """Return the requested chain IDs in canonical order."""
    return tuple(normalize_target_id(target)[4:])


def materialize_target_pdb(target, source_path=None, target_dir=None):
    """Write a chain-filtered PDB containing all chains requested by ``target``."""
    target_dir = target_dir or TARGET_DIR
    canonical = normalize_target_id(target)
    source_path = source_path or os.path.join(target_dir, f"{canonical[:4]}.pdb")
    output_path = os.path.join(target_dir, f"{canonical}.pdb")
    if not target_chain_ids(canonical):
        return source_path
    if os.path.exists(output_path):
        raw_value = str(target).strip()
        legacy_path = os.path.join(target_dir, f"{raw_value}.pdb")
        if len(raw_value) == 5 and legacy_path != output_path and not os.path.exists(legacy_path):
            shutil.copyfile(output_path, legacy_path)
        return output_path

    wanted = set(target_chain_ids(canonical))
    os.makedirs(target_dir, exist_ok=True)
    with open(source_path, "r") as source_handle, open(output_path, "w") as output_handle:
        for line in source_handle:
            if line.startswith(("ATOM", "HETATM", "TER")) and line[21].strip() in wanted:
                output_handle.write(line)
        output_handle.write("END\n")

    # Keep the historical spelling available for callers that used the raw
    # single-chain token as a filename, while canonical consumers use output_path.
    raw_value = str(target).strip()
    legacy_path = os.path.join(target_dir, f"{raw_value}.pdb")
    if len(raw_value) == 5 and legacy_path != output_path and not os.path.exists(legacy_path):
        shutil.copyfile(output_path, legacy_path)
    return output_path

def download_pdb_file(pdb_name, pdb_dir):
    final_pdb = f"{pdb_dir}/{pdb_name}.pdb"
    gz_file = f"{pdb_dir}/{pdb_name}.ent.gz"

    try:
        url = f"https://files.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{pdb_name[:4].lower()}.ent.gz"
        response = urllib.request.urlopen(url)

        with open(gz_file, "wb") as fh:
            fh.write(response.read())

        with gzip.open(gz_file, "rb") as f_in, open(final_pdb, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(gz_file)

        return True
    except Exception as e:
        print(f"PDB download failed for {pdb_name}: {e}")
        return False

def materialize_chain_pdb(target):
    """Backward-compatible wrapper for the old single-chain helper."""
    return materialize_target_pdb(target)

def pdb_downloader(inputs_csv=None):
    """Download targets declared by *inputs_csv* or the configured default."""
    receptor_targets = []
    ligand_targets = []
    raw_targets = []

    # Read pair list
    # Local change: read from env-configurable CSV path.
    df = pd.read_csv(inputs_csv or INPUTS_CSV)
    for _, row in df.iterrows():
        try:
            receptor_targets.append(normalize_target_id(row["Receptor"]))
            ligand_targets.append(normalize_target_id(row["Ligand"]))
            raw_targets.extend((str(row["Receptor"]).strip(), str(row["Ligand"]).strip()))
        except ValueError as exc:
            print(f"Skipping target pair {row['Receptor']}, {row['Ligand']}: {exc}")
            continue
    
    # Download and process PDBs
    for target in list(set(receptor_targets + ligand_targets)):
        if not os.path.exists(f"{TARGET_DIR}/{target[:4].lower()}.pdb"):
            if not download_pdb_file(target[:4].lower(), TARGET_DIR):
                print(f"Failed to download PDB {target}")
                continue
        else:
            print(f"PDB {target} already exists")

    for targeSt in sorted(set(receptor_targets + ligand_targets)):
        if os.path.exists(f"{TARGET_DIR}/{target[:4].lower()}.pdb"):
            materialize_target_pdb(target)
    for raw_target in raw_targets:
        source = f"{TARGET_DIR}/{normalize_target_id(raw_target)[:4]}.pdb"
        if os.path.exists(source):
            materialize_target_pdb(raw_target, source_path=source)
    return list(set(receptor_targets)), list(set(ligand_targets))

if __name__ == "__main__":
    pdb_downloader()
