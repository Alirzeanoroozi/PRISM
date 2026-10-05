import os
import csv
import shutil

PDB_DIR = "benchmark/processed/pdbs/rigid"               # full 4-letter pdbs are here (e.g. 1fgn.pdb)
CHAINWISE_DIR = "benchmark/processed/pdbs/rigid/chainwise"    # output root
CSV_FILE = "benchmark/data/T_Rigid.csv"

os.makedirs(CHAINWISE_DIR, exist_ok=True)

def parse_pdb_and_chains(token):
    """
    token example: '2FJF_HL' or '4KZN_AB' or '1ABC'
    returns (pdb4_lower, chains_string, token_clean)
      e.g. ('2fjf', 'HL', '2FJF_HL')
    """
    token = (token or "").strip()
    if not token:
        return None, None, None

    if "_" in token:
        pdb_part, chain_part = token.split("_", 1)
    else:
        pdb_part, chain_part = token, ""

    pdb4 = pdb_part[:4].lower()
    chains = "".join([c for c in chain_part.strip() if c.isalnum()])  # keep letters/digits
    return pdb4, chains, token

def write_multi_chain_pdb(in_pdb_path, out_pdb_path, chains_set):
    """
    Keep ATOM/HETATM/ANISOU/TER lines only if chain in chains_set.
    Keep non-coordinate records as-is (HEADER/TITLE/REMARK/etc.).
    """
    with open(in_pdb_path, "r") as fin, open(out_pdb_path, "w") as fout:
        for line in fin:
            rec = line[:6]
            if rec in ("ATOM  ", "HETATM", "ANISOU", "TER   "):
                if len(line) > 21 and line[21] in chains_set:
                    fout.write(line)
            else:
                fout.write(line)

def make_chainwise_from_csv():
    missing_full = []
    written = 0

    with open(CSV_FILE, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            for col in ("PDB ID 1", "PDB ID 2"):
                token = (row.get(col) or "").strip()
                if not token:
                    continue

                pdb4, chains, token_clean = parse_pdb_and_chains(token)
                if not pdb4:
                    continue

                in_pdb = os.path.join(PDB_DIR, f"{pdb4}.pdb")
                if not os.path.exists(in_pdb):
                    missing_full.append((token_clean, in_pdb))
                    continue

                # Folder like: processed/chainwise/2FJF_HL/   (or processed/chainwise/1ABC/)
                out_folder = os.path.join(CHAINWISE_DIR, token_clean)
                os.makedirs(out_folder, exist_ok=True)

                if chains:
                    # One combined file containing all specified chains
                    out_pdb = os.path.join(out_folder, f"{pdb4}_{chains}.pdb")
                    if not os.path.exists(out_pdb):
                        write_multi_chain_pdb(in_pdb, out_pdb, set(chains))
                        written += 1
                else:
                    # No chain specified -> copy whole PDB
                    out_pdb = os.path.join(out_folder, f"{pdb4}.pdb")
                    if not os.path.exists(out_pdb):
                        shutil.copy2(in_pdb, out_pdb)
                        written += 1

    print(f"Chainwise files written/copied: {written}")

    if missing_full:
        print("\n[WARN] Missing full PDB files (not found in processed/pdbs):")
        for token, path in missing_full[:20]:
            print(f"  {token} -> expected {path}")
        if len(missing_full) > 20:
            print(f"  ... and {len(missing_full)-20} more")

if __name__ == "__main__":
    make_chainwise_from_csv()
