import os
import socket
import csv
import gzip
import shutil
import urllib.request

#replace TARGET_DIR and read .csb file

TARGET_DIR = "benchmark/processed/pdbs/midium"
os.makedirs(TARGET_DIR, exist_ok=True)

def download_pdb_file(pdb_name, pdb_dir):
    pdb4 = pdb_name[:4].lower()
    final_pdb = f"{pdb_dir}/{pdb4}.pdb"
    gz_file = f"{pdb_dir}/{pdb4}.ent.gz"

    url = f"https://files.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{pdb4}.ent.gz"

    try:
        print(f"[DOWNLOAD] {pdb4}  ->  {url}")
        response = urllib.request.urlopen(url, timeout=30)

        with open(gz_file, "wb") as fh:
            fh.write(response.read())

        # Decompress
        with gzip.open(gz_file, "rb") as f_in, open(final_pdb, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(gz_file)

        # Sanity check: empty file means something went wrong
        if os.path.getsize(final_pdb) == 0:
            raise RuntimeError("Downloaded PDB is empty after decompression.")

        return True

    except urllib.error.HTTPError as e:
        # Very common: 404 means that entry isn't available at this mirror
        print(f"[HTTP ERROR] {pdb4} status={e.code} reason={e.reason} url={url}")
        return False

    except urllib.error.URLError as e:
        print(f"[URL ERROR] {pdb4} reason={e.reason} url={url}")
        return False

    except (socket.timeout, TimeoutError) as e:
        print(f"[TIMEOUT] {pdb4} error={e} url={url}")
        return False

    except OSError as e:
        # Includes gzip errors like "Not a gzipped file"
        print(f"[OS/GZIP ERROR] {pdb4} error={e} url={url}")
        return False

    except Exception as e:
        print(f"[OTHER ERROR] {pdb4} error={e} url={url}")
        return False

def pdb_downloader():
    targets = []

    # Read CSV and collect PDB IDs from columns: "PDB ID 1" and "PDB ID 2"
    with open("benchmark/data/T_medium.csv", "r", newline="") as filehnd:
        reader = csv.DictReader(filehnd)
        for row in reader:
            pdb1 = (row.get("PDB ID 1") or "").strip()
            pdb2 = (row.get("PDB ID 2") or "").strip()

            # Expect 4-letter PDB codes in these columns (e.g., 2FJF, 4KZN)
            if len(pdb1) >= 4:
                targets.append(pdb1[:4])
            else:
                if pdb1:
                    print(f"Skipping PDB ID 1 value '{pdb1}' (too short).")

            if len(pdb2) >= 4:
                targets.append(pdb2[:4])
            else:
                if pdb2:
                    print(f"Skipping PDB ID 2 value '{pdb2}' (too short).")

    failed = []
    # Download and process PDBs
    for target in list(set(targets)):
        target4 = target[:4].lower()
        outpath = f"{TARGET_DIR}/{target4}.pdb"

        if not os.path.exists(outpath):
            ok = download_pdb_file(target4, TARGET_DIR)
            if not ok:
                print(f"[FAIL] {target4}")
                failed.append(target4)
                continue
        else:
            print(f"[SKIP] already exists: {target4}")

    # write failures
    if failed:
        with open("failed_downloads.txt", "w") as f:
            for x in sorted(set(failed)):
                f.write(x + "\n")
        print(f"Saved failures to failed_downloads.txt ({len(set(failed))} unique).")
    return list(set(targets))

if __name__ == "__main__":
    pdbs = pdb_downloader()
    print(f"Done. Unique PDB IDs collected: {len(pdbs)}")