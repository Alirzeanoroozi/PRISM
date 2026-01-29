# written by Alper Baspinar
# updated with RCSB API by Konuralp ilim

import os
import subprocess, glob
import os
import gzip
import shutil
import urllib.request

def convert_mmcif_to_pdb(mmcif_filepath, protein):
    beem_exe = os.path.abspath("external_tools/BeEM-master/BeEM")

    mmcif_filepath = os.path.abspath(mmcif_filepath)
    pdb_dir = os.path.dirname(mmcif_filepath)

    ret = subprocess.call([beem_exe, mmcif_filepath], cwd=pdb_dir)
    if ret != 0:
        raise RuntimeError("BeEM failed")

    merge_script = os.path.abspath("../../run_files/merge_bundles.py")
    output_pdb = "{}.pdb".format(protein)

    ret = subprocess.call(
        ["python", merge_script, "*-bundle*.pdb", output_pdb],
        cwd=pdb_dir
    )
    if ret != 0:
        raise RuntimeError("merge_bundles failed")

    # cleanup
    for f in glob.glob(os.path.join(pdb_dir, "*-bundle*.pdb")):
        os.remove(f)

def download_pdb_file(pdb_id):
    pdb_dir = "pdb"
    final_pdb = os.path.join(pdb_dir, "{}.pdb".format(pdb_id))
    gz_file = os.path.join(pdb_dir, "pdb{}.ent.gz".format(pdb_id))

    # --------------------------------------------------
    # 1) Try PDB format
    # --------------------------------------------------
    try:
        url = (
            "https://files.pdbj.org/pub/pdb/data/structures/all/pdb/"
            "pdb{}.ent.gz".format(pdb_id)
        )
        response = urllib.request.urlopen(url)

        with open(gz_file, "wb") as fh:
            fh.write(response.read())

        with gzip.open(gz_file, "rb") as f_in, open(final_pdb, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(gz_file)

        return True
    except Exception as e:
        print("PDB download failed for {}: {}".format(pdb_id, e))
    # --------------------------------------------------
    # 2) Fall back to mmCIF
    # --------------------------------------------------
    try:
        mmcif_gz = os.path.join(pdb_dir, "{}.cif.gz".format(pdb_id))
        mmcif_file = os.path.join(pdb_dir, "{}.cif".format(pdb_id))

        mmcif_url = (
            "https://files.pdbj.org/pub/pdb/data/structures/all/mmCIF/"
            "{}.cif.gz".format(pdb_id)
        )

        response = urllib.request.urlopen(mmcif_url)

        with open(mmcif_gz, "wb") as fh:
            fh.write(response.read())

        with gzip.open(mmcif_gz, "rb") as f_in, open(mmcif_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(mmcif_gz)

        convert_mmcif_to_pdb(mmcif_file, pdb_id)

        if not os.path.exists(final_pdb):
            return None

        with open(final_pdb, "rb") as fh:
            return fh.read()

    except Exception as e:
        print("Error downloading {}: {}".format(pdb_id, e))
        return None

def pdb_downloader():
    targets = []

    # Read pair list
    with open("input_pair_list.txt", "r") as filehnd:
        for line in filehnd:
            line = line.strip().split()
            if len(line) == 2:
                if len(line[0]) >= 4 and len(line[1]) >= 4:
                    targets.append(line[0])
                    targets.append(line[1])

    # Process PDB IDs
    sorted_targets = list(set([pdb[:4].lower() + ''.join(sorted(set(pdb[4:]))) for pdb in targets]))

    # Download and process PDBs
    for target in sorted_targets:
        if not os.path.exists(os.path.join("pdb", target[:4] + ".pdb")):
            download_pdb_file(target[:4])

    return sorted_targets
