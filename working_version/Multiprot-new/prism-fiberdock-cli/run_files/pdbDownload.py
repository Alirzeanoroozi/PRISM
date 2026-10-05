#!/usr/bin/env python
#Written by Alper Baspinar
#Updated with PDBj HTTPS + mmCIF(BeEM) fallback by Konuralp Ilim
#CLI/fiberdock port: dropped unused `requests`/`StringIO` hard deps (fetch uses urllib2)
import os
import ConfigParser
import gzip
import shutil
import subprocess, glob
import urllib2


class PDBdownload:
    def __init__(self, workPath):
        self.currentPath = os.getcwd()
        self.workPath = os.path.abspath(workPath)
        os.chdir(self.workPath)

        config = ConfigParser.ConfigParser()
        config.read('prism.ini')
        self.pdbPath = os.path.abspath(config.get('Pdb_Folder', 'pdb_path'))

    def convert_mmcif_to_pdb(self, mmcif_filepath, protein):
        beem_exe = os.path.abspath("../../external_tools/BeEM-master/BeEM")
        if not os.path.isfile(beem_exe):
            raise RuntimeError("BeEM executable not found: " + beem_exe)

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

        # cleanup BeEM intermediates
        for f in glob.glob(os.path.join(pdb_dir, "*-bundle*.pdb")):
            os.remove(f)
        for f in glob.glob(os.path.join(pdb_dir, "*chain-id-mapping*")):
            os.remove(f)

    def download_pdb_file(self, pdb_id):
        """
        Download a structure and return legacy-PDB content (bytes) or None.
        - Try legacy PDB format first (files.pdbj.org)
        - Fall back to mmCIF + BeEM conversion for entries with no legacy PDB
        """
        pdb_id = pdb_id.lower()
        pdb_dir = self.pdbPath

        if not os.path.isdir(pdb_dir):
            os.makedirs(pdb_dir)

        final_pdb = os.path.join(pdb_dir, "{}.pdb".format(pdb_id))
        gz_file = os.path.join(pdb_dir, "pdb{}.ent.gz".format(pdb_id))

        # 1) legacy PDB format
        try:
            url = (
                "https://files.pdbj.org/pub/pdb/data/structures/all/pdb/"
                "pdb{}.ent.gz".format(pdb_id)
            )
            response = urllib2.urlopen(url)

            with open(gz_file, "wb") as fh:
                fh.write(response.read())

            with gzip.open(gz_file, "rb") as f_in, open(final_pdb, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

            os.remove(gz_file)

            with open(final_pdb, "rb") as fh:
                return fh.read()

        except Exception as e:
            print("PDB download failed for {}: {}".format(pdb_id, e))

        # 2) mmCIF fallback (large / PDB-format-unavailable entries)
        try:
            mmcif_gz = os.path.join(pdb_dir, "{}.cif.gz".format(pdb_id))
            mmcif_file = os.path.join(pdb_dir, "{}.cif".format(pdb_id))

            mmcif_url = (
                "https://files.pdbj.org/pub/pdb/data/structures/all/mmCIF/"
                "{}.cif.gz".format(pdb_id)
            )

            response = urllib2.urlopen(mmcif_url)

            with open(mmcif_gz, "wb") as fh:
                fh.write(response.read())

            with gzip.open(mmcif_gz, "rb") as f_in, open(mmcif_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

            os.remove(mmcif_gz)

            self.convert_mmcif_to_pdb(mmcif_file, pdb_id)

            if os.path.exists(mmcif_file):
                os.remove(mmcif_file)

            if not os.path.exists(final_pdb):
                return None

            with open(final_pdb, "rb") as fh:
                return fh.read()

        except Exception as e:
            print("Error downloading {}: {}".format(pdb_id, e))
            return None

    def PDBdownloader(self):
        preLeftTarget = []
        preRightTarget = []
        preTemplateList = []

        # Check for required files
        pairListPath = "lists/pair_list"
        templateListPath = "lists/template_list"
        if not (os.path.exists(pairListPath)) or not (os.path.exists(templateListPath)):
            print "pairlist or templatelist does not exists..."
            os.chdir(self.currentPath)
            return preLeftTarget, preRightTarget, preTemplateList

        # Read pair list
        with open(pairListPath, "r") as filehnd:
            for line in filehnd:
                line = line.strip().split()
                if len(line) == 2:
                    if len(line[0]) >= 4 and len(line[1]) >= 4:
                        preLeftTarget.append(line[0])
                        preRightTarget.append(line[1])

        # Read template list
        with open(templateListPath, "r") as filehnd:
            for line in filehnd:
                preTemplateList.append(line.strip()[:6])

        # Normalise PDB IDs (lowercase id + sorted unique chains)
        leftTarget = []
        for pdb in preLeftTarget:
            newPdb = pdb[:4].lower()
            extra = sorted(set(pdb[4:]))
            leftTarget.append(newPdb + ''.join(extra))

        rightTarget = []
        for pdb in preRightTarget:
            newPdb = pdb[:4].lower()
            extra = sorted(set(pdb[4:]))
            rightTarget.append(newPdb + ''.join(extra))

        if os.path.exists("preprocess"):
            os.chdir(self.currentPath)
            return leftTarget, rightTarget, preTemplateList

        if len(leftTarget + rightTarget) == 0:
            print "Pdblist does not contain any appropriate pdb"
            os.chdir(self.currentPath)
            return leftTarget, rightTarget, preTemplateList

        # Remove duplicate (unordered) pairs
        l, r, tt = [], [], []
        for i in range(len(leftTarget)):
            t = str(leftTarget[i]) + str(rightTarget[i])
            tr = str(rightTarget[i]) + str(leftTarget[i])
            if not (t in tt):
                tt.append(t)
                tt.append(tr)
                l.append(leftTarget[i])
                r.append(rightTarget[i])

        leftTarget, rightTarget = l, r

        # Local per-job pdb directory
        localPdb = os.path.join(self.workPath, "pdb")
        if not os.path.exists(localPdb):
            os.mkdir(localPdb, 0777)

        # Fetch each protein: local job cache -> shared pdb_path -> download
        leftList = []
        rightList = []
        for index in range(len(leftTarget)):
            proteinNameLeft = leftTarget[index]
            proteinNameRight = rightTarget[index]
            check1 = check2 = False

            # left protein
            if not os.path.exists(os.path.join(localPdb, proteinNameLeft[:4] + ".pdb")):
                if os.path.exists(os.path.join(self.pdbPath, proteinNameLeft[:4] + ".pdb")):
                    os.system("cp %s/%s.pdb %s/%s.pdb" % (self.pdbPath, proteinNameLeft[:4], localPdb, proteinNameLeft[:4]))
                    check1 = True
                else:
                    pdb_content = self.download_pdb_file(proteinNameLeft[:4])
                    if pdb_content:
                        with open(os.path.join(localPdb, proteinNameLeft[:4] + ".pdb"), 'wb') as f:
                            f.write(pdb_content)
                        check1 = True
            else:
                check1 = True

            # right protein
            if not os.path.exists(os.path.join(localPdb, proteinNameRight[:4] + ".pdb")):
                if os.path.exists(os.path.join(self.pdbPath, proteinNameRight[:4] + ".pdb")):
                    os.system("cp %s/%s.pdb %s/%s.pdb" % (self.pdbPath, proteinNameRight[:4], localPdb, proteinNameRight[:4]))
                    check2 = True
                else:
                    pdb_content = self.download_pdb_file(proteinNameRight[:4])
                    if pdb_content:
                        with open(os.path.join(localPdb, proteinNameRight[:4] + ".pdb"), 'wb') as f:
                            f.write(pdb_content)
                        check2 = True
            else:
                check2 = True

            if check1 and check2:
                leftList.append(proteinNameLeft)
                rightList.append(proteinNameRight)

        os.chdir(self.currentPath)
        return leftList, rightList, preTemplateList
