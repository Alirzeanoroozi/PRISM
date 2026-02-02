import os

from naccess_utils import get_asa_complex

RSATHRESHOLD = 15.0
SCFFTHRESHOLD = 5.0

TARGET_RSA_DIR = "templates/rsas"
os.makedirs(TARGET_RSA_DIR, exist_ok=True)

def extract_surface(queries):
    for protein in queries:
        if run_naccess(protein) == 1:
            asaDict = extract_surface(protein)
            if len(asaDict) != 0:
                with open(f"surface_extraction/{protein}.asa.pdb", "w") as asahnd:
                    for key in sorted(asaDict.keys()):
                        asahnd.writelines(asaDict[key])
                    asahnd.writelines("END")

def extract_surface(protein):
    asa_complex = get_asa_complex(protein, TARGET_RSA_DIR)
    rsa_residues = [key for key in asa_complex if asa_complex[key] > RSATHRESHOLD]

    for asaline in asahandler.readlines():
        if asaline[:3] == "END":
            break
        if asaline[:4] == "ATOM":
            atom = asaline[13:15]
            resNo = asaline[22:26].strip() 
            chain = asaline[21]
            serialNumber = asaline[6:11]
            for reschn, resseq in rsalist:
                if atom == "CA" and resNo == resseq and chain == reschn:
                    asaDict[serialNumber] = asaline

    for line in asahandler.readlines():
        if line[:3] == "END":
            break
        serialNumber = line[6:11]
        if line[:4] == "ATOM" and line[13:15] == "CA" and serialNumber not in asaDict:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            for key in asaDict:
                xcoor = float(asaDict[key][30:38])
                ycoor = float(asaDict[key][38:46])
                zcoor = float(asaDict[key][46:54])
                dist = ((xcoor-x)**2 + (ycoor-y)**2 + (zcoor-z)**2)**0.5
                if dist <= SCFFTHRESHOLD:
                    if serialNumber not in asaDict:
                        asaDict[serialNumber] = line
    return asaDict
