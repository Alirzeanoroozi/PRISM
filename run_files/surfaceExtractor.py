import os

RSATHRESHOLD = 15.0
SCFFTHRESHOLD = 5.0

os.makedirs("surface_extraction", exist_ok=True)

def extract_surface(queries):
    for protein in queries:
        if run_naccess(protein) == 1:
            asaDict = extract_surface(protein)
            output_asa_file(protein, asaDict)

def run_naccess(protein):
    protein_path = f"pdb/{protein}.pdb"

    if os.path.exists(protein_path):
        os.system(f"external_tools/naccess/naccess {protein_path}")
    else:
        print(f"protein {protein} does not exist in the protein path, check previous stages ")
        return -1
    
    rsa_file = protein + ".rsa"
    asa_file = protein + ".asa"
    logfile = protein + ".log"
    naccess_output_file = f"surface_extraction/{protein}.rsa"
    
    if os.path.exists(rsa_file):
        os.system("mv %s %s"%(rsa_file, naccess_output_file))
    if os.path.exists(asa_file):
        os.system("rm %s" % (asa_file))
    if os.path.exists(logfile):
        os.system("rm %s" % (logfile))
    return 1

def extract_surface(protein):
    rsaPath = f"surface_extraction/{protein}.rsa"
    proteinPath = f"pdb/{protein}.pdb" 
    rsalist = []
    asaDict = {}

    with open(rsaPath, "r") as rsahandler, open(proteinPath, "r") as asahandler:
        for rsaline in rsahandler.readlines():
            if rsaline[:3] == "RES":
                rsaline = rsaline.strip()
                resR = rsaline[4:7]
                standard = standard_data(resR)
                if standard != -1:
                    reschn = rsaline[8]
                    resseq = rsaline[9:13].strip()
                    try:
                        absoluteacc = float(rsaline[14:22])
                    except:
                        print(f"naccess output handle error for protein {protein}")
                        return asaDict
                    relativeacc = absoluteacc * 100 / standard
                    if relativeacc > RSATHRESHOLD:
                        rsalist.append((reschn,resseq))

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

    return find_scaffold(asaDict, proteinPath)

def find_scaffold(asaDict, proteinPath):
    with open(proteinPath, "r") as pdbhnd:
        for line in pdbhnd.readlines():
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

def standard_data(resname):
    return {
        'ALA':107.95,
        'CYS':134.28,
        'ASP':140.39,
        'GLU':172.25,
        'PHE':199.48,
        'GLY':80.1,
        'HIS':182.88,
        'ILE':175.12,
        'LYS':200.81,
        'LEU':178.63,
        'MET':194.15,
        'ASN':143.94,
        'PRO':136.13,
        'GLN':178.5,
        'ARG':238.76,
        'SER':116.5,
        'THR':139.27,
        'VAL':151.44,
        'TRP':249.36,
        'TYR':212.76,
    }.get(resname, -1)

def output_asa_file(protein, asaDict):
    if len(asaDict) != 0:
        with open(f"surface_extraction/{protein}.asa.pdb", "w") as asahnd:
            for key in sorted(asaDict.keys()):
                asahnd.writelines(asaDict[key])
            asahnd.writelines("END")
            asahnd.close()

