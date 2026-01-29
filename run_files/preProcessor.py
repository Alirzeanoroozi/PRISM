# a class to do preprocessing on the proteins before surface extraction
# written by Alper Baspinar

import os

amino = [
    'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
    'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'
]

def preprocess_input_proteins(targets):    
    for target in targets:
        path = "pdbs/%s.pdb" % target[0:4]
        extra = target[4:] #if pdb also contains chains pdb files should be splitted
        chainWriter(target, splitChains(extra, path) if len(extra) != 0 else splitNotChain(path))
        
def splitNotChain(path):
    residueList = []
    alternateLocationDic = {}
    icodeDic = {}
    with open(path, "r") as filehnd:
        for line in filehnd.readlines():
            if line[:3] == "END":
                break
            if line[0:4] == "ATOM":
                chain = line[21]
                atomName = line[12:16]
                resName = line[17:20]
                residueSeq = line[22:26]
                altLoc = line[16] #icode and alternate location handled
                icode = line[26]
                key = chain+residueSeq
                key2 = key+atomName
                check1 = alternateLocationDic.has_key(key2)
                check2 = icodeDic.has_key(key)
                check3 = resName in amino
                #these lines are needed to detect alternateLocation or icode differences, the first one is taken...
                if check1 == False and check2 == False and check3 == True:
                    alternateLocationDic[key2] = altLoc
                    icodeDic[key] = icode
                    residueList.append(line)
                elif check1 == False and check2 == True and check3 == True:
                    alternateLocationDic[key2] = altLoc
                    if icodeDic[key] == icode:
                        residueList.append(line)
    return residueList
    
def splitChains(extra, path):
    residueList = []
    alternateLocationDic = {}
    icodeDic = {}
    with open(path, "r") as filehnd:
        for line in filehnd.readlines():
            if line[:3] == "END":
                break
            if line[0:4] == "ATOM":
                chain = line[21]
                if chain in extra:
                    atomName = line[12:16]
                    residueSeq = line[22:26]
                    resName = line[17:20]
                    altLoc = line[16]
                    icode = line[26]
                    key = chain+residueSeq
                    key2 = key+atomName
                    check1 = key2 in alternateLocationDic
                    check2 = key in icodeDic
                    check3 = resName in amino
                    #these lines are needed to detect alternateLocation or icode differences, the first one is taken...
                    if check1 == False and check2 == False and check3 == True:
                        alternateLocationDic[key2] = altLoc
                        icodeDic[key] = icode
                        residueList.append(line)
                    elif check1 == False and check2 == True and check3 == True:
                        alternateLocationDic[key2] = altLoc
                        if icodeDic[key] == icode:
                            residueList.append(line)

    return residueList
    
def chainWriter(protein, residueList):
    os.makedirs("preprocess", exist_ok=True)
    writePath = "preprocess/"+protein+".pdb"
    with open(writePath, "w") as filehnd:
        for residue in residueList:
            filehnd.writelines(residue)
        filehnd.writelines("END")
