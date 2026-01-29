import os

def three2One(resName):
        RESIDUE_LIST = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q',\
                        'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'}
        if resName in RESIDUE_LIST:
            return RESIDUE_LIST[resName]
        else:
            return 'X'

def dictionaryVDWRadii():
    VDW_RADII = { 'C':1.76, 'N':1.65, 'O':1.40, 'CA':1.87, 'H':1.20, 'S':1.85, 'CB':1.87, 'CZ':1.76, 'NZ':1.50, 'CD':1.81, 'CE':1.81, 'CG':1.81, 'C1':1.80, 'P':1.90 }
    return VDW_RADII

def dictionaryVDWRadiiExtended(resName):
    VDW_RADII_EXTENDED = {  'ALA': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OXT': 1.40  },\
        'ARG': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'NE': 1.65, 'CZ': 1.76, 'NH1': 1.65, 'NH2': 1.65, 'OXT': 1.40 },\
        'ASP': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'OD1': 1.40, 'OD2': 1.40, 'OXT': 1.40  },\
        'ASN': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'OD1': 1.40, 'ND2': 1.65, 'OXT': 1.40 },\
        'CYS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'SG': 1.85, 'OXT': 1.40 },\
        'GLU': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.76, 'OE1': 1.40, 'OE2': 1.40 , 'OXT': 1.40 },\
        'GLN': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.76, 'OE1': 1.40, 'NE2': 1.65, 'OXT': 1.40  },\
        'GLY': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'OXT': 1.40 },\
        'HIS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'ND1': 1.65, 'CD2': 1.76, 'CE1': 1.76, 'NE2': 1.65, 'OXT': 1.40  },\
        'ILE': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG1': 1.87, 'CG2': 1.87, 'CD1': 1.87 , 'OXT': 1.40 },\
        'LEU': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD1': 1.87, 'CD2': 1.87, 'OXT': 1.40 },\
        'LYS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'CE': 1.87, 'NZ': 1.50 , 'OXT': 1.40 },\
        'MET': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'SD': 1.85, 'CE': 1.87, 'OXT': 1.40 },\
        'PHE': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CE1': 1.76, 'CE2': 1.76, 'CZ': 1.76, 'OXT': 1.40  },\
        'PRO': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'OXT': 1.40  },\
        'SER': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OG': 1.40, 'OXT': 1.40 },\
        'THR': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OG1': 1.40, 'CG2': 1.87, 'OXT': 1.40 },\
        'TRP': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'NE1': 1.65, 'CE2': 1.76, 'CE3': 1.76, 'CZ2': 1.76, 'CZ3': 1.76, 'CH2': 1.76, 'OXT': 1.40  },\
        'TYR': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CE1': 1.76, 'CE2': 1.76, 'CZ': 1.76, 'OH': 1.40, 'OXT': 1.40 },\
        'VAL': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG1': 1.87, 'CG2': 1.87, 'OXT': 1.40  }, 'ASX': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'AD1': 1.50, 'AD2': 1.50 }, 'GLX': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD': 1.87, 'AE1': 1.50, 'AE2': 1.50, 'OXT': 1.40  },\
        'ACE': { 'C': 1.76, 'O': 1.40, 'CA': 1.87 }, 'PCA':dictionaryVDWRadii(), 'UNK':dictionaryVDWRadii() }

    return VDW_RADII_EXTENDED[resName]

def createInteractionList(structure,target):
    chainList = []
    alternateLocationDic = {}
    icodeDic = {}
    for line in structure:
        if line[:3] == "END":
            break
        elif line[0:4] == "ATOM":
            chain = line[21]
            atomName = line[12:16].strip()
            resName = line[17:20].strip()
            resSeq = line[22:26].strip()
            coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            altLoc = line[16] #icode and alternate location handled
            icode = line[26]
            key = chain+resSeq
            key2 = key+atomName
            check1 = alternateLocationDic.has_key(key2)
            check2 = icodeDic.has_key(key)
            #these lines are needed to detect alternateLocation or icode differences, the first one is taken...
            if check1 == False and check2 == False:
                alternateLocationDic[key2] = altLoc
                icodeDic[key] = icode
                try:
                    d1 = dictionaryVDWRadiiExtended(resName)
                    if d1.has_key(atomName):
                                d1 = d1[atomName]
                    else:
                        d1 = 0
                    check3 = three2One(resName)
                    check4 = atomName[0]
                    chainList.append([target, chain, coordinates, d1,check3, check4,resName, resSeq])
                except:
                    continue
            elif check1 == False and check2 == True:
                alternateLocationDic[key2] = altLoc
                if icodeDic[key] == icode:
                    try:
                        d1 = dictionaryVDWRadiiExtended(resName)
                        if atomName in d1:
                            d1 = d1[atomName]
                        else:
                            d1 = 0
                        check3 = three2One(resName)
                        check4 = atomName[0]
                        chainList.append([target, chain, coordinates, d1,check3, check4,resName,resSeq])
                    except:
                        continue
    return chainList

def fiberdockInterfaceExtractor(fileName,fileOut,structure1,structure2):
    splited = fileName.split("/")[-1].split("_")
    target1 = splited[1]
    target2 = splited[3]
    chainList1 = createInteractionList(structure1,target1)
    chainList2 = createInteractionList(structure2,target2)
    interact1 = {}
    interact2 = {}
    contact = {}
    for atom1 in chainList1:
        for atom2 in chainList2:
            if atom1[4] != 'X' and atom1[5] != 'H' and atom2[4] != 'X' and atom2[5] != 'H':
                cutoff = atom1[3] + atom2[3] + 0.5
                coordinates1 = atom1[2]
                coordinates2 = atom2[2]
                distance = ((coordinates1[0] - coordinates2[0])**2 + (coordinates1[1] - coordinates2[1])**2 + (coordinates1[2] - coordinates2[2])**2)**0.5
                if distance <= cutoff:
                    key = str(atom1[7])+str(atom2[7])
                    if key not in contact:
                        contact[key] = str(atom1[0])+"_"+str(atom1[1])+"_"+str(atom1[6])+"_"+str(atom1[7])+"\t<-->\t"+str(atom2[0])+"_"+str(atom2[1])+"_"+str(atom2[6])+"_"+str(atom2[7])
                        interact1[atom1[6]+"_"+atom1[7]+"_"+atom1[1]] = 1
                        interact2[atom2[6]+"_"+atom2[7]+"_"+atom2[1]] = 1
                
    ##write contacts             
    filehnd = open(fileOut,"w")
    filehnd.write("Interface Residues Contacts\t\t\n")
    filehnd.write("target1_chain_resName_resNo\t<-->\ttarget2_chain_resName_resNo\n\n")
    for key in contact:
        filehnd.write(contact[key]+"\n")
    filehnd.close()
    os.system("chmod 775 %s" % (fileOut))    
    try:
        ##write flexible refinemnet result in updated way to fix visualization
        filehnd = open(fileName,"w")
        for line in structure1:
            chain = line[21]
            atomName = line[12:16].strip()
            resName = line[17:20].strip()
            resSeq = line[22:26].strip()
            key = resName+"_"+resSeq+"_"+chain
            if key in interact1:
                line = line[:60]+"  3.00"+line[66:]
            else:
                line = line[:60]+"  1.00"+line[66:]
            filehnd.writelines(line)
        filehnd.writelines("TER\n")
        for line in structure2:
            chain = line[21]
            atomName = line[12:16].strip()
            resName = line[17:20].strip()
            resSeq = line[22:26].strip()
            key = resName+"_"+resSeq+"_"+chain
            if key in interact2:
                line = line[:60]+"  2.00"+line[66:]
            else:
                line = line[:60]+"  4.00"+line[66:]
            filehnd.writelines(line)
        filehnd.writelines("END\n")
        filehnd.close()    
        os.system("chmod 775 %s" % (fileName))  
    except Exception as e: 
        print(e)  
