# written by Alper Baspinar
# this class gets the interface list from user then checks if it exists or not. If it doesn't exist, it creates interface pdb and hotspot_data for interfaces.

import os

def template_generator():
    calculated_templates = []
    for template in [line.strip()[:6] for line in open("templates.txt", "r").readlines()]:
        try:
            print(f"Interface Generation for {template} started!!")
            create_interface(template)
            print(f"Interface Generation for {template} Finished!!")
            print(f"HotSpot Generation for {template} Started!!")
            hotspot_creator(template)
            print(f"HotSpot Generation for {template} Finished!!")
            print(f"Contact Generation for {template} Started!!")
            contact_writer(template)
            print(f"Contact Generation for {template} Finished!!")
            calculated_templates.append(template)
        except Exception as e:
            print(f"Error in template generation for {template}: {e}")
            continue
    return calculated_templates

def distance_calculator(coordinates1, coordinates2):
    return ((float(coordinates1[0]) - float(coordinates2[0]))**2 + (float(coordinates1[1]) - float(coordinates2[1]))**2 + (float(coordinates1[2]) - float(coordinates2[2]))**2)**0.5

def create_interface(template):
    protein, chain_id1, chain_id2 = template[:4].lower(), template[4], template[5]
    pdb_path = f"pdb/{protein}.pdb"
    left = set()
    right = set()
    chain_list1 = []
    chain_list2 = []
    chain_interact1_res = {}
    chain_interact2_res = {}
    chain_ca1 = {}
    chain_ca2 = {}

    if os.path.exists(pdb_path):
        os.system(f"external_tools/naccess/naccess {pdb_path}")
    else:
        print(f"template {template} does not exist")
    ######################################################################
    for line in lineList:
        if line[:3] == "END":
            break
        elif line[:4] == "ATOM":
            atomName = line[12:16].strip()
            resName = line[17:20]
            chainId = line[21]
            resSeq = int(line[22:26].strip())
            coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if chainId == chainId1:
                d1 = self.dictionaryVDWRadiiExtended(resName)
                try:
                    d1 = d1[atomName]
                except:
                    d1 = 0
                check1 = self.three2One(resName)
                check2 = atomName[0]
                leftValue = self.three2One(resName)+str(resSeq)+chainId1
                leftContact = chainId1+"."+self.three2One(resName)+"."+str(resSeq)+"\t"
                chainList1.append([line, resSeq, coordinates, d1,check1, check2, leftValue,leftContact])
                if atomName == 'CA':
                    chainCa1[resSeq] = [line,coordinates]
            elif chainId == chainId2:
                d1 = self.dictionaryVDWRadiiExtended(resName)
                try:
                    d1 = d1[atomName]
                except:
                    d1 = 0
                check1 = self.three2One(resName)
                check2 = atomName[0]
                rightValue = self.three2One(resName)+str(resSeq)+chainId2
                rightContact = chainId2+"."+self.three2One(resName)+"."+str(resSeq)
                chainList2.append([line, resSeq, coordinates, d1, check1, check2, rightValue,rightContact])
                if atomName == 'CA':
                    chainCa2[resSeq] = [line,coordinates]

    for atom1 in chainList1:
        for atom2 in chainList2:
            if atom1[4] != 'X' and atom1[5] != 'H' and atom2[4] != 'X' and atom2[5] != 'H':
                cutoff = atom1[3] + atom2[3] + 0.5
                distance = self.distanceCalculator(atom1[2],atom2[2])
                if distance <= cutoff:
                    self.left.add(atom1[6])
                    self.right.add(atom2[6])
                    key = str(atom1[1])+str(atom2[1])
                    if self.contact.has_key(key) == False:
                        self.contact[key] = atom1[7]+atom2[7]
                    try:
                        chainInteract1Res[atom1[1]] = chainCa1[atom1[1]]
                    except KeyError:
                        continue
                    try:
                        chainInteract2Res[atom2[1]] = chainCa2[atom2[1]]
                    except KeyError:
                        continue
    self.left = list(self.left)
    self.left.sort()
    self.right = list(self.right)
    self.right.sort()
    #adding nearby residues
    #first chain
    atom1Iterable = chainCa1.keys()
    atom2Iterable = chainInteract1Res.keys()
    for key1 in atom1Iterable:
        if key1 not in atom2Iterable:
                for key2 in atom2Iterable:
                    distance = self.distanceCalculator(chainCa1[key1][1],chainInteract1Res[key2][1])
                    if distance <= self.cutOff:
                        try:
                            chainInteract1Res[key1] = chainCa1[key1]
                        except KeyError:
                            continue

    #second chain
    atom1Iterable = chainCa2.keys()
    atom2Iterable = chainInteract2Res.keys()
    for key1 in atom1Iterable:
        if key1 not in atom2Iterable:
                for key2 in atom2Iterable:
                    distance = self.distanceCalculator(chainCa2[key1][1],chainInteract2Res[key2][1])
                    if distance <= self.cutOff:
                        try:
                            chainInteract2Res[key1] = chainCa2[key1]
                        except KeyError:
                            continue

    if len(chain_dict1) == 0 or len(chain_dict2) == 0:
        print(f"Interface {interface} does not exist!!")
        return 0
    with open(f"template/interfaces/{interface}_{interface[4]}.int","w") as interface_writer_left:
        for key in sorted(chain_dict1.iterkeys()):
            interface_writer_left.writelines(chain_dict1[key][0])
    with open(f"template/interfaces/{interface}_{interface[5]}.int","w") as interface_writer_right:
        for key in sorted(chain_dict2.iterkeys()):
            interface_writer_right.writelines(chain_dict2[key][0])
    
def three2one(res_name):
    return {
        'ALA':'A',
        'CYS':'C',
        'ASP':'D',
        'GLU':'E',
        'PHE':'F',
        'GLY':'G',
        'HIS':'H',
        'ILE':'I',
        'LYS':'K',
        'LEU':'L',
        'MET':'M',
        'ASN':'N',
        'PRO':'P',
        'GLN':'Q',
        'ARG':'R',
        'SER':'S',
        'THR':'T',
        'VAL':'V',
        'TRP':'W',
        'TYR':'Y'
    }.get(res_name, 'X')

def vdw_radii():
    return { 'C':1.76, 'N':1.65, 'O':1.40, 'CA':1.87, 'H':1.20, 'S':1.85, 'CB':1.87, 'CZ':1.76, 'NZ':1.50, 'CD':1.81, 'CE':1.81, 'CG':1.81, 'C1':1.80, 'P':1.90 }

def vdw_radii_extended(res_name):
    return {  
        'ALA': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OXT': 1.40  },\
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
        'ACE': { 'C': 1.76, 'O': 1.40, 'CA': 1.87 },
        'PCA':dictionary_vdw_radii(),
        'UNK':dictionary_vdw_radii()
    }

def contact_writer(interface):
    contactWrite = open(f"template/contact/{interface}.txt","w")
    for key in contact.keys():
        contactWrite.writelines(contact[key]+"\n")
    contact = {}

def contact_potentials(pdb_name, chain1,chain2):
    contactDict = self.contactingResidues(pdbName,chain1,chain2)
    ddPP = {}
    allResidues = self.Ca_CoordinatesFetch_BothChain(pdbName,chain1,chain2)
    matrixPotentials = self.PairPot()
    for residues in contactDict.keys():
        totalContact = 0
        res1 = residues
        for res2 in contactDict[res1]:
            aa1 = res1[4]
            aa2 = res2[4]
            if aa1 != 'X' and aa2 != 'X':
                temp = [aa1,aa2]
                temp.sort()
                contact = temp[0]+"-"+temp[1]
                totalContact += matrixPotentials[contact]
        if len(contactDict[res1]) != 0:
            ddPP[res1] = totalContact
        if len(contactDict[res1]) == 0:
            ddPP[res1] = 0.0
    for item in allResidues.keys():
        try:
            ddPP[pdbName+self.three2One(item[0:3])+item[3:0]]
        except KeyError:
            ddPP[pdbName+self.three2One(item[0:3])+item[3:0]] = 0.0
    return ddPP

def PairPot(self):
    PairPotentialDict = {'Q-Y': -4.75, 'Q-Q': -4.71, 'Q-S': -3.06, 'Q-R': -4.13, 'H-R': -3.24, 'Q-T': -3.44, 'Q-W': -4.87, 'Q-V': -4.29, 'C-F': -4.81, 'S-W': -2.86, 'E-E': -2.85, 'C-N': -1.76, 'E-G': -1.58, 'C-L': -6.13, 'C-K': -2.28, 'C-I': -5.53, 'C-H': -4.79, 'C-G': -3.43, 'E-L': -4.35, 'C-E': -3.02, 'C-D': -3.56, 'E-I': -3.99, 'E-H': -3.15, 'E-K': -2.83, 'E-T': -2.34, 'E-W': -3.84, 'E-V': -3.2, 'E-Q': -3.45, 'E-P': -1.7, 'C-Y': -4.71, 'E-R': -4.05, 'C-W': -4.7, 'C-V': -4.82, 'C-T': -2.92, 'C-S': -3.41, 'C-R': -4.41, 'C-Q': -4.37, 'C-P': -2.97, 'M-Y': -5.01, 'H-Y': -3.78, 'H-T': -2.73, 'M-T': -3.33, 'H-V': -4.14, 'M-V': -5.16, 'M-Q': -4.46, 'M-P': -3.11, 'M-S': -3.33, 'M-R': -4.13, 'H-L': -4.85, 'H-M': -4.47, 'H-N': -3.05, 'M-N': -3.05, 'H-H': -3.08, 'H-I': -4.55, 'H-K': -2.14, 'S-Y': -3.43, 'F-Y': -5.33, 'F-V': -5.33, 'S-V': -2.95, 'K-M': -2.36, 'K-L': -3.24, 'F-R': -4.51, 'F-S': -4.01, 'F-P': -3.46, 'F-Q': -4.89, 'F-N': -3.89, 'K-V': -2.19, 'F-L': -6.65, 'F-M': -5.58, 'K-S': -1.91, 'F-K': -2.7, 'F-H': -3.82, 'F-I': -6.11, 'F-F': -6.45, 'F-G': -3.76, 'F-T': -3.85, 'K-Y': -2.64, 'S-T': -2.22, 'N-W': -3.31, 'S-S': -2.14, 'V-Y': -4.54, 'V-V': -4.86, 'V-W': -4.7, 'K-W': -3.12, 'K-T': -1.64, 'N-Y': -3.14, 'K-N': -1.63, 'K-R': -2.11, 'A-I': -4.45, 'A-H': -3.29, 'A-K': -1.91, 'A-M': -4.42, 'A-L': -5.02, 'A-N': -2.06, 'A-A': -3.51, 'A-C': -3.99, 'K-P': -0.5, 'A-E': -2.69, 'A-D': -2.74, 'A-G': -2.24, 'A-F': -4.43, 'A-Y': -3.96, 'E-F': -4.1, 'A-Q': -3.65, 'A-P': -1.78, 'A-S': -2.49, 'A-R': -3.26, 'A-T': -2.38, 'A-W': -4.66, 'A-V': -3.86, 'G-K': -1.32, 'G-I': -3.5, 'G-H': -2.47, 'G-N': -1.63, 'G-M': -3.02, 'G-L': -3.79, 'I-I': -5.97, 'I-K': -2.88, 'I-M': -5.37, 'I-L': -6.67, 'I-N': -3.42, 'I-Q': -4.65, 'I-P': -3.09, 'I-S': -3.47, 'I-R': -4.29, 'I-T': -3.61, 'I-W': -5.79, 'I-V': -5.31, 'I-Y': -5.39, 'G-R': -2.71, 'G-Q': -3.02, 'G-P': -1.01, 'G-W': -3.77, 'G-V': -2.77, 'G-T': -1.81, 'R-R': -3.98, 'C-M': -5.07, 'E-N': -2.43, 'C-C': -7.23, 'K-K': -1.02, 'P-P': -0.83, 'P-Q': -2.67, 'K-Q': -2.56, 'F-W': -5.72, 'D-M': -3.69, 'E-S': -2.57, 'W-Y': -5.12, 'Y-Y': -4.58, 'E-Y': -3.62, 'W-W': -5.16, 'P-T': -1.24, 'R-S': -2.78, 'P-V': -2.43, 'P-W': -3.96, 'R-V': -3.7, 'R-W': -4.54, 'R-T': -2.43, 'P-S': -1.56, 'G-G': -1.77, 'R-Y': -4.32, 'P-Y': -2.92, 'D-D': -2.47, 'D-E': -2.35, 'D-F': -3.52, 'D-G': -1.84, 'D-H': -2.93, 'D-I': -3.46, 'D-K': -2.47, 'D-L': -4.16, 'G-Y': -3.34, 'D-N': -2.5, 'D-P': -1.24, 'D-Q': -3.15, 'D-R': -3.92, 'D-S': -2.34, 'D-T': -2.41, 'D-V': -3.28, 'D-W': -3.81, 'D-Y': -3.63, 'M-W': -5.49, 'N-V': -2.76, 'L-Y': -5.67, 'N-T': -2.05, 'H-W': -4.5, 'N-R': -2.67, 'N-S': -2.13, 'N-P': -1.01, 'N-Q': -3.0, 'L-P': -3.75, 'L-Q': -5.36, 'L-R': -5.01, 'L-S': -4.12, 'L-T': -4.28, 'L-V': -6.03, 'L-W': -6.4, 'H-P': -1.8, 'G-S': -1.68, 'L-L': -7.16, 'L-M': -6.24, 'L-N': -3.94, 'E-M': -3.76, 'N-N': -1.99, 'H-S': -3.12, 'M-M': -5.89, 'T-T': -2.12, 'T-V': -2.87, 'T-W': -3.48, 'T-Y': -3.33, 'H-Q': -4.2, 'P-R': -2.43}
    return PairPotentialDict

def contactingResidues(self,pdbName,chain1, chain2):
    """
    ###########Contacting residues according to center of mass#########
    This function extracts contacting residues from each from one chain. The cutoff is defined as 7 Angstrom.
    Contacting pairs of each residue is dumped in a dictionary of each individual residue.
    This dictionary will be used to calculate delta pair potentials for each residue.
    """
    centeredResidues = self.center_of_mass(pdbName,chain1,chain2)
    contactDict = {}
    temp = []
    for res1 in centeredResidues.keys():
        coor1 = centeredResidues[res1]
        for res2 in centeredResidues.keys():
            coor2 = centeredResidues[res2]
            if res1 != res2:
                if res1[-1] != res2[-1]:
                    distance = self.distanceCalculator(coor1,coor2)
                    if distance <= 7.0:
                        try:
                            if contactDict.has_key(res1):
                                temp = contactDict[res1]
                            temp.append(res2)
                            contactDict.update({res1: temp})
                            temp = []
                        except KeyError:
                            continue
                elif res1[-1] == res2[-1]:
                    if abs(int(res1[5:len(res1)-1])-int(res2[5:len(res2)-1])) > 3:
                        distance = self.distanceCalculator(coor1,coor2)
                        if distance <= 7.0:
                            try:
                                if contactDict.has_key(res1):
                                    temp = contactDict[res1]
                                temp.append(res2)
                                contactDict.update({res1: temp})
                                temp = []
                            except KeyError:
                                continue
    return contactDict
##############################################Center of Mass###############################################################
## Specify the residue atoms of interest to extract the center of mass coordinates. Dump these coordinates into a dictionary,
## where keys are residue ID(residue name + residue number + chain), the values x,y,z coordinates.
## Use this funtion in "contactingResidues" function to extract contacting pairs.
def center_of_mass(self,pdbName,chain1,chain2):
    alphaCarbons = self.Ca_CoordinatesFetch_BothChain(pdbName,chain1,chain2)

    ##dictionary of interested atoms which are sent by okeskin.
    dictAtom = {'ACE':['CA'], 'ALA':['CB'], 'GLY':['CA'], 'SER':['OG'], 'ASN':['OD1', 'ND2'], 'ASP':['OD1', 'OD2'], 'CYS':['SG'], 'PRO':['CB', 'CG', 'CD'], \
                'THR':['OG1'], 'VAL':['CG1', 'CG2'], 'ARG':['NE', 'NH1', 'NH2'], 'GLN':['OE1', 'NE2'], 'GLU':['OE1', 'OE2'], 'HIS':['CG', 'ND1', 'CD2', 'CE1', 'NE2'], 'ILE':['CD1'], \
                'LEU':['CD1', 'CD2'], 'LYS':['NZ'], 'MET':['SD'], 'PHE':['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'TRP':['CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'], 'TYR':['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']}

    centerCoordinates = {}
    for residues in alphaCarbons.keys():
        resName = residues[0:3]
        allAtoms = alphaCarbons[residues]
        x_total = 0
        y_total = 0
        z_total = 0
        coordinateDict = {}
        for ind_atom in allAtoms:
            coordinateDict[ind_atom[0]] = ind_atom[1]
        atoms = dictAtom[resName]
        try:
            for item in atoms:
                if coordinateDict.has_key(item):
                    coordinates = coordinateDict[item]
                else:
                    coordinates = coordinateDict['CA']
                x = coordinates[0]
                y = coordinates[1]
                z = coordinates[2]
                x_total += x
                y_total += y
                z_total += z
            real_coordinates = [x_total/len(atoms), y_total/len(atoms), z_total/len(atoms)]
            centerCoordinates[pdbName+self.three2One(resName)+residues[3:]] = real_coordinates
        except KeyError:
            continue
    return centerCoordinates
##############################################Atom Fetch###################################################################
## This function extracts all atoms and their coordinates for given pdb name and chain names.
## It returns the dictionary of atoms for individual residues in the protein chains.
## This function is used by "center_of_mass"
def Ca_CoordinatesFetch_BothChain(self,pdbName,chain1,chain2):
    alphaCarbons = {}
    atomList = []
    filehnd = open(self.pdbPath+"/%s.pdb" % pdbName, 'r')
    for line in filehnd.readlines():
        if line[:3] == "END":
            break
        elif line[:4] == "ATOM" and line[26] == " ":
            atomName = line[12:16].replace(" ","")
            resName = line[17:20]
            chainId = line[21]
            resSeq = int(line[22:26].replace(" ",""))
            coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if chainId == chain1 or chainId == chain2:
                try:
                    key = resName+str(resSeq)+chainId
                    if alphaCarbons.has_key(key):
                        atomList = alphaCarbons[key]
                    atomList.append([atomName,coordinates])
                    alphaCarbons.update({key: atomList})
                    atomList = []
                except KeyError:
                    continue
    filehnd.close()
    return alphaCarbons

##############################################ASA FILE READING#############################################################

def asaComplex(self,interface):
    #ASA for complex
    ASAcomplex = {}
    filehnd = open("%s.rsa" % interface, 'r')
    if self.external_tool_choice == 0:
        for rsaline in filehnd.readlines():
            standard = self.StandardData(rsaline[:3])
            if standard != -1:
                rsatemp = rsaline.strip().split()
                resName = self.three2One(rsaline[:3])
                chain = rsatemp[1]
                resPosition = rsatemp[2]
                try:
                    absoluteacc = float(rsatemp[5])
                except:
                    print("pops output handle error for protein %s" % interface)
                    return ASAcomplex
                relativeacc = absoluteacc*100/standard
                ASAcomplex[resName+resPosition+chain] = relativeacc
    else:
        for rsaline in filehnd.readlines():
            if rsaline[:3] == "RES":
                rsaline = rsaline.strip()
                resR = rsaline[4:7]
                standard = self.StandardData(resR)
                if standard != -1:
                    resName = self.three2One(resR)
                    chain = rsaline[8]
                    resPosition = rsaline[9:13].strip()
                    try:
                        absoluteacc = float(rsaline[14:22])
                    except:
                        print("naccess output handle error for protein %s" % proteinname)
                        return ASAcomplex
                    relativeacc = absoluteacc*100/standard
                    ASAcomplex[resName+resPosition+chain] = relativeacc    
        
    filehnd.close()
    return ASAcomplex

def StandardData(self,resname):
    residueDict = {'ALA':107.95,'CYS':134.28,'ASP':140.39,'GLU':172.25,'PHE':199.48,'GLY':80.1,'HIS':182.88,'ILE':175.12,\
                    'LYS':200.81,'LEU':178.63,'MET':194.15,'ASN':143.94,'PRO':136.13,'GLN':178.5,'ARG':238.76,'SER':116.5,\
                    'THR':139.27,'VAL':151.44,'TRP':249.36,'TYR':212.76}

    if residueDict.has_key(resname):
        return residueDict[resname]
    else:
        return -1
##############################################Hotspot Creation##############################################################
def hot_spot_writer(interface, hot_spot_dict1, hot_spot_dict2):
    hotspotWriter = open(f"template/hotspot/{interface}.txt","w")
    hots = str(len(hot_spot_dict1)+len(hot_spot_dict2))
    hotspotWriter.writelines("# "+interface+".pdb "+ hots+" "+hots+"\n")
    for key in sorted(hot_spot_dict1.iterkeys()):
        hotspotWriter.writelines(hot_spot_dict1[key]+"\n")
    for key in sorted(hot_spot_dict2.iterkeys()):
        hotspotWriter.writelines(hot_spot_dict2[key]+"\n")

    hotspotWriter.close()

def hotspot_creator(interface):
        hot_spot_dict1 = {}
        hot_spot_dict2 = {}
        #naccess calls for all interfaceList
        try:
            asa_complex = asa_complex(interface)
        except KeyError:
            print("Error in creating hotspot for "+interface)
            return
        pdbId = interface[0:4]
        chain1 = interface[4]
        chain2 = interface[5]
        pair_potential = contact_potentials(pdbId, chain1, chain2)
        for item in self.left:
            try:
                rel_comp_asa = asa_complex[item]
                potential = pair_potential[pdbId+item]

                if rel_comp_asa <= 20.0 and abs(potential) >= 18.0:
                    key = int(item[1:len(item)-1])
                    hot_spot_dict1[key] = item[-1]+"."+item[0]+"."+item[1:len(item)-1]+" "+item[0]
            except KeyError:
                continue
        for item in self.right:
            try:
                rel_comp_asa = asa_complex[item]
                potential = pair_potential[pdbId+item]
                if rel_comp_asa <= 20.0 and abs(potential) >= 18.0:
                    key = int(item[1:len(item)-1])
                    hot_spot_dict2[key] = item[-1]+"."+item[0]+"."+item[1:len(item)-1]+" "+item[0]
            except KeyError:
                continue

        hot_spot_writer(interface, hot_spot_dict1, hot_spot_dict2)