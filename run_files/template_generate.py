import os
from Bio.PDB import PDBParser

from pdb_download import download_pdb_file

def template_generator():
    calculated_templates = []
    for template in [line.strip()[:6] for line in open("templates.txt", "r").readlines()]:
        pdb_path = f"pdb/{template[:4]}.pdb"
        if not os.path.exists(pdb_path):
            print(f"Template {template} does not exist, downloading...")
            download_pdb_file(template[:4].lower())
            print(f"Template {template} downloaded")
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
    contact = {}

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
                d1 = vdw_radii_extended(resName).get(atomName, 0)
                check1 = three2one(resName)
                check2 = atomName[0]
                leftValue = three2one(resName)+str(resSeq)+chainId1
                leftContact = chainId1+"."+three2one(resName)+"."+str(resSeq)+"\t"
                chainList1.append([line, resSeq, coordinates, d1,check1, check2, leftValue,leftContact])
                if atomName == 'CA':
                    chainCa1[resSeq] = [line,coordinates]
            elif chainId == chainId2:
                d1 = vdw_radii_extended(resName).get(atomName, 0)
                check1 = three2one(resName)
                check2 = atomName[0]
                rightValue = three2one(resName)+str(resSeq)+chainId2
                rightContact = chainId2+"."+three2one(resName)+"."+str(resSeq)
                chainList2.append([line, resSeq, coordinates, d1, check1, check2, rightValue,rightContact])
                if atomName == 'CA':
                    chainCa2[resSeq] = [line,coordinates]

    for atom1 in chainList1:
        for atom2 in chainList2:
            if atom1[4] != 'X' and atom1[5] != 'H' and atom2[4] != 'X' and atom2[5] != 'H':
                cutoff = atom1[3] + atom2[3] + 0.5
                distance = distance_calculator(atom1[2], atom2[2])
                if distance <= cutoff:
                    left.add(atom1[6])
                    right.add(atom2[6])
                    key = str(atom1[1])+str(atom2[1])
                    if key not in contact:
                        contact[key] = atom1[7]+atom2[7]
                    try:
                        chain_interact1_res[atom1[1]] = chainCa1[atom1[1]]
                    except KeyError:
                        continue
                    try:
                        chain_interact2_res[atom2[1]] = chainCa2[atom2[1]]
                    except KeyError:
                        continue

    #adding nearby residues
    #first chain
    atom1Iterable = chainCa1.keys()
    atom2Iterable = chainInteract1Res.keys()
    for key1 in atom1Iterable:
        if key1 not in atom2Iterable:
                for key2 in atom2Iterable:
                    distance = distance_calculator(chainCa1[key1][1],chain_interact1_res[key2][1])
                    if distance <= cut_off:
                        try:
                            chain_interact1_res[key1] = chainCa1[key1]
                        except KeyError:
                            continue

    #second chain
    atom1Iterable = chainCa2.keys()
    atom2Iterable = chainInteract2Res.keys()
    for key1 in atom1Iterable:
        if key1 not in atom2Iterable:
                for key2 in atom2Iterable:
                    distance = distance_calculator(chainCa2[key1][1],chain_interact2_res[key2][1])
                    if distance <= cut_off:
                        try:
                            chain_interact2_res[key1] = chainCa2[key1]
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

def contact_writer(interface):
    contactWrite = open(f"template/contact/{interface}.txt","w")
    for key in contact.keys():
        contactWrite.writelines(contact[key]+"\n")
    contact = {}

def contact_potentials(pdb_name, chain1, chain2):
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



def contactingResidues(protein, chain1, chain2):
    centeredResidues = center_of_mass(protein, chain1, chain2)
    contactDict = {}
    temp = []
    for res1 in centeredResidues.keys():
        coor1 = centeredResidues[res1]
        for res2 in centeredResidues.keys():
            coor2 = centeredResidues[res2]
            if res1 != res2:
                if res1[-1] != res2[-1]:
                    distance = distance_calculator(coor1,coor2)
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
                        distance = distance_calculator(coor1,coor2)
                        if distance <= 7.0:
                            try:
                                if res1 in contactDict:
                                    temp = contactDict[res1]
                                temp.append(res2)
                                contactDict.update({res1: temp})
                                temp = []
                            except KeyError:
                                continue
    return contactDict

def center_of_mass(protein, chain1, chain2):
    alphaCarbons = fetch_ca_coordinates(protein, chain1, chain2)

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
                if item in coordinateDict:
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
            centerCoordinates[protein + three2one(resName)+residues[3:]] = real_coordinates
        except KeyError:
            continue
    return centerCoordinates

def fetch_ca_coordinates(protein, chain1, chain2):
    alphaCarbons = {}
    pdb_path = f"pdb/{protein}.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(protein, pdb_path)
    
    for model in structure:
        for chain in model:
            if chain.id == chain1 or chain.id == chain2:
                for residue in chain:
                    if not hasattr(residue, "id"):
                        continue
                    resName = residue.get_resname()
                    resSeq = residue.id[1]
                    key = f"{resName}{resSeq}{chain.id}"
                    atomList = []
                    for atom in residue:
                        atomName = atom.get_name()
                        coordinates = list(atom.get_coord())
                        atomList.append([atomName, coordinates])
                    if atomList:
                        alphaCarbons[key] = atomList
    return alphaCarbons

##############################################ASA FILE READING#############################################################

def asa_complex(interface):
    asa_complex = {}
    with open(f"surface_extraction/{interface}.rsa", 'r') as f:
        for rsaline in f.readlines():
            if rsaline[:3] == "RES":
                rsaline = rsaline.strip()
                resR = rsaline[4:7]
                standard = standard_data(resR)
                if standard != -1:
                    resName = three2one(resR)
                    chain = rsaline[8]
                    resPosition = rsaline[9:13].strip()
                    try:
                        absoluteacc = float(rsaline[14:22])
                    except:
                        print(f"naccess output handle error for protein {interface}")
                        return asa_complex
                    relativeacc = absoluteacc * 100 / standard
                    asa_complex[resName + resPosition + chain] = relativeacc    
    
    return asa_complex



##############################################Hotspot Creation##############################################################
def hot_spot_writer(interface, hot_spot_dict1, hot_spot_dict2):
    hotspotWriter = open(f"template/hotspot/{interface}.txt","w")
    hots = str(len(hot_spot_dict1)+len(hot_spot_dict2))
    hotspotWriter.writelines("# "+interface+".pdb "+ hots+" "+hots+"\n")
    for key in sorted(hot_spot_dict1.iterkeys()):
        hotspotWriter.writelines(hot_spot_dict1[key]+"\n")
    for key in sorted(hot_spot_dict2.iterkeys()):
        hotspotWriter.writelines(hot_spot_dict2[key]+"\n")

def hotspot_creator(template):
    protein, chain_id1, chain_id2 = template[:4].lower(), template[4], template[5]
    hot_spot_dict1 = {}
    hot_spot_dict2 = {}
    asa_complex = asa_complex(template)
    pair_potential = contact_potentials(protein, chain_id1, chain_id2)
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