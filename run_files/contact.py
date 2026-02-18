from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from .utils import vdw_radii_extended, three_to_one, DEFAULT_VDW

def get_contacts(protein):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(protein, f"templates/pdbs/{protein[:4].lower()}.pdb")
    chain_data = {}
    chain_ca = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if not is_aa(residue, standard=True):
                    continue
                res_seq = residue.id[1]
                res_name = residue.get_resname()
                vdw_map = vdw_radii_extended(res_name)
                atoms = []
                ca_coords = None
                for atom in residue:
                    name = atom.get_name()
                    coords = list(atom.get_coord())
                    vdw = vdw_map.get(name, DEFAULT_VDW)
                    atoms.append((name, coords, vdw))
                    if name == "CA":
                        ca_coords = coords
                if atoms and ca_coords is not None:
                    chain_data[chain.id][res_seq] = (res_name, atoms)
                    chain_ca[chain.id][res_seq] = (res_name, ca_coords)
    return chain_data, chain_ca

def create_interaction_list(structure, target):
    chain_list = []
    alternate_location_dict = {}
    icode_dict = {}
    for line in structure:
        if line[:3] == "END":
            break
        elif line[0:4] == "ATOM":
            chain = line[21]
            atomName = line[12:16].strip()
            resName = line[17:20].strip()
            resSeq = line[22:26].strip()
            coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            altLoc = line[16]  # icode and alternate location handled
            icode = line[26]
            key = chain + resSeq
            key2 = key + atomName
            check1 = key2 in alternate_location_dict
            check2 = key in icode_dict
            # these lines are needed to detect alternateLocation or icode differences,
            # the first one is taken...
            if not check1 and not check2:
                alternate_location_dict[key2] = altLoc
                icode_dict[key] = icode
                try:
                    d1 = vdw_radii_extended(resName)
                    d1 = d1.get(atomName, 0)
                    check3 = three_to_one(resName)
                    check4 = atomName[0]
                    chain_list.append([target, chain, coordinates, d1, check3, check4, resName, resSeq])
                except:
                    continue
            elif not check1 and check2:
                alternate_location_dict[key2] = altLoc
                if icode_dict[key] == icode:
                    try:
                        d1 = vdw_radii_extended(resName)
                        d1 = d1.get(atomName, 0)
                        check3 = three_to_one(resName)
                        check4 = atomName[0]
                        chain_list.append([target, chain, coordinates, d1, check3, check4, resName, resSeq])
                    except:
                        continue
    return chain_list

def fiberdock_interface_extractor(file_name, file_out, structure1, structure2):
    splited = file_name.split("/")[-1].split("_")
    target1 = splited[1]
    target2 = splited[3]
    chain_list1 = create_interaction_list(structure1, target1)
    chain_list2 = create_interaction_list(structure2, target2)
    interact1 = {}
    interact2 = {}
    contact = {}
    for atom1 in chain_list1:
        for atom2 in chain_list2:
            if atom1[4] != 'X' and atom1[5] != 'H' and atom2[4] != 'X' and atom2[5] != 'H':
                cutoff = atom1[3] + atom2[3] + 0.5
                coordinates1 = atom1[2]
                coordinates2 = atom2[2]
                distance = ((coordinates1[0] - coordinates2[0])**2 + (coordinates1[1] - coordinates2[1])**2 + (coordinates1[2] - coordinates2[2])**2)**0.5
                if distance <= cutoff:
                    key = str(atom1[7]) + str(atom2[7])
                    if key not in contact:
                        contact[key] = (
                            f"{atom1[0]}_{atom1[1]}_{atom1[6]}_{atom1[7]}"
                            "\t<-->\t"
                            f"{atom2[0]}_{atom2[1]}_{atom2[6]}_{atom2[7]}"
                        )
                        interact1[f"{atom1[6]}_{atom1[7]}_{atom1[1]}"] = 1
                        interact2[f"{atom2[6]}_{atom2[7]}_{atom2[1]}"] = 1
                
    # write contacts
    filehnd = open(file_out, "w")
    filehnd.write("Interface Residues Contacts\t\t\n")
    filehnd.write("target1_chain_resName_resNo\t<-->\ttarget2_chain_resName_resNo\n\n")
    for key in contact:
        filehnd.write(contact[key] + "\n")
    filehnd.close()
    os.system("chmod 775 %s" % (file_out))
    try:
        # write flexible refinement result in updated way to fix visualization
        filehnd = open(file_name, "w")
        for line in structure1:
            chain = line[21]
            atomName = line[12:16].strip()
            resName = line[17:20].strip()
            resSeq = line[22:26].strip()
            key = resName + "_" + resSeq + "_" + chain
            if key in interact1:
                line = line[:60] + "  3.00" + line[66:]
            else:
                line = line[:60] + "  1.00" + line[66:]
            filehnd.writelines(line)
        filehnd.writelines("TER\n")
        for line in structure2:
            chain = line[21]
            atomName = line[12:16].strip()
            resName = line[17:20].strip()
            resSeq = line[22:26].strip()
            key = resName + "_" + resSeq + "_" + chain
            if key in interact2:
                line = line[:60] + "  2.00" + line[66:]
            else:
                line = line[:60] + "  4.00" + line[66:]
            filehnd.writelines(line)
        filehnd.writelines("END\n")
        filehnd.close()    
        os.system("chmod 775 %s" % (file_name))
    except Exception as e: 
        print(e)  
