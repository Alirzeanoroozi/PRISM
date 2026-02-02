# Keys in this file are like this: GLN_682_A --> residue name, residue number, chain id

import os
import pickle
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from utils import distance_calculator, three2one, ATOM_DICT, PAIR_POTENTIAL
from naccess_utils import get_asa_complex

RELATIVE_ASA_THRESHOLD = 20.0
CONTACT_POTENTIAL_THRESHOLD = 18.0
RESIDUE_DISTANCE_THRESHOLD = 3
CONTACT_DISTANCE_THRESHOLD = 7.0

HOTSPOT_DIR = "templates/hotspots"
os.makedirs(HOTSPOT_DIR, exist_ok=True)

RSA_DIR = "templates/rsas"
os.makedirs(RSA_DIR, exist_ok=True)

def hotspot_creator(template):
    protein, chain_id1, chain_id2 = template[:4].lower(), template[4], template[5]
    asa_complex = get_asa_complex(template, RSA_DIR)
    contact_potentials = get_contact_potentials(protein, chain_id1, chain_id2)

    hotspot_dict = {chain_id1: [], chain_id2: []}
    for item in asa_complex:
        chain = item.split("_")[2]
        if asa_complex[item] <= RELATIVE_ASA_THRESHOLD and abs(contact_potentials[item]) >= CONTACT_POTENTIAL_THRESHOLD:
            hotspot_dict[chain].append((item.split("_")[1], item.split("_")[0])) # (residue number, residue name)

    with open(f"{HOTSPOT_DIR}/{template}.pkl", "wb") as f:
        pickle.dump(hotspot_dict, f)

def get_contact_potentials(protein, chain1, chain2):
    contact_dict = contacting_residues(protein, chain1, chain2)
    contact_potentials = {}
    for res1 in contact_dict:
        total_contact = 0
        for res2 in contact_dict[res1]:
            aa1 = three2one(res1.split("_")[0])
            aa2 = three2one(res2.split("_")[0])
            if aa1 != 'X' and aa2 != 'X':
                temp = [aa1, aa2]
                temp.sort()
                total_contact += PAIR_POTENTIAL[f"{temp[0]}-{temp[1]}"]
        contact_potentials[res1] = total_contact if len(contact_dict[res1]) != 0 else 0.0
    return contact_potentials

def contacting_residues(protein, chain1, chain2):
    cms = center_of_mass(protein, chain1, chain2)
    contact_dict = {}
    for res1, coor1 in cms.items():
        for res2, coor2 in cms.items():
            if res1 != res2:
                if res1[-1] != res2[-1]: # Check if the residues are in different chains
                    if distance_calculator(coor1, coor2) <= CONTACT_DISTANCE_THRESHOLD:
                        if res1 in contact_dict:
                            contact_dict[res1].append(res2)
                        else:
                            contact_dict[res1] = [res2]
                elif res1[-1] == res2[-1]:
                    if abs(int(res1.split("_")[1]) - int(res2.split("_")[1])) > RESIDUE_DISTANCE_THRESHOLD:
                        if distance_calculator(coor1, coor2) <= CONTACT_DISTANCE_THRESHOLD:
                            if res1 in contact_dict:
                                contact_dict[res1].append(res2)
                            else:
                                contact_dict[res1] = [res2]
        if res1 not in contact_dict:
            contact_dict[res1] = []
    return contact_dict

def center_of_mass(protein, chain1, chain2):
    all_atoms = fetch_all_atoms_coordinates(protein, chain1, chain2)

    center_coordinates = {}
    for residue in all_atoms:
        atoms_coordinates = all_atoms[residue].values()
        center_coordinates[residue] = [sum(coord[i] for coord in atoms_coordinates) / len(atoms_coordinates) for i in range(3)]
    return center_coordinates

def fetch_all_atoms_coordinates(protein, chain1, chain2):
    pdb_path = f"pdbs/{protein}.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(protein, pdb_path)
    model = structure[0]

    all_atoms = {}    
    for chain in model:
        if chain.id == chain1 or chain.id == chain2:
            for residue in chain:
                if not is_aa(residue, standard=True): # Check for only ATOM residues (standard amino acids) only include standard amino acids (ATOM, not HETATM)
                    continue
                residue_name = residue.get_resname()
                residue_number = residue.id[1]
                chain_id = chain.id
                ###### TODO: fix the ATOM_DICT issue ###########
                # all_atoms[f"{residue_name}_{residue_number}_{chain_id}"] = {atom.get_name(): list(atom.get_coord()) for atom in residue if atom.get_name() in ATOM_DICT[residue_name]}
                all_atoms[f"{residue_name}_{residue_number}_{chain_id}"] = {atom.get_name(): list(atom.get_coord()) for atom in residue}
                # assert len(all_atoms[f"{residue_name}_{residue_number}_{chain_id}"]) == len(ATOM_DICT[residue_name])
    return all_atoms


if __name__ == "__main__":
    template = "1a28AB"
    protein = template[:4].lower()
    # asa_complex = get_asa_complex(protein)
    # for k, v in list(asa_complex.items())[0:5]:
    #     print(k, v)
    # print(len(asa_complex))
    
    # all_atoms = fetch_all_atoms_coordinates(protein, template[4], template[5])
    # for k, v in all_atoms.items():
    #     print(k, v)
    
    cm = center_of_mass(protein, template[4], template[5])
    for k, v in list(cm.items())[0:5]:
        print(k, v)
    print(len(cm))
        
    contact_dict = contacting_residues(protein, template[4], template[5])
    for k, v in list(contact_dict.items())[0:5]:
        print(k, v)
    print(len(contact_dict))
    
    contact_potential_dict = contact_potentials(protein, template[4], template[5])
    for k, v in list(contact_potential_dict.items())[0:5]:
        print(k, v)
    print(len(contact_potential_dict))
    
    hotspot_creator(template)
    with open(f"{HOTSPOT_DIR}/{template}.pkl", "rb") as f:
        hotspot_dict = pickle.load(f)
    for k, v in hotspot_dict.items():
        print(k, v)
    print(len(hotspot_dict))
