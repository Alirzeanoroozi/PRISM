import os
from Bio.PDB import PDBParser

from pdb_download import download_pdb_file
from analyse_pdbs import analyse_pdb
from hotspot import hotspot_creator
from interface import generate_interface

def template_generator():
    calculated_templates = []
    for template in [line.strip()[:6] for line in open("templates.txt", "r").readlines()]:
        pdb_path = f"pdb/{template[:4]}.pdb"
        if not os.path.exists(pdb_path):
            print(f"Template {template} does not exist, downloading...")
            download_pdb_file(template[:4].lower())
            print(f"Template {template} downloaded")
            
        if analyse_pdb(pdb_path):
            print(f"Template {template} is not valid, skipping...")
            continue
        try:
            print(f"Interface Generation for {template} started!!")
            generate_interface(template)
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