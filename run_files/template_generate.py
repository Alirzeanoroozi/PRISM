import os
from tqdm import tqdm

from .pdb_download import download_pdb_file
from .hotspot import hotspot_creator
from .interface import generate_interface, get_contacts

def template_generator():
    calculated_templates = []
    with open("processed/templates.txt", "r") as f:
        templates = [line.strip() for line in f.readlines()]
    for template in tqdm(templates):
        pdb_path = f"templates/pdbs/{template[:4].lower()}.pdb"
        if not os.path.exists(pdb_path):
            print(f"Template {template} does not exist, start downloading...")
            download_pdb_file(template[:4].lower(), "templates/pdbs")
            if not os.path.exists(pdb_path):
                print(f"Failed to download template {template}, skipping...")
                continue
        try:
            # print(f"Interface Generation for {template} started!!")
            generate_interface(template)
            # print(f"Interface Generation for {template} Finished!!")
            # print(f"HotSpot Generation for {template} Started!!")
            hotspot_creator(template)
            # print(f"HotSpot Generation for {template} Finished!!")
            # print(f"Contact Generation for {template} Started!!")
            # get_contacts(template)
            # print(f"Contact Generation for {template} Finished!!")
            calculated_templates.append(template)
        except Exception as e:
            print(f"Error in template generation for {template}: {e.with_traceback()}")
            continue
    return calculated_templates
