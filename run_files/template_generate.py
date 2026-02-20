from tqdm import tqdm

from .hotspot import hotspot_creator
from .interface import generate_interface
from .contact import get_contacts

def template_generator():
    calculated_templates = []
    with open("templates/checked_templates.txt", "r") as f:
        templates = [line.strip() for line in f.readlines()]
    for template in tqdm(templates):
        try:
            # print(f"Interface Generation for {template} started!!")
            generate_interface(template)
            # print(f"Interface Generation for {template} Finished!!")
            # print(f"HotSpot Generation for {template} Started!!")
            hotspot_creator(template)
            # print(f"HotSpot Generation for {template} Finished!!")
            # print(f"Contact Generation for {template} Started!!")
            get_contacts(template)
            # print(f"Contact Generation for {template} Finished!!")
            calculated_templates.append(template)
        except Exception as e:
            print(f"Error in template generation for {template}: {e}")
            continue
    with open("templates/calculated_templates.txt", "w") as f:
        for template in calculated_templates:
            f.write(f"{template}\n")
    return calculated_templates
