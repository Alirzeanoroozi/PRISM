from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from .hotspot import hotspot_creator
from .interface import generate_interface
from .contact import get_contacts

def process_template(template):
    try:
        generate_interface(template)
        hotspot_creator(template)
        get_contacts(template)
        return template, None
    except Exception as e:
        return None, f"Error in template generation for {template}: {e}"

def template_generator():
    calculated_templates = []

    with open("templates/checked_templates.txt", "r") as f:
        templates = [line.strip() for line in f.readlines()]

    futures = []
    with ThreadPoolExecutor() as executor:
        for template in templates:
            futures.append(executor.submit(process_template, template))

        for future in tqdm(as_completed(futures), total=len(templates)):
            result, error = future.result()
            if result is not None:
                calculated_templates.append(result)
            if error:
                print(error)

    with open("templates/calculated_templates.txt", "w") as f:
        for template in calculated_templates:
            f.write(f"{template}\n")
    return calculated_templates
