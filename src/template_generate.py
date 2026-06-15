"""Per-template artifact generation: interfaces, hotspots, contacts."""

import os
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
    except Exception as exc:
        return None, f"{template}: {exc}"


def template_generator(input_list="templates/checked_templates.txt",
                      output_list="templates/calculated_templates.txt",
                      max_workers=None):
    if not os.path.exists(input_list):
        raise FileNotFoundError(
            f"Required template list not found: {input_list}. "
            "Run the analyse_pdbs stage first or pass an alternative path."
        )
    with open(input_list) as f:
        templates = [line.strip() for line in f if line.strip()]

    calculated = []
    workers = max_workers or min(8, os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(process_template, t): t for t in templates}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Templates"):
            ok, err = fut.result()
            if ok is not None:
                calculated.append(ok)
            elif err:
                print(err)

    with open(output_list, "w") as f:
        for t in calculated:
            f.write(f"{t}\n")
    return calculated
