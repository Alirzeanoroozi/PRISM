import argparse

from run_files.pdb_download import pdb_downloader
from run_files.analyse_pdbs import run_analysis
from run_files.template_generate import template_generator
from run_files.surface_extract import extract_surfaces
from run_files.alignment import align
from run_files.transformation import transformer
from run_files.rosetta_refinement import FlexibleRefinement

def main(args):
    print("PDB download stage started...")
    targets = pdb_downloader()
    print("targets", targets)
    print("PDB download stage finished...")

    if args.generate_templates:
        print("Filtering templates stage started...")
        results, filtered_count, failed_count = run_analysis()
        print(f"Analysed {len(results)} template files.")
        print(f"Filtered {filtered_count} templates")
        print(f"Failed {failed_count} templates")
        print("Filtering templates stage finished...")

        print("Template generation stage started...")
        templates = template_generator()
        print("Templates generated, templates length", len(templates))
        print("Template generation stage finished...")
    else:
        with open("processed/templates.txt", "r") as f:
            templates = [line.strip() for line in f.readlines()]
        print("Templates loaded, templates length", len(templates))

    print("Surface extraction stage started...")
    extract_surfaces(targets)
    print("Surface extraction stage finished...")

    print("Structural alignment stage started...")
    align(targets, templates)
    print("Structural alignment stage finished...")

    print("Transformation filtering stage started...")
    passed_pairs = transformer(templates)
    print("Passed pairs", len(passed_pairs))
    for pair in passed_pairs:
        print(pair)
    print("Transformation filtering stage finished...")

    # print("Flexible refinement stage started...")
    # energy_structure = FlexibleRefinement().refiner()
    # print("Flexible refinement stage finished...")
    # return energy_structure

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate_templates", type=bool, default=False)
    args = parser.parse_args()
    main(args)
