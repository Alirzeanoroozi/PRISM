import argparse

from src.pdb_download import pdb_downloader
from src.analyse_pdbs import run_analysis
from src.template_generate import template_generator
from src.surface_extract import extract_surfaces
from src.alignment import align
from src.alignment_gtalign import align_gtalign
from src.transformation import transformer
from src.rosetta_refinement import refiner

def main(args):
    alignment_output_dir = "processed/alignment"
    print("PDB download stage started...")
    receptor_targets, ligand_targets = pdb_downloader()
    targets = receptor_targets + ligand_targets
    for r, l in zip(receptor_targets, ligand_targets):
        print(r, "->", l)
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
        with open("templates/calculated_templates.txt", "r") as f:
            templates = [line.strip() for line in f.readlines()]
        print("Templates loaded, templates length", len(templates))

    print("Surface extraction stage started...")
    extract_surfaces(targets)
    print("Surface extraction stage finished...")

    print("Structural alignment stage started...")
    # Local change: optional backend switch; default remains original TMalign pipeline.
    if args.aligner == "tmalign":
        align(targets, templates)
    else:
        alignment_output_dir = "processed/alignment_gtalign"
        align_gtalign(
            targets,
            templates,
            gtalign_path=args.gtalign_path,
            output_dir=alignment_output_dir,
            dev_min_length=args.gtalign_dev_min_length,
            pre_score=args.gtalign_pre_score,
            speed=args.gtalign_speed,
            refinement=args.gtalign_refinement,
        )
    print("Structural alignment stage finished...")

    print("Transformation filtering stage started...")
    passed_pairs = transformer(templates, alignment_dir=alignment_output_dir)
    print("Passed pairs", len(passed_pairs))
    for pair in passed_pairs:
        print(pair)
    print("Transformation filtering stage finished...")

    print("Rosetta refinement stage started...")
    refiner(passed_pairs)
    print("Rosetta refinement stage finished...")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate_templates", type=bool, default=False)
    # Local change: keep TMalign as default; allow GTalign on user request.
    parser.add_argument("--aligner", choices=["tmalign", "gtalign"], default="tmalign")
    parser.add_argument("--gtalign_path", type=str, default="gtalign")
    parser.add_argument("--gtalign_dev_min_length", type=int, default=3)
    parser.add_argument("--gtalign_pre_score", type=float, default=0.0)
    parser.add_argument("--gtalign_speed", type=int, default=0)
    parser.add_argument("--gtalign_refinement", type=int, default=3)
    args = parser.parse_args()
    main(args)
