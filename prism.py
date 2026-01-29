from run_files.pdb_download import pdb_downloader
from run_files.template_generate import template_generator
from run_files.surfaceExtractor import SurfaceExtractor
from run_files.structuralAlignmentTM import StructuralAligner
from run_files.transformationFiltering import TransformFilter
from run_files.flexibleRefinementRosetta import FlexibleRefinement

print("PDB download stage started...")
targets = pdb_downloader()
print('targets', targets)
print("PDB download stage finished...")

# check template
print("Template generation stage started...")
templates = template_generator()
print("Template generation stage finished...")

print("SurfaceExtraction stage started...")
SurfaceExtractor(targets).surfaceExtractor()
print("SurfaceExtraction stage finished...")

print("Structural Alignment stage started...")
StructuralAligner(targets, templates)
print("Structural Alignment stage finished...")

print("Transformation Filtering stage started...")
TransformFilter(targets, templates).transformer()
print("Transformation Filtering stage finished...")

print("Flexible Refinement stage started...")
energy_Structure = FlexibleRefinement().refiner()
print("Flexible Refinement stage finished...")

