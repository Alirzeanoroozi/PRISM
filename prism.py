# written by Alper Baspinar
from run_files.pdbDownload import PDBdownloader # pdbdownloader
from run_files.templateGenerator import TemplateGenerator # templateGenerator
from run_files.preProcessor import preprocess_input_proteins # preProcessor
from run_files.surfaceExtractor import SurfaceExtractor # surfaceExtractor
from run_files.structuralAlignmentTM import StructuralAligner # structuralAligner
from run_files.transformationFiltering import TransformFilter # transformFiltering
from run_files.flexibleRefinementRosetta import FlexibleRefinement # flexible refinement step

print("PDB download stage started...")
left_targets, right_targets = PDBdownloader()
print('left_targets', left_targets)
print('right_targets', right_targets)
print("PDB download stage finished...")

#checks template
print("Template generation stage started...")
templates = TemplateGenerator().generator()
print("Template generation stage finished...")

print("PreProcess stage started...")
all_targets = list(set(left_targets + right_targets))
preprocess_input_proteins(all_targets)
print('all_targets', all_targets)
print("PreProcess stage finished...")

print("SurfaceExtraction stage started...")
SurfaceExtractor(all_targets).surfaceExtractor()
print("SurfaceExtraction stage finished...")

print("Structural Alignment stage started...")
StructuralAligner(all_targets, templates)
print("Structural Alignment stage finished...")

print("Transformation Filtering stage started...")
TransformFilter(all_targets, templates).transformer()
print("Transformation Filtering stage finished...")

print("Flexible Refinement stage started...")
energy_Structure = FlexibleRefinement().refiner()
print("Flexible Refinement stage finished...")

