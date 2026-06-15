# PRISM (Protein Interactions by Structural Matching)

End-to-end pipeline:

1. Download receptor/ligand PDBs listed in `inputs.csv`.
2. Optionally filter and pre-compute template artifacts (interface CA PDBs,
   hotspot tables, residue contact maps).
3. Compute target surface residues with FreeSASA.
4. Structurally align target surfaces to every template interface chain with
   TM-align (default) or GTalign.
5. For each receptor-ligand pair, find templates that match both sides,
   transform each input into the template frame, filter by HotPoint-style
   hotspot support, inter-template contact support, and CA-CA steric clash.
6. Optionally run a Rosetta prepack + local docking refinement on accepted
   complexes.

## Inputs

`inputs.csv` accepts multi-chain identifiers; each value is a 4-character PDB
id followed by one or more chain letters:

```csv
Receptor,Ligand
3i6eE,3i6eF
1fgnHL,1tfhA
```

`1fgnHL` means the receptor is chains H and L of `1fgn` (e.g. an antibody
Fv); `1tfhA` means a single-chain ligand. Multi-chain entries flow through
SASA, alignment, transformation, and refinement.

## External tools

```bash
cd external_tools
g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp
chmod 755 TMalign
cd ..
```

For GTalign:
```bash
conda install minmarg::gtalign_mp   # CPU
# or
conda install minmarg::gtalign_gpu  # GPU
```

For Rosetta refinement install [PyRosetta](https://www.pyrosetta.org/).

## Templates

Download the template archive and extract it before running the template
generation stage:

```bash
# Place templates.zip in the repo root, then:
unzip templates.zip
```

## Run

```bash
# First run: rebuild template artifacts then align/transform
python prism.py --generate_templates

# Subsequent runs: reuse cached templates
python prism.py --inputs_csv inputs.csv

# Use GTalign instead of TM-align
python prism.py --aligner gtalign --gtalign_path /path/to/gtalign

# Restrict to top N templates and run Rosetta refinement
python prism.py --template_limit 100 --refine
```

Accepted complexes land in `processed/output/`; per-stage intermediates live
under `processed/`.

## Tests

```bash
pip install pytest
python -m pytest tests/ -v
```

The test suite is hermetic (no internet, no NACCESS, no Rosetta) and covers
the pdb download, SASA, surface extraction, interface/contact/hotspot
generation, TM-align output parsing, transformation filters, Rosetta helper
utilities, and benchmark chain inference.

## Benchmark

See `benchmark/README.md`. The single-pair scorer supports multi-chain
receptor/ligand and infers partner groups from `TER` markers when chain
flags are omitted.
