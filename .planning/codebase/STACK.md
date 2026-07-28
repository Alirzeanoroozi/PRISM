# Stack

## Overview

PRISM-prescript is a Python command-line pipeline for protein–protein docking,
with a substantial Python benchmark and reproducibility layer. The maintained
entry point is `prism.py`; the repository also retains legacy Python 2 tooling
under `working_version/` for compatibility investigations, not routine changes.

## Runtime and language

- Primary language: Python 3.
- Current validated operational interpreter: `/home/rshadi25/.conda/envs/gtalign_env/bin/python` (Python 3.11, Biopython, GTalign, and PyRosetta on the project host).
- Reproducible base recipes are `environment.yaml` (Python 3.11.13, NumPy 1.26.4, pandas 2.3.3, Biopython 1.84) and the older `environment.yml` (Python 3.10, unpinned scientific packages).
- `environment-pyrosetta.yaml` describes a separate Python 3.11 environment; the licensed PyRosetta distribution is intentionally installed separately.
- There is no `pyproject.toml`, `setup.py`, or package build step. Imports assume execution from the repository root.

## Key Python libraries

- Biopython: PDB parsing, chain/residue access, sequence/alignment utilities.
- NumPy: coordinate arrays, transformations, numerical metrics.
- pandas: CSV input and tabular benchmark processing.
- `tqdm`: template-generation progress reporting.
- FreeSASA: optional surface-area backend (`freesasa==2.2.1` in the maintained recipe).
- DockQ: optional scoring dependency (`dockq==2.1.3` in the maintained recipe).
- matplotlib and reportlab: analysis/plotting and report-generation utilities.
- PyRosetta: optional licensed refinement backend, loaded lazily by `src/pyrosetta_refinement.py`.

## External executables

The pipeline shells out to tools that are not Python packages: TMalign,
GTalign CPU/GPU, NACCESS, Rosetta 2022.42, MultiProt, and FiberDock. Their
paths, versions, licenses, architecture, and runtime libraries are part of
the execution environment rather than a lockfile. Benchmark scripts add
Slurm, optional headless PyMOL, and evaluator-specific subprocesses.

## Configuration model

Small CLI surfaces are defined with `argparse` in `prism.py` and benchmark
scripts. Many operational controls are environment variables, including
`PRISM_INPUTS_CSV`, `PRISM_SURFACE_BACKEND`, aligner/refiner paths, filtering
thresholds, audit destinations, and Slurm runner settings. Defaults are
source-defined and often relative to the current working directory.

## Version and reproducibility notes

The repository contains multiple historical recipes and installed environments;
do not infer that every environment directory is current. For production-like
runs, use the documented `gtalign_env` interpreter, explicit GTalign paths,
Rosetta 2022.42 module loading, isolated work directories, and recorded stage
and provenance manifests. Diagnostic threshold overrides must remain separate
from stable defaults.
