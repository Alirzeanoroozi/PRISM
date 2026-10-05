# Rosetta-Output-Specific Benchmark Scripts

This folder contains the PRISM/Rosetta-output-specific benchmark pipeline.

Use these files when the input models come from PRISM-style `rosetta_output_*`
folders and the model filenames encode `PDB ID 1` / `PDB ID 2` partners.

Main entry points:

- `analyze_prism_all_benchmarks.py`
- `analyze_prism_rigid_results.py`
- `fix_model_chain_names.py`
- `submit_prism_analysis_all.sbatch`

Rosetta-specific helper/reporting scripts are also kept here so the top-level
`benchmark/scripts` folder stays focused on general scoring utilities such as
`dockq.py`, `irmsd.py`, and `score_single_prism_pair.py`.

