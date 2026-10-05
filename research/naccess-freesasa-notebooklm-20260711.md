# NACCESS and FreeSASA Research Handoff

NotebookLM notebook: `PRISM Pipeline Protein-DNA Tool Integration Research`

NotebookLM searched and imported ten web sources on 2026-07-11. The most
relevant sources were:

- NACCESS official page: https://www.bioinf.manchester.ac.uk/naccess/nacwelcome.html
- FreeSASA paper: https://arxiv.org/pdf/1601.06764
- FreeSASA CLI documentation: https://freesasa.github.io/1.1/CLI.html
- FreeSASA Python documentation: https://freesasa.github.io/python/intro.html
- FreeSASA RSA output documentation: https://freesasa.github.io/python/functions.html
- FreeSASA implementation documentation: https://github.com/mittinatten/freesasa/blob/master/doc/doxy-main.md

## Implementation Findings

- NACCESS uses Lee-Richards solvent-accessible surface calculations with a
  typical 1.4 Angstrom water probe.
- FreeSASA supports Lee-Richards and Shrake-Rupley calculations and provides
  Python and command-line interfaces.
- FreeSASA can emit an NACCESS-like RSA format, but relative values depend on
  the selected radii/reference table. Missing relative values may be written
  as `N/A` instead of NACCESS's `-99.9`.
- FreeSASA's default radii are not automatically identical to NACCESS. The
  backend therefore uses explicit Lee-Richards/1.4 Angstrom settings and then
  applies PRISM's existing `standard_data` residue normalization.
- FreeSASA and NACCESS have documented differences in CA main-chain/side-chain
  classification and nucleic-acid classification. These must be considered in
  protein-DNA comparisons.

## Local Decision

NACCESS remains the default backend for reproducibility with existing PRISM
results. FreeSASA is an explicit opt-in backend through
`--surface_backend freesasa` or `PRISM_SURFACE_BACKEND=freesasa`. The selected
FreeSASA interpreter is passed through `--freesasa_python` or
`PRISM_FREESASA_PYTHON` so the main PRISM environment does not need to import
FreeSASA directly.

## Validation Requirements

- Compare NACCESS and FreeSASA absolute residue areas on the same PDB.
- Compare the resulting residue surface sets after the PRISM threshold.
- Record backend, interpreter, algorithm, probe radius, and radii/reference
  configuration with every benchmark.
- Treat prediction counts as backend-dependent until equivalence is measured.
