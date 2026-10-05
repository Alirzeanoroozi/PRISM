# Refinement-tool research handoff

This note records the NotebookLM query used during the PRISM pipeline
implementation. It is a literature/tool-capability handoff, not evidence
that a tool is installed or functional locally.

NotebookLM sources queried:

- `PRISM benchmarking - main project`
- `Structure Alignment - main project`
- `Thesis_papers_prism_histone - main project`

The query asked for deployable protein-protein energy/refinement tools,
including CPU/GPU requirements, Python compatibility, licensing, multichain
support, and suitability for controlled benchmarking.

## Source-grounded findings

- FiberDock is the refinement/scoring component described by the PRISM
  protocol. It consumes rigid-body solutions, uses CHARMM-based energies,
  side-chain rotamers, and normal modes, and produces refined structures and
  global energies.
- The PRISM protocol describes Python 2-era execution plus external NACCESS,
  MultiProt, and FiberDock dependencies.
- Rosetta is a well-established docking/refinement suite and is already
  available as the cluster `rosetta/2022.42` module in this project’s observed
  environment.
- PyRosetta was not documented by the queried NotebookLM sources; its local
  installation and licensing/API behavior must be verified independently
  before use.
- HADDOCK and FireDock appear as alternatives in the sources, but the queried
  material does not establish a reproducible local installation for this
  benchmark.

## Decision for this repository

Use the existing Rosetta module for the primary refinement arm. Treat
PyRosetta, HADDOCK, and FireDock as exploratory until their versions,
installation, input contract, output contract, and benchmark behavior are
independently frozen. Treat FiberDock as primary only if its complete
hydrogen/NMA/refinement/export chain passes the local capability contract;
the current evidence proves only energy-only execution.
