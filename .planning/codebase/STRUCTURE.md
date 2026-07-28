# Structure

## Repository map

- `prism.py` — maintained command-line orchestration entry point.
- `src/` — current Python pipeline modules and evaluation helpers.
  - `alignment*.py` — TMalign, GTalign, and MultiProt adapters.
  - `pdb_download.py`, `surface_extract.py`, `naccess_utils.py`,
    `freesasa_runner.py` — input and surface stages.
  - `transformation.py`, `template_filtering.py`, `candidate_audit.py` —
    filtering, transformation, and provenance contracts.
  - `*_refinement.py` — Rosetta, PyRosetta, and FiberDock backends.
  - `compare.py`, `eval/` — DockQ/iRMSD evaluation.
  - `candidate_ranker.py`, `candidate_selector.py`, `ranking_*.py` — optional
    deterministic ranking and offline ranking analyses.
  - `eda/` — exploratory sequence and plotting utilities.
- `benchmark/` — benchmark data conventions, job manifests, scoring tools,
  Slurm templates, and retained processed evidence.
  - `scripts/` — manifest builders, collectors, validators, scorers, audits,
    replays, and benchmark-specific tests.
  - `jobs/` — batch configuration and job templates.
  - `originals/` — archived benchmark inputs and source material.
- `tests/` — pytest-oriented unit and contract tests for current code and
  benchmark tooling; some specialized environment tests are conditional.
- `external_tools/` — checked-in or staged third-party binaries and tool
  payloads, including NACCESS, TMalign, MultiProt, and FiberDock.
- `working_version/` — legacy/reference trees; compatibility-only.
- `new_template/`, `template_old/`, `template_files/` — template assets and
  historical staging views.
- `docs/` — stable run instructions, validation reports, execution plans,
  chronology, and decision evidence.
- `tmp/agent/`, `processed/`, `templates/`, `logs/` — generated run state and
  artifacts. They are operational evidence or scratch, not normal source
  modules; use isolated run subdirectories.
- `.planning/` — learnship project state and codebase reference documents.

## Naming and path conventions

Python modules use lowercase `snake_case`; functions and local variables use
the same style, while constants are uppercase. CLI options retain historical
underscore spellings in `prism.py` (`--generate_templates`,
`--surface_backend`, `--gtalign_path`) alongside newer hyphenated options.

PDB selectors preserve a four-character PDB ID plus chain suffix. The canonical
normalizer sorts and de-duplicates chain IDs, while raw selectors remain in
benchmark inputs for provenance. Generated PDB filenames encode partner,
template, orientation, and sometimes `_L`/`_R` side markers.

## Where to make changes

Change current pipeline behavior in `prism.py` and the corresponding `src/`
adapter/module. Add benchmark collection and audit behavior under
`benchmark/scripts/` with focused tests under `tests/`. Update stable operating
instructions in `docs/STABLE_PIPELINE.md` when defaults or contracts change.
Do not edit generated outputs, raw benchmark inputs, project memory, or the
legacy compatibility tree without an explicit, scoped reason.
