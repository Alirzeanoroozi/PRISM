# Add an optional PyRosetta refinement adapter and environment probe

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` current while implementing.

## Purpose / Big Picture

Add an opt-in PyRosetta refinement boundary that can report an unavailable installation without importing PyRosetta at module import time, and can refine an explicitly supplied PDB when PyRosetta is available. Add a standalone probe script that emits JSON suitable for environment auditing. The existing external-CLI Rosetta backend remains the default and unchanged by this work.

## Progress

- [x] Route to `PRISM-prescript`, read project guidance and memory, and inspect the existing Rosetta backend and environment files.
- [x] Write focused failing tests for lazy import, fail-closed status, provenance metadata, and probe CLI behavior.
- [x] Implement `src/pyrosetta_refinement.py` without adding a dependency or changing `src/rosetta_refinement.py`.
- [x] Implement `benchmark/scripts/probe_pyrosetta_environment.py`.
- [x] Run focused tests and review the diff for isolation from unrelated edits.
- [x] Update project memory with the durable optional-backend decision and validation result.

## Surprises & Discoveries

- The repository contains both `environment.yml` and `environment.yaml`; the requested `environment.yaml` currently has no PyRosetta dependency, and neither file will be modified.
- `prism.py` imports the existing CLI backend directly, so an independent adapter module is sufficient to keep the default backend unchanged.
- The current worktree contains unrelated staged, modified, and untracked files, including changes to `src/rosetta_refinement.py`; edits must be additive and must not normalize or revert that file.

## Decision Log

- Decision: Keep the adapter opt-in and expose no automatic fallback to CLI Rosetta.
  Rationale: A missing or unusable PyRosetta installation must remain an explicit unavailable/failed result rather than silently changing scientific backends.
  Date/Author: 2026-07-14 / Codex
- Decision: Use lazy `importlib.import_module("pyrosetta")` only inside the probe/refinement call path.
  Rationale: The standard benchmark environment does not declare or guarantee PyRosetta.
  Date/Author: 2026-07-14 / Codex
- Decision: Return JSON-serializable result dictionaries and write a sidecar metadata record for refinement attempts.
  Rationale: This matches the repository's provenance-oriented scripts and makes package/version/import errors, hashes, command metadata, and environment metadata inspectable without requiring PyRosetta.
  Date/Author: 2026-07-14 / Codex

## Outcomes & Retrospective

The adapter and CLI probe are implemented and the focused suite passes under `gtalign_env` (`5 passed`), with two existing Rosetta helper tests also passing. The current standard Python lacks `pytest`, so verification used `/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest`. The live probe reports `status="unavailable"` and `ModuleNotFoundError: No module named 'pyrosetta'`, as expected. A repository-wide `git diff --check` remains noisy because unrelated pre-existing edits in `src/rosetta_refinement.py` and `src/surface_extract.py` contain trailing whitespace; those files were not changed by this task.

## Context and Orientation

The default Rosetta path is `src/rosetta_refinement.py`, imported by `prism.py` and using external `docking_prepack_protocol.static.linuxgccrelease` and `docking_protocol.static.linuxgccrelease` commands. The new module will not be imported by `prism.py`. The probe script lives at `benchmark/scripts/probe_pyrosetta_environment.py` and adds the repository root to `sys.path` using the same pattern as other benchmark scripts.

## Plan of Work

First define the public result contract in tests. Then implement a lazy environment probe and an explicit adapter with `refine`, `refine_pair`, and a small batch helper. The adapter will hash input files before execution, hash successful output files afterward, capture safe command/environment metadata, write a sidecar JSON record, and return `status="unavailable"` when import fails. It will never import or invoke `src.rosetta_refinement` as a fallback. Finally, the CLI probe will print/write the same probe record and optionally fail its process only when `--require-available` is requested.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, add focused tests and run them once to confirm the expected failures.
2. From `/scratch/rshadi25/GitHub/PRISM-prescript`, implement `src/pyrosetta_refinement.py` with no top-level PyRosetta import and no environment-file change.
3. From `/scratch/rshadi25/GitHub/PRISM-prescript`, implement the probe CLI and run it in the current environment; preserve its JSON as a temporary validation artifact only if needed.
4. From `/scratch/rshadi25/GitHub/PRISM-prescript`, run the focused tests, inspect `git diff --stat` and critical regions, and confirm `prism.py` and `src/rosetta_refinement.py` are not changed by this task.

## Validation and Acceptance

Acceptance requires focused tests showing: importing the adapter does not import PyRosetta; a missing or import-error PyRosetta returns `available=false` and `status="unavailable"`; no CLI Rosetta fallback is called; successful injected/fake refinement records input/output SHA-256 hashes and metadata; the probe CLI emits valid JSON and leaves the default environment declaration free of PyRosetta. A real PyRosetta run is not required or claimed because the dependency is intentionally unverified and absent from the environment declaration.

## Idempotence and Recovery

The adapter writes only its explicitly requested output and sidecar metadata paths. Re-running with the same paths overwrites only those opt-in artifacts. If implementation needs to be backed out, remove the new module, probe script, tests, plan, and memory note; do not revert unrelated worktree changes or the existing Rosetta backend.

## Artifacts and Notes

Primary files: `src/pyrosetta_refinement.py`, `benchmark/scripts/probe_pyrosetta_environment.py`, and focused tests under `tests/`. No PyRosetta package, lockfile, or environment declaration will be added.

## Interfaces and Dependencies

The standard library is sufficient for probing, metadata, hashing, JSON, and subprocess-free CLI operation. PyRosetta is an optional runtime discovered only through a lazy import. The adapter's real execution path expects the conventional `pyrosetta.init`, `pose_from_pdb`, `get_fa_scorefxn`, and docking protocol APIs; any runtime/API error is returned as `status="failed"` with no CLI fallback.
