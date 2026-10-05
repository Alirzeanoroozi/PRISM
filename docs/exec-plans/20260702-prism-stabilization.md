# Stabilize PRISM-prescript smoke execution

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` up to date as work proceeds.

## Purpose / Big Picture

Make `PRISM-prescript` safe to test as a current-pipeline repo instead of only as a benchmark repo. Success means:
- bounded local checks still pass
- a repo-local smoke runner can stage a simple `1FGNH` vs `1TFHA` run
- known current-pipeline helper bugs do not block the first successful passed pair

## Progress

- [x] Reload project memory and inspect current runnable entry points
- [x] Identify concrete current-pipeline blockers in `src/rosetta_refinement.py`
- [x] Patch current-pipeline helper bugs and add a bounded smoke runner
- [x] Re-run stable checks and the smoke runner
- [x] Summarize remaining risks

## Surprises & Discoveries

- Observation: `prism.py --generate_templates False` is not safe because `argparse` with `type=bool` treats many non-empty strings as true.
  Evidence: `prism.py` currently defines `parser.add_argument("--generate_templates", type=bool, default=False)`.
- Observation: Rosetta partner parsing in `src/rosetta_refinement.py` uses filename character offsets instead of real chain IDs.
  Evidence: the current code uses `passed0[4]` and `passed1[4]`.
- Observation: the Rosetta stage also calls `get_contacts()` with a signature that does not exist in this repo.
  Evidence: `src/contact.py` exports `get_contacts(template)` only, while `src/rosetta_refinement.py` passes four arguments.
- Observation: the current TMalign path wrote only per-target CSV summaries, while `src/transformation.py` loaded per-pair JSONs.
  Evidence: the first smoke rerun failed at `processed/alignment/1FGNH_1kcaCH_C.json` even though `processed/alignment/1FGNH.csv` and `processed/alignment/1TFHA.csv` had been written.
- Observation: the bundled NACCESS shell wrapper still hardcoded `/scratch/rshadi25/GitHub/PRISM/external_tools/naccess`.
  Evidence: the smoke run failed with `unable to assign a vdw radii file` until the wrapper and Python launcher were corrected.

## Decision Log

- Decision: keep the stabilization additive and local to `PRISM-prescript`.
  Rationale: the repo is mid-rebase and already has many unrelated changes, so this pass should avoid broad refactors.
  Date/Author: 2026-07-02 / Codex
- Decision: add a repo-local smoke runner instead of changing benchmark scripts into a current-pipeline launcher.
  Rationale: benchmark validation and root pipeline execution serve different purposes.
  Date/Author: 2026-07-02 / Codex

## Outcomes & Retrospective

The stabilization pass succeeded. `PRISM-prescript` now has:
- a stable benchmark-side validation bundle (`benchmark/scripts/run_stable_checks.sh`)
- a bounded repo-local current-pipeline smoke runner (`benchmark/scripts/run_prism_pipeline_smoke.sh`)
- patched current-pipeline helpers for boolean CLI parsing, template-list fallback loading, NACCESS staging, TMalign JSON emission, and Rosetta partner/contact handling

The validated smoke case still ends with `Passed pairs 0`, so the remaining question is biological/filtering behavior rather than bootstrap failure.

## Context and Orientation

Key files for this stabilization pass:
- `prism.py`
- `src/pdb_download.py`
- `src/transformation.py`
- `src/contact.py`
- `src/rosetta_refinement.py`
- `benchmark/scripts/run_stable_checks.sh`

Relevant assets already present in the repo:
- `inputs.csv` with `1FGNH,1TFHA`
- `templates/interfaces/1kcaCH_C_int.pdb`
- `templates/interfaces/1kcaCH_H_int.pdb`
- `templates/interfaces_lists/1kcaCH.json`
- `processed/pdbs/1fgn.pdb`
- `processed/pdbs/1tfh.pdb`

## Plan of Work

Patch the current-pipeline helper layer first, then add a small smoke runner that stages a temporary working directory with one input pair and one template. Keep the smoke runner independent from the benchmark result tree and make it succeed even when no pair passes filtering, as long as the pipeline stages execute cleanly.

## Concrete Steps

1. Working directory: `/scratch/rshadi25/GitHub/PRISM-prescript`
   Patch `prism.py`, `src/pdb_download.py`, `src/transformation.py`, `src/contact.py`, and `src/rosetta_refinement.py`.
2. Working directory: `/scratch/rshadi25/GitHub/PRISM-prescript`
   Add `benchmark/scripts/run_prism_pipeline_smoke.sh` and a small helper-test script.
3. Working directory: `/scratch/rshadi25/GitHub/PRISM-prescript`
   Run:
   - `bash benchmark/scripts/run_stable_checks.sh`
   - `python benchmark/scripts/test_prism_pipeline_helpers.py`
   - `bash benchmark/scripts/run_prism_pipeline_smoke.sh`

## Validation and Acceptance

Acceptance criteria:
- stable checks still pass
- helper tests pass
- the smoke runner exits `0`, produces alignment outputs, and reaches the transformation stage without crashing
- no direct Rosetta helper crash remains in the current code path

## Idempotence and Recovery

- All new validation artifacts should go under `tmp/agent/` or a temporary smoke workspace.
- The smoke runner should use a fresh temporary directory each run.
- No existing benchmark outputs should be modified in place.

## Artifacts and Notes

- Stable-check entry point: `benchmark/scripts/run_stable_checks.sh`
- New smoke-run artifacts: under `tmp/agent/prism-pipeline-smoke-*`
- Most recent successful smoke workspace: `tmp/agent/prism-pipeline-smoke-kqowe4`

## Interfaces and Dependencies

- Python interpreter should be a repo-compatible env with `numpy`, `pandas`, and `biopython`.
- NACCESS is invoked via `external_tools/naccess/naccess`.
- TM-align is invoked via `external_tools/TMalign`.
- Rosetta command names remain external dependencies, but the smoke case is expected to finish without using them when no pair passes thresholds.
