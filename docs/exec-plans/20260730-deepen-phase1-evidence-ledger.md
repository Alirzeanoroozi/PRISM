# Deepen the Phase 1 run evidence ledger

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`,
`Decision Log`, and `Outcomes & Retrospective` current while implementing.

## Purpose / Big Picture

Give Phase 1 one canonical deep module for declared run identity, execution
attempt metadata, row-aware artifact observations, closeout, and the
pre-consumption validator. Researchers should be able to run the isolated CLI
fixture, mutate or duplicate an artifact, and see a precise rejection before a
mock consumer is allowed to score it. Existing benchmark script imports and
CLIs remain usable through thin adapters; stable pipeline defaults do not
change.

## Progress

- [x] Inspect current contracts and establish the smallest canonical module
- [x] Add `src/provenance` core and compatibility adapters
- [x] Add the end-to-end Phase 1 fixture and focused regressions
- [x] Run validation and review the critical changed regions
- [ ] Record outcomes and project learnings

## Surprises & Discoveries

- Observation: the repository has no `CONTEXT.md`; the domain vocabulary is
  currently in `.planning/PROJECT.md` and the Phase 1 context.
  Evidence: `rg --files -g 'CONTEXT.md'` returned no file.
- Observation: provenance, artifact writing, candidate lineage, score
  contracts, and completion classification are currently split across five
  benchmark modules.
  Evidence: `benchmark/scripts/investigation_provenance.py`,
  `investigation_artifacts.py`, `investigation_lineage.py`,
  `investigation_contracts.py`, and `pipeline_completion_contract.py`.
- Observation: the worktree contains extensive user changes and untracked
  benchmark evidence.
  Evidence: `git status --short` on 2026-07-30.
- Observation: importing the canonical package from a path-executed legacy
  script requires a repository-root fallback because Python sets `sys.path[0]`
  to `benchmark/scripts`.
  Evidence: direct `python benchmark/scripts/investigation_provenance.py --help`
  first failed with `ModuleNotFoundError: No module named 'src'`; the narrow
  fallback now makes the same command exit 0.

## Decision Log

- Decision: limit this change to Phase 1 run identity, manifest, artifact
  ledger, closeout, and consumer gate.
  Rationale: candidate lifecycle, evaluation joins, HPC recovery, and evidence
  bundles are explicitly deferred to later phases.
  Date/Author: 2026-07-30 / user-confirmed architecture grilling.
- Decision: make `src/provenance` the canonical module location.
  Rationale: production orchestration and benchmark consumers need one import
  surface without making script paths the owner of the contract.
  Date/Author: 2026-07-30 / user-confirmed architecture grilling.
- Decision: keep existing `benchmark/scripts/investigation_*` imports and CLIs
  as thin compatibility adapters.
  Rationale: preserve retained benchmark workflows and reduce migration risk in
  the dirty worktree.
  Date/Author: 2026-07-30 / user-confirmed architecture grilling.
- Decision: enter execution through the reproducible CLI/wrapper first, then
  add an opt-in `prism.py` hook.
  Rationale: prove the seam in isolation before coupling it to orchestration.
  Date/Author: 2026-07-30 / user-confirmed architecture grilling.
- Decision: make the minimum proof an end-to-end temporary fixture covering
  symlinks, missing artifacts, retry linkage, mutation, duplicate rows, and
  mock-consumer rejection.
  Rationale: the interface is the test surface; helper-only tests would miss
  identity/ledger/closeout/gate wiring.
  Date/Author: 2026-07-30 / user-confirmed architecture grilling.

## Outcomes & Retrospective

The Phase 1 evidence semantics now live in `src/provenance/run_evidence.py`.
The existing provenance script delegates canonical JSON and file hashing while
retaining its current CLI and template-preflight behavior. The new fixture
proves contract/attempt identity, symlink hashing, mutation rejection, missing
records, duplicate keys, row mismatch, closeout, and CLI output.

The full repository pytest invocation was attempted but returned no captured
pytest summary in the shell; the focused Phase 1 and compatibility subset is
the reliable validation result. A future phase should decide whether the
remaining template-preflight and benchmark artifact writers should migrate
fully into `src/provenance` or remain format-specific adapters.

## Context and Orientation

The project is a file-oriented Python 3 protein–protein docking pipeline. The
Phase 1 domain terms are run, contract, execution attempt, artifact
observation, row identity, scientific role, closeout, and consumer gate. The
stable production path is NACCESS + TMalign + external Rosetta.

Relevant current code:

- `benchmark/scripts/investigation_provenance.py` contains deterministic JSON,
  hashing, Git, environment redaction, executable/package, Slurm, and template
  preflight helpers.
- `benchmark/scripts/investigation_artifacts.py` writes several TSV views and
  freezes pose artifacts, but imports candidate/score contracts.
- `benchmark/scripts/investigation_lineage.py` owns append-only candidate
  lineage and ranking helpers; its candidate-specific fields must not become
  Phase 1 artifact fields.
- `benchmark/scripts/investigation_contracts.py` owns evaluation and ranking
  foreign-key contracts and remains outside this change.
- `benchmark/scripts/pipeline_completion_contract.py` classifies current run
  directories and must not become the owner of Phase 1 identity.
- `.planning/phases/01-run-identity-and-manifest/01-CONTEXT.md` and
  `docs/adr/0001-contract-and-attempt-identity.md`/
  `0002-row-aware-artifact-ledger-and-consumer-gate.md` define the accepted
  contract/attempt and row-aware ledger decisions.

## Plan of Work

First create a standard-library-only `src/provenance` package with a deep
run-evidence module. It will own canonical serialization, contract/attempt
records, expected-artifact observations, row-key collision checks, symlink
metadata, immutable closeout validation, and the reusable consumer gate. Keep
format-specific writing and legacy command names at the edge.

Then add a compatibility layer so current provenance callers can import the
canonical helpers without changing their CLI behavior. Do not merge candidate
lineage, ranking, evaluator joins, or stage lifecycle into the package. Add a
temporary-directory fixture that exercises the whole path and proves failures
remain explicit. Only after that proof passes should an opt-in pipeline hook be
considered; this plan does not change stable `prism.py` defaults.

## Concrete Steps

1. Working directory `/scratch/rshadi25/GitHub/PRISM-prescript`: inspect the
   current module/test contracts and record any incompatible field names before
   editing.
2. Add `src/provenance/run_evidence.py` and `src/provenance/__init__.py` with
   the Phase 1 core and stable typed records.
3. Update the smallest relevant benchmark script imports to delegate to the
   canonical core while preserving current CLI paths and output formats.
4. Add `tests/test_run_evidence.py` for success, mutation, symlink, missing,
   duplicate-key, row-mismatch, retry, and secret-redaction behavior.
5. Run the focused tests from the repository root with
   `/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q
   tests/test_run_evidence.py tests/test_investigation_provenance.py
   tests/test_investigation_artifacts.py`.
6. Run the broader contract subset only if the focused suite is green, then
   inspect `git diff --stat`, `git diff --check`, and the critical regions.

## Validation and Acceptance

Acceptance is observable when:

- identical declared contract inputs produce identical contract hashes while
  distinct attempts retain distinct `run_id` values;
- a dirty-tree/source/configuration inventory is represented without requiring
  a clean worktree or hashing unrelated generated evidence;
- every expected artifact has an explicit row identity, scientific role,
  relative path, existence/status, size, and SHA256; symlinks retain link
  metadata while hashing target bytes;
- closeout rejects changed bytes, duplicate keys, missing rows, and mismatched
  dataset identities before a mock consumer;
- retries create a new attempt linked to the original and never rewrite its
  ledger;
- sentinel secret values do not appear in serialized provenance;
- existing focused compatibility tests remain green.

## Idempotence and Recovery

All fixtures use `tmp_path` or a fresh `tmp/agent/<run-id>/` directory and write
derived evidence only. Re-running the focused tests is safe. Do not mutate raw
benchmark inputs, validated outputs, project memory, or the existing dirty
artifacts. If an adapter migration fails, preserve the new core module and
restore only the touched compatibility file from the patch after review; never
use `git reset --hard` or broad cleanup.

This is standard-library/local validation work; no Slurm job, external binary,
download, package installation, or compute-node execution is required.

## Artifacts and Notes

- ExecPlan: `docs/exec-plans/20260730-deepen-phase1-evidence-ledger.md`
- Canonical code: `src/provenance/`
- Focused regression: `tests/test_run_evidence.py`
- No temporary outputs should remain outside pytest-managed temporary paths.

## Interfaces and Dependencies

The package must use Python 3.11-compatible standard-library types and
`pathlib.Path`; it must not add dependencies. Public records should be
serializable to canonical JSON and deterministic TSV. Existing benchmark
adapters may continue to expose their current function names and CLI flags.
The package must not import candidate ranking, DockQ, Rosetta, Slurm, or network
clients. Slurm fields are captured as nullable runtime data only when present.
