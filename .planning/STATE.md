# Project State

## Project Reference

See: `.planning/PROJECT.md` (updated 2026-07-28)

**Core value:** Researchers must be able to trust what a PRISM result means, how it was produced, and whether the evidence is complete enough to compare scientifically.
**Current focus:** Phase 1 — Run Identity and Manifest, with stage/backend contract hardening

## Current Position

Phase: 1 of 5 (Run Identity and Manifest)
Plan: Cross-repository boundary and contract hardening
Status: In progress; software validation complete, scientific/backend gates remain open
Last activity: 2026-09-17 — DockQ contract correction, standalone scorer and
Slurm CPU propagation fixes, and independent two-repository validation

Progress: [██████░░░░] 60%

## Performance Metrics

**Velocity:**
- Total plans completed: 1
- Average duration: ~30 min
- Total execution time: ~3 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1 — Run Identity and Manifest | 1 | ~3h | ~3h |
| 2 — Stage and Candidate Contracts | 0 | 0 | - |
| 3 — Evaluation and Mapping Audit | 0 | 0 | - |
| 4 — Reproducible HPC Operations | 0 | 0 | - |
| 5 — Evidence Bundle and Regression Gate | 0 | 0 | - |

**Recent Trend:**
- Last 5 plans: 1 (Phase 1)
- Trend: Active hardening; broad validation remains incomplete

## Accumulated Context

### Decisions

Decisions are logged in `PROJECT.md` Key Decisions.

- Preserve stable NACCESS + TMalign + external-Rosetta defaults.
- Keep legacy MultiProt/FiberDock evidence separate from current-pipeline claims.
- Make manifest, status, identity, evaluation, and evidence contracts precede workflow-engine migration or learned ranking.

### Pending Todos

- Reconcile the verified repository-local DockQ runtime (Python 3.9.23) with
  the `environment.yaml` Python 3.11.13 recipe before batch scoring.
- Regenerate corrected benchmark scores and EDA from the canonical staged
  scorer; existing transformed-run CSV/JSON remains diagnostic evidence.
- Explain CPU/GPU GTalign raw-output and parameter differences before combining arms.
- Complete MultiProt true-TM, external-Rosetta observability, and orientation/threshold gates.
- Make the orientation notebook provenance inventory bounded and rerun its
  post-processing after the fix.
- Resolve the six evidence-CSV whitespace diffs and one source trailing-space diff
  only after the dirty worktree classification is reviewed; do not normalize
  benchmark evidence silently.

### Blockers/Concerns

- Existing worktree contains extensive user changes and staged benchmark files; implementation must remain path-scoped and non-destructive.
- External binaries, shared filesystems, Slurm QOS, and live structure downloads require explicit runtime/provenance validation.
- Slurm evidence from 2026-09-12 is run-scoped under tmp/agent/: stable
  TMalign job 1658925, GTalign CPU 1658928, GTalign GPU 1658929, and
  optional backend job 1658926, plus explicit DockQ failure replay 1658930.
- Orientation notebook job 1658931 completed the three no-refinement arms but
  was canceled during broad asset-hash post-processing; partial evidence is
  preserved under tmp/agent/20260912-orientation-study/.
- Full software suite: 344 passed, 6 skipped; one SciPy/NumPy warning.
- Runtime manifest and bounded run-identity smoke passed under
  tmp/agent/20260912-run-identity-smoke/.
- Repository-local DockQ replay job 1658942 passed in
  tmp/agent/20260912-dockq-repo-env/; raw and adapter scores agree, but this
  remains one-pair evaluator evidence rather than a benchmark result.
- Retained BM55 GTalign raw DockQ audit over 14,470 JSON records passed with
  zero parser errors and zero corrected values outside [0,1]. The old
  `best_dockq` sum reached 13.0026 for a representative multichain record;
  corrected complete-complex scoring uses `GlobalDockQ` 0.6843489.
- Focused correction suites passed: 45 `PRISM-prescript` tests and 29 PRISM
  runtime/helper tests; repository-wide pytest collection remains invalid due
  pre-existing duplicate-module and package-layout conflicts.
- After the final correction pass, focused suites passed: 57
  `PRISM-prescript` tests and 50 PRISM tests; a separate PRISM scorer/Slurm
  subset passed 16 tests before the final no-align regressions. Changed Python
  and Slurm syntax checks passed.
- PRISM standalone benchmark scoring now reads raw JSON, uses bounded
  `GlobalDockQ`, records `valid_unscored`, hashes/provenance, and passes an
  explicit DockQ CPU count. Notebook scoring wrappers propagate that count.
- Low-level DockQ helpers default to one CPU when no caller override is
  supplied; final focused tests, compilation, Slurm syntax, whitespace, and
  pinned-help checks passed.
- Canonical pairwise recovery rows now use explicit `scored_cross_only` status
  when complete DockQ mapping crashes; complete-complex `GlobalDockQ` remains
  unavailable on those rows.
- Independent Terra validation returned `PASS WITH CAVEATS`; it found no
  defect in the final scope/status claims and confirmed the production rerun
  and corrected BM55/EDA regeneration remain intentionally pending.
- Graphify refresh was attempted after the code change and interrupted after
  two minutes while traversing the large dirty/evidence tree; no benchmark
  artifacts were modified.

## Session Continuity

Last session: 2026-09-17
Stopped at: DockQ correction plan post-Green validation checkpoint; corrected
score regeneration and EDA remain open
Resume file: `docs/exec-plans/20260917-dockq-contract-correction.md`

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 001 | PRISM repository guided tour | 2026-07-31 | `86bc911ed9e` | `.planning/quick/001-prism-repository-guided-tour/` |
| 002 | PRISM showcase test case | 2026-07-31 | `058217cad27` | `.planning/quick/002-prism-showcase-test-case/` |
| 003 | Refresh CLI instructions | 2026-07-31 | `09c70f9b06a` | `.planning/quick/003-refresh-cli-instructions/` |
| 004 | Current-pipeline function inventory | 2026-08-01 | pending | `.planning/quick/004-current-pipeline-function-inventory/` |
