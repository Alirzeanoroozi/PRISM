# Project State

## Project Reference

See: `.planning/PROJECT.md` (updated 2026-07-28)

**Core value:** Researchers must be able to trust what a PRISM result means,
how it was produced, and whether the evidence is complete enough to compare
scientifically.
**Current focus:** Phase 1 — Run Identity and Manifest

## Current Position

Phase: 1 of 5 (Run Identity and Manifest)
Plan: 0 of 0 in current phase
Status: Domain model refined; ready to plan
Last activity: 2026-07-29 — Phase 1 grilling decisions and ADRs captured

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity:**
- Total plans completed: 0
- Average duration: 0 min
- Total execution time: 0.0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1 — Run Identity and Manifest | 0 | 0 | - |
| 2 — Stage and Candidate Contracts | 0 | 0 | - |
| 3 — Evaluation and Mapping Audit | 0 | 0 | - |
| 4 — Reproducible HPC Operations | 0 | 0 | - |
| 5 — Evidence Bundle and Regression Gate | 0 | 0 | - |

**Recent Trend:**
- Last 5 plans: none
- Trend: Not established

## Accumulated Context

### Decisions

Decisions are logged in `PROJECT.md` Key Decisions.

- Preserve stable NACCESS + TMalign + external-Rosetta defaults.
- Keep legacy MultiProt/FiberDock evidence separate from current-pipeline
  claims.
- Make manifest, status, identity, evaluation, and evidence contracts precede
  workflow-engine migration or learned ranking.

### Pending Todos

None yet.

### Blockers/Concerns

- Existing worktree contains extensive user changes and staged benchmark files;
  implementation must remain path-scoped and non-destructive.
- External binaries, shared filesystems, Slurm QOS, and live structure downloads
  require explicit runtime/provenance validation.

## Session Continuity

Last session: 2026-07-29
Stopped at: Phase 1 domain model refined; ready for plan-phase 1.
Resume file: None
