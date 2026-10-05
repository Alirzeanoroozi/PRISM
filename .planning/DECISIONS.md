# Project Decisions

## DEC-001: Separate declared contract identity from execution attempt

**Date:** 2026-07-29
**Type:** architecture
**Status:** accepted

Use a canonical contract hash for declared source/tool/config/input identity and a separate readable `run_id` for each operational attempt. Runtime facts, status, and evolving artifact observations are linked provenance, not contract identity.

**ADR:** `docs/adr/0001-contract-and-attempt-identity.md`

## DEC-002: Row-aware artifact ledger and explicit consumer gate

**Date:** 2026-07-29
**Type:** architecture
**Status:** accepted

Use an expected-inventory-driven TSV ledger keyed by `(dataset_row_id, scientific_role, run_relative_path)`, with target-byte hashes plus symlink metadata, per-row/whole-ledger digests, mandatory validation for new manifest-aware runs, and an explicit `--legacy-unverified` bypass for historical roots.

**ADR:** `docs/adr/0002-row-aware-artifact-ledger-and-consumer-gate.md`

## DEC-003: Use GlobalDockQ for complete mappings

**Date:** 2026-09-17
**Type:** scientific evaluation
**Status:** accepted

For a complete multichain mapping, report DockQ's bounded `GlobalDockQ` as
the scalar quality metric. Preserve `best_dockq` only as a diagnostic sum of
interface scores; never expose that sum as an unqualified DockQ value. Do not
promote the first interface's fnat/iRMSD/LRMSD fields to a multichain result.

**Evidence:** DockQ 2.1.3 raw replay over 14,470 retained JSON records; focused
parser tests and the representative correction from 13.0026 to 0.6843489.

## DEC-004: Separate complete-complex and requested cross-interface scores

**Date:** 2026-09-17
**Type:** scientific evaluation
**Status:** accepted

Keep complete-complex `GlobalDockQ`, internal-interface diagnostics, and the
requested receptor-ligand cross-interface scores in separate fields and
separate score scopes. Cross-interface failures and missing native interfaces
remain explicit. Every DockQ invocation records the mapping, raw JSON path and
hash, exact argv, and explicit CPU count.

**Consequence:** Existing benchmark artifacts are preserved as diagnostic
evidence. Corrected BM55 scores and downstream EDA must be regenerated into a
new output root using the canonical staged scorer before scientific comparison.

## DEC-005: Successful execution is distinct from a usable complete score

**Date:** 2026-09-17
**Type:** evaluation provenance
**Status:** accepted

A zero exit code and parseable DockQ JSON do not by themselves establish a
complete-complex quality score. If a valid document has no usable bounded
`GlobalDockQ`, record `valid_unscored` with the reason and retain the raw JSON;
do not count the row as scored or fall back to human-readable `--short`
output.

**Consequence:** EDA and cross-model summaries must filter on score status and
score scope, while valid-unscored and execution-failure rows remain in audit
tables.

## DEC-006: Low-level DockQ helpers default to one CPU

**Date:** 2026-09-17
**Type:** resource provenance
**Status:** accepted

Low-level DockQ parser helpers pass `--n_cpu 1` when callers omit a CPU
setting. Slurm-aware callers may override this with the actual allocated
count, which is recorded in row provenance.

**Rationale:** DockQ's internal default can request more parallelism than a
job allocation. A one-CPU default makes direct and legacy helper calls safe
without preventing an explicitly resource-matched batch invocation.
