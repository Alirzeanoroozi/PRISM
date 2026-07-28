# Phase 1 Discussion Log: Run Identity and Manifest

**Date:** 2026-07-29
**Mode:** standard
**Outcome:** Context ready for planning

This log records the options presented and the user's selected choices. The
recommended option was shown first for each decision, but the selected option
was user-confirmed.

## Gray areas selected

The user selected **“All four (Recommended)”**, covering identity and schema,
artifact hashing, enforcement, and runtime/configuration capture.

## Run identity and schema

### Identity model

Options:
- **“Run ID + manifest hash (Recommended)”** — readable unique attempt ID plus
  canonical manifest hash for immutable declared identity.
- “Content-addressed only” — identical declared reruns share identity.
- “Timestamp/PID only” — retain the existing operational directory identity.

User choice: **“Run ID + manifest hash (Recommended)”**.

Rationale captured: preserve readable per-attempt directories while separating
operational attempt identity from the stable scientific contract identity.

### Manifest format

Options:
- **“Canonical JSON + TSV ledger (Recommended)”** — JSON for run/config
  identity and TSV for inspectable row-oriented artifact records.
- “Single canonical JSON” — nest metadata and artifacts in one file.
- “TSV-first package” — use flat TSV files with minimal JSON.

User choice: **“Canonical JSON + TSV ledger (Recommended)”**.

Rationale captured: retain machine-stable JSON and human/CLI-friendly rows.

### Git state

Options:
- **“Allow, record exact state (Recommended)”** — capture HEAD, porcelain
  status, and diff hashes/content references without blocking dirty worktrees.
- “Require clean tree” — block runs until changes are committed or removed.
- “Record HEAD only” — permit dirty trees but identify only by commit.

User choice: **“Allow, record exact state (Recommended)”**.

Rationale captured: the current worktree is intentionally dirty with existing
user changes, so provenance must distinguish those changes without destructive
cleanup or a clean-tree requirement.

### Selectors

Options:
- **“Preserve raw + normalized (Recommended)”** — retain the declared selector
  and canonical row/chain form consumed by execution.
- “Normalized only” — store only the execution form.
- “Raw only” — leave normalization implicit in code/logs.

User choice: **“Preserve raw + normalized (Recommended)”**.

Rationale captured: preserve both user intent and the exact execution input.

## Artifact ledger

### Artifact coverage

Options:
- **“All materialized scientific artifacts (Recommended)”** — staged inputs,
  templates, intermediates, models, natives, evaluator artifacts, and explicit
  missing/unavailable records.
- “Declared outputs only” — hash only launcher/stage-declared files.
- “Final models and scores only” — defer source/intermediate provenance.

User choice: **“All materialized scientific artifacts (Recommended)”**.

Rationale captured: downstream scientific interpretation needs provenance for
the complete materialized chain, not just final files.

### Symlink handling

Options:
- **“Hash target bytes + record link metadata (Recommended)”** — preserve
  logical path, target, resolved path, size, and target digest.
- “Hash link text only” — treat the target string as content.
- “Reject symlinks” — require copied regular files.

User choice: **“Hash target bytes + record link metadata (Recommended)”**.

Rationale captured: preserve the staging view while making consumed bytes
verifiable.

### Ledger timing

Options:
- **“Append during stages + finalize (Recommended)”** — record artifacts as
  they appear and re-hash/validate at terminal closeout.
- “Final snapshot only” — hash after completion.
- “Start and final snapshots” — capture declared inputs and final outputs only.

User choice: **“Append during stages + finalize (Recommended)”**.

Rationale captured: failed or timed-out work must retain partial evidence while
still receiving a final completeness check.

### Artifact identity

Options:
- **“Row + role + relative path (Recommended)”** — durable row identity,
  scientific role, and run-relative path; reject path-only joins.
- “Relative path only” — use run-relative path as the key.
- “Content hash only” — deduplicate identical bytes across roles.

User choice: **“Row + role + relative path (Recommended)”**.

Rationale captured: a path or byte hash alone cannot prove which benchmark row
or scientific role consumed an artifact.

## Enforcement and retry behavior

### Hash mismatch

Options:
- **“Fail closed with audit evidence (Recommended)”** — block consumption,
  emit mismatch details, and retain the original ledger.
- “Warn and continue” — score with a warning.
- “Re-hash automatically” — replace the expected digest.

User choice: **“Fail closed with audit evidence (Recommended)”**.

Rationale captured: changed bytes must not silently enter scientific scoring.

### Row mismatch

Options:
- **“Reject before consumption (Recommended)”** — fail with explicit
  duplicate/missing/mismatch diagnostics while retaining observed records.
- “Mark audit-invalid and continue” — generate partial downstream output.
- “Silently normalize” — resolve aliases or duplicates automatically.

User choice: **“Reject before consumption (Recommended)”**.

Rationale captured: identity defects are safer to diagnose before scoring than
to repair implicitly downstream.

### Retry policy

Options:
- **“New immutable attempt, linked (Recommended)”** — preserve the original
  and link a corrected attempt with parent/supersedes provenance.
- “In-place repair” — update the original ledger.
- “New unrelated run” — restart without an explicit relation.

User choice: **“New immutable attempt, linked (Recommended)”**.

Rationale captured: retry history is part of reproducibility and must not erase
the original failure evidence.

### Gate location

Options:
- **“Reusable library + CLI gate (Recommended)”** — one validator for code and
  a reproducible command for scoring/collectors.
- “Library only” — keep validation internal to Python callers.
- “CLI only” — use shell-level checks without a shared programmatic contract.

User choice: **“Reusable library + CLI gate (Recommended)”**.

Rationale captured: both pipeline code and independent downstream consumers
need the same pre-consumption contract.

## Runtime and configuration capture

### Environment capture

Options:
- **“Allowlist + versions + redaction (Recommended)”** — relevant environment,
  package/tool versions, seeds, and Slurm data without unrelated secrets.
- “Full environment” — capture every variable.
- “Minimal environment” — record only interpreter and tool paths.

User choice: **“Allowlist + versions + redaction (Recommended)”**.

Rationale captured: retain reproducibility-critical information while honoring
the project's secret-safety boundary.

### Command record

Options:
- **“Structured argv + cwd/env + launcher hash (Recommended)”** — replayable
  argument tokens, working directory, effective relevant environment, and
  launcher identity.
- “Shell command string” — one copy-pastable shell line.
- “Launcher path only” — rely on a script path and repository state.

User choice: **“Structured argv + cwd/env + launcher hash (Recommended)”**.

Rationale captured: structured records avoid shell quoting ambiguity while
retaining enough information to reproduce the invocation.

### Execution context

Options:
- **“Explicit context with nullable Slurm fields (Recommended)”** — one schema
  for local and Slurm runs with scheduler fields populated when present.
- “Slurm-only schema” — require scheduler metadata.
- “Separate local and Slurm manifests” — use different schemas by mode.

User choice: **“Explicit context with nullable Slurm fields (Recommended)”**.

Rationale captured: local smoke checks and Slurm runs should remain comparable
without pretending local execution has scheduler metadata.

### Tool identity

Options:
- **“Resolved path + version/probe + hash when possible (Recommended)”** —
  record path, reported version/probe, and executable/script digest when
  accessible.
- “Path and version only” — omit executable hashing.
- “Tool names only” — retain portable names without runtime verification.

User choice: **“Resolved path + version/probe + hash when possible (Recommended)”**.

Rationale captured: external binaries and scripts are a material part of the
  scientific runtime and should be fingerprinted when feasible.

## Deferred ideas

- Stage/candidate lifecycle contracts, evaluation joins, HPC recovery, and
  evidence bundles remain in later roadmap phases.
- Content-addressed reuse, provenance querying, workflow-engine migration, and
  learned ranking quality remain deferred V2/later concerns.

---
*Phase: 01-run-identity-and-manifest*
*Discussion logged: 2026-07-29*
