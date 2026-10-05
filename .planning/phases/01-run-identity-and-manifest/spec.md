# Feature Specification: Run Identity and Manifest

**Phase**: 1 — Run Identity and Manifest
**Milestone**: v1.0 — Pipeline Reliability and Provenance

---

## Overview

Define an immutable, inspectable identity and provenance contract for every PRISM run before adding more execution evidence. The phase covers:
- Run manifest (declared scientific contract)
- Execution attempt identity (operational run with runtime context)
- Materialized artifact ledger (row/role/path-keyed TSV)
- Pre-consumption validation gate (artifact hash + row identity verification)

---

## User Stories (Priority Order)

### US1: Declared Contract Identity
**As a** researcher, **I want** a stable scientific identity for my pipeline declaration **so that** re-running with the same inputs/configuration produces the same contract hash even on different hosts/Slurm allocations.

**Acceptance Criteria**:
- Two runs with identical selectors, templates, thresholds, backends, tool versions → same `declared_contract_hash`
- Changing any declared source/config file → different `declared_contract_hash`
- Run ID, timestamps, host, Slurm fields, artifact observations do NOT affect contract hash
- Schema: `contracts/declared-contract.schema.json`

**Priority**: P1 (Foundational)

---

### US2: Execution Attempt Identity
**As a** researcher, **I want** a readable unique `run_id` for each operational attempt **so that** retries and corrections are distinguishable while preserving the declared contract.

**Acceptance Criteria**:
- Each attempt gets `run_id` pattern: `prism-YYYYMMDD-HHMMSS-PID-HASH8`
- Retry links to parent via `parent_run_id` + `supersedes_reason`
- Runtime context captured: argv, cwd, env allowlist, seeds, packages, executables, Slurm metadata, Git provenance
- Schema: `contracts/run-manifest.schema.json` (execution attempt record)

**Priority**: P1 (Foundational)

---

### US3: Artifact Ledger with Row/Role Identity
**As a** pipeline stage, **I want** to append artifact observations to a TSV ledger keyed by `(dataset_row_id, scientific_role, run_relative_path)` **so that** downstream consumers can unambiguously select scientific outputs.

**Acceptance Criteria**:
- Ledger records: run_id, stage, row_id, role, path, path_kind, size, sha256, target_sha256 (for symlinks), link metadata, materialized_at, produced_by, status
- Symlinks: hash target bytes, retain link target text and resolved path
- Missing/unavailable artifacts get explicit records (not implicit absence)
- Primary key uniqueness enforced by validation gate
- Schema: `contracts/artifact-ledger.schema.json`

**Priority**: P1 (Foundational)

---

### US4: Pre-Consumption Validation Gate
**As a** scoring consumer, **I want** a validation gate that rejects mutated artifacts, duplicate row identities, and missing expected artifacts **before** computing DockQ/iRMSD **so that** results are only computed on verified evidence.

**Acceptance Criteria**:
- Validates every ledger row: current file hash matches recorded hash
- Rejects `status=fail` if any artifact hash mismatch or duplicate primary key
- Returns `status=warn` if artifacts missing/unavailable or rows incomplete
- Scans serialized provenance for secret sentinels
- Output schema: `contracts/validation-gate.schema.json`
- CLI: `python -m src.validation_gate --run-root <dir> --output <file>`

**Priority**: P1 (Foundational)

---

### US5: Pipeline Integration
**As a** pipeline operator, **I want** `prism.py` to initialize run identity and emit artifacts at each stage boundary **so that** the contract is automatic for every run without code changes to stages.

**Acceptance Criteria**:
- `prism.py` entry point creates `RunIdentity`, emits `declared-contract.json`, `run-manifest.json`, starts `artifact_ledger.tsv`
- Stages append to ledger at materialization boundaries (opt-in via `PRISM_STAGE_STATUS_PATH`)
- Closeout re-hashes expected artifacts and writes `validation_gate.json`
- Stable NACCESS+TMalign+external-Rosetta defaults unchanged

**Priority**: P2 (Integration)

---

## Technical Constraints

- Python 3.11, stdlib only for core identity modules (hashlib, json, pathlib, subprocess)
- Reuse `benchmark/scripts/investigation_provenance.py` for hashing, env capture, Git provenance
- Reuse `candidate_audit.py` append-only pattern for ledger writing
- Reuse `collect_pipeline_verification_baseline.py` duplicate rejection logic
- No new dependencies for Phase 1
- Isolated run roots: `tmp/agent/prism-*/`
- Must work on login nodes (orchestration) and Slurm compute (stages)

---

## Out of Scope (Deferred to Later Phases)

- Stage/candidate lifecycle taxonomy (Phase 2)
- Identity-safe score joins, chain mappings, evaluation audits (Phase 3)
- Reproducible Slurm launchers, bounded task recovery (Phase 4)
- Linked evidence bundles, full regression gate (Phase 5)

---

## References

- ADR-0001: `docs/adr/0001-contract-and-attempt-identity.md`
- ADR-0002: `docs/adr/0002-row-aware-artifact-ledger-and-consumer-gate.md`
- Context: `.planning/phases/01-run-identity-and-manifest/01-CONTEXT.md`
- Research: `.planning/phases/01-run-identity-and-manifest/01-RESEARCH.md`
- Data Model: `.planning/phases/01-run-identity-and-manifest/data-model.md`
- Quickstart: `.planning/phases/01-run-identity-and-manifest/quickstart.md`
