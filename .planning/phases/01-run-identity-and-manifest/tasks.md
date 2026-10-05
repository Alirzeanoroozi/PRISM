# Tasks: Run Identity and Manifest

**Input**: Design documents from `.planning/phases/01-run-identity-and-manifest/`

**Prerequisites**: plan.md, spec.md, data-model.md, contracts/, quickstart.md

**Tests**: Contract tests included per spec.md (US1-US4 are P1 foundational with test criteria)

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Create new module files and test directories

- [ ] T001 Create `src/run_identity.py` module stub with public API signatures
- [ ] T002 Create `src/artifact_ledger.py` module stub with public API signatures
- [ ] T003 Create `src/validation_gate.py` module stub with public API signatures
- [ ] T004 Create test directories: `tests/unit/`, `tests/contract/`, `tests/integration/`
- [ ] T005 [P] Add `__init__.py` exports for new modules in `src/__init__.py`

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core utilities that ALL user stories depend on - MUST complete before any story

- [ ] T006 [P] Implement canonical JSON serialization in `src/run_identity.py:_canonical_json()` (sorted keys, compact separators, UTF-8, reject NaN)
- [ ] T007 [P] Implement SHA256 file hashing in `src/run_identity.py:sha256_file()` (chunked read, reuse `investigation_provenance.py` pattern)
- [ ] T008 [P] Implement environment allowlist + redaction in `src/run_identity.py:_redact_env()`, `_redact_argv()` (reuse `investigation_provenance.py` patterns)
- [ ] T009 [P] Implement Git provenance capture in `src/run_identity.py:capture_git_provenance()` (HEAD, status, diff hash, declared untracked)
- [ ] T010 [P] Implement tool fingerprinting in `src/run_identity.py:probe_executable()` (resolved path, version output, file SHA256)
- [ ] T011 [P] Implement run_id generation in `src/run_identity.py:generate_run_id()` (pattern: `prism-YYYYMMDD-HHMMSS-PID-HASH8`)
- [ ] T012 [P] Implement TSV writer/reader in `src/artifact_ledger.py:ArtifactLedger.append()`, `.read_all()` (tab-separated, header row, append-only)
- [ ] T013 [P] Implement symlink handling in `src/artifact_ledger.py:_resolve_symlink_fields()` (link_target, resolved_path, target_sha256, is_broken_link)
- [ ] T014 [P] Implement primary key validation in `src/artifact_ledger.py:_check_primary_key()` (dataset_row_id, scientific_role, run_relative_path)
- [ ] T015 [P] Implement validation gate core logic in `src/validation_gate.py:ValidationGate.validate()` (artifact hash checks, row identity checks, secret scan)

---

## Phase 3: User Story 1 — Declared Contract Identity (P1) 🎯 MVP

**Goal**: Stable scientific identity for pipeline declaration; identical inputs/config → identical contract hash

**Independent Test**: `python -m pytest tests/contract/test_declared_contract.py -v` + quickstart Scenario 1

### Contract Tests (Write First, Must Fail)

- [ ] T016 [P] [US1] Contract test: declared contract schema validation in `tests/contract/test_declared_contract.py`
- [ ] T017 [P] [US1] Contract test: canonical JSON determinism in `tests/contract/test_declared_contract.py`
- [ ] T018 [P] [US1] Contract test: hash stability across runs in `tests/contract/test_declared_contract.py`

### Implementation

- [ ] T019 [P] [US1] Implement `DeclaredContract` dataclass in `src/run_identity.py` (all fields from schema)
- [ ] T020 [P] [US1] Implement `DeclaredContract.from_cli_args()` factory in `src/run_identity.py`
- [ ] T021 [P] [US1] Implement `DeclaredContract.canonical_json()` and `contract_hash()` in `src/run_identity.py`
- [ ] T022 [US1] Implement `DeclaredContract.write()` → `run_root/declared-contract.json` in `src/run_identity.py`
- [ ] T023 [US1] Implement `TemplateSelector` and `QuerySelector` dataclasses in `src/run_identity.py` (resolve raw→normalized)
- [ ] T024 [US1] Implement `ResourceRequest` dataclass in `src/run_identity.py` (CPU, memory, time, partition, GPU)
- [ ] T025 [US1] Implement `SourceInventory` dataclass in `src/run_identity.py` (git_head, git_diff_hash, declared_untracked)
- [ ] T026 [US1] Implement `ToolFingerprint` dataclass in `src/run_identity.py` (name, version, sha256)

---

## Phase 4: User Story 2 — Execution Attempt Identity (P1)

**Goal**: Readable unique `run_id` per operational attempt; captures runtime context; links retries to parent

**Independent Test**: `python -m pytest tests/contract/test_run_manifest.py -v` + quickstart Scenario 1, 6, 7

### Contract Tests (Write First, Must Fail)

- [ ] T027 [P] [US2] Contract test: execution attempt schema validation in `tests/contract/test_run_manifest.py`
- [ ] T028 [P] [US2] Contract test: run_id format validation in `tests/contract/test_run_manifest.py`
- [ ] T029 [P] [US2] Contract test: parent_run_id linkage in `tests/contract/test_run_manifest.py`

### Implementation

- [ ] T030 [P] [US2] Implement `RunIdentity` dataclass in `src/run_identity.py` (run_id, created_at, host, user, status, parent_run_id, supersedes_reason)
- [ ] T031 [P] [US2] Implement `RuntimeContext` dataclass in `src/run_identity.py` (argv, cwd, env, seeds, packages, executables, slurm, git_provenance)
- [ ] T032 [US2] Implement `RunManifest` dataclass in `src/run_identity.py` (run_identity, declared_contract_hash, runtime_context)
- [ ] T033 [US2] Implement `RunManifest.write()` → `run_root/run-manifest.json` in `src/run_identity.py`
- [ ] T034 [US2] Implement `initialize_run()` entry point in `src/run_identity.py` (creates run_root, writes declared-contract.json, run-manifest.json, initializes artifact_ledger.tsv)
- [ ] T035 [US2] Implement `retry_run()` in `src/run_identity.py` (new run_id, same declared_contract_hash, parent_run_id link)

---

## Phase 5: User Story 3 — Artifact Ledger with Row/Role Identity (P1)

**Goal**: Append-only TSV ledger keyed by (dataset_row_id, scientific_role, run_relative_path); handles symlinks

**Independent Test**: `python -m pytest tests/contract/test_artifact_ledger.py -v` + quickstart Scenario 2

### Contract Tests (Write First, Must Fail)

- [ ] T036 [P] [US3] Contract test: artifact ledger schema validation in `tests/contract/test_artifact_ledger.py`
- [ ] T037 [P] [US3] Contract test: primary key uniqueness enforcement in `tests/contract/test_artifact_ledger.py`
- [ ] T038 [P] [US3] Contract test: symlink fields populated correctly in `tests/contract/test_artifact_ledger.py`

### Implementation

- [ ] T039 [P] [US3] Implement `ArtifactRecord` dataclass in `src/artifact_ledger.py` (all fields from schema)
- [ ] T040 [P] [US3] Implement `ArtifactRecord.record_id` property (primary key concatenation) in `src/artifact_ledger.py`
- [ ] T041 [US3] Implement `ArtifactLedger.append()` with primary key check in `src/artifact_ledger.py`
- [ ] T042 [US3] Implement `ArtifactLedger.read_all()` with type coercion in `src/artifact_ledger.py`
- [ ] T043 [US3] Implement `ArtifactLedger.closeout()` (re-hash expected artifacts, validate current bytes) in `src/artifact_ledger.py`
- [ ] T044 [US3] Implement stage integration helpers: `append_alignment()`, `append_refinement()`, etc. in `src/artifact_ledger.py`

---

## Phase 6: User Story 4 — Pre-Consumption Validation Gate (P1)

**Goal**: CLI + library validator rejecting mutated artifacts, duplicate row identities, missing expected artifacts

**Independent Test**: `python -m pytest tests/contract/test_validation_gate.py -v` + quickstart Scenarios 3, 4, 5

### Contract Tests (Write First, Must Fail)

- [ ] T045 [P] [US4] Contract test: validation gate schema validation in `tests/contract/test_validation_gate.py`
- [ ] T046 [P] [US4] Contract test: artifact hash mismatch → fail in `tests/contract/test_validation_gate.py`
- [ ] T047 [P] [US4] Contract test: duplicate primary key → fail in `tests/contract/test_validation_gate.py`
- [ ] T048 [P] [US4] Contract test: missing artifact → warn in `tests/contract/test_validation_gate.py`
- [ ] T049 [P] [US4] Contract test: secret sentinel detection in `tests/contract/test_validation_gate.py`

### Implementation

- [ ] T050 [P] [US4] Implement `ArtifactCheck` and `RowIdentityCheck` dataclasses in `src/validation_gate.py`
- [ ] T051 [US4] Implement `ValidationGate.validate_artifacts()` (hash comparison per ledger row) in `src/validation_gate.py`
- [ ] T052 [US4] Implement `ValidationGate.validate_row_identities()` (duplicate PK, missing roles) in `src/validation_gate.py`
- [ ] T053 [US4] Implement `ValidationGate.secret_scan()` (sentinel in env/argv/config) in `src/validation_gate.py`
- [ ] T054 [US4] Implement `ValidationGate.validate()` orchestration + `ValidationResult` dataclass in `src/validation_gate.py`
- [ ] T055 [US4] Implement CLI entry point `python -m src.validation_gate` in `src/validation_gate.py` (args: --run-root, --output, --fail-on-warn)
- [ ] T056 [US4] Implement `ValidationResult.write()` → `run_root/validation_gate.json` in `src/validation_gate.py`

---

## Phase 7: User Story 5 — Pipeline Integration (P2)

**Goal**: `prism.py` initializes run identity, emits artifacts at stage boundaries, runs validation at closeout

**Independent Test**: `python -m pytest tests/integration/test_prism_integration.py -v` + quickstart Scenario 7

### Integration Tests (Write First, Must Fail)

- [ ] T057 [P] [US5] Integration test: prism.py creates run identity on startup in `tests/integration/test_prism_integration.py`
- [ ] T058 [P] [US5] Integration test: stage events appended to ledger in `tests/integration/test_prism_integration.py`
- [ ] T059 [P] [US5] Integration test: closeout validation runs automatically in `tests/integration/test_prism_integration.py`

### Implementation

- [ ] T060 [US5] Modify `prism.py:main()` to call `initialize_run()` at entry in `prism.py`
- [ ] T061 [US5] Add `PRISM_STAGE_STATUS_PATH` opt-in hook in `prism.py` (append stage events to ledger)
- [ ] T062 [US5] Add ledger append calls at each stage materialization boundary in `prism.py` (input staging, alignment output, transformation output, refinement output, evaluation output)
- [ ] T063 [US5] Add closeout call to `ValidationGate.validate()` at pipeline end in `prism.py`
- [ ] T064 [US5] Ensure stable defaults unchanged (NACCESS+TMalign+external-Rosetta) in `prism.py`

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, quickstart validation, cleanup

- [ ] T065 [P] Run quickstart.md Scenario 1 (local run identity + manifest)
- [ ] T066 [P] Run quickstart.md Scenario 2 (symlink handling)
- [ ] T067 [P] Run quickstart.md Scenario 3 (changed artifact → fail)
- [ ] T068 [P] Run quickstart.md Scenario 4 (duplicate row identity → fail)
- [ ] T069 [P] Run quickstart.md Scenario 5 (missing artifact → warn)
- [ ] T070 [P] Run quickstart.md Scenario 6 (Slurm capture)
- [ ] T071 [P] Run quickstart.md Scenario 7 (secret redaction)
- [ ] T072 [P] Run contract test suite: `python -m pytest tests/contract/ -v`
- [ ] T073 [P] Run unit test suite: `python -m pytest tests/unit/ -v`
- [ ] T074 Update `docs/STABLE_PIPELINE.md` with Phase 1 usage notes
- [ ] T075 Add module docstrings and type hints to `src/run_identity.py`, `src/artifact_ledger.py`, `src/validation_gate.py`

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No deps - start immediately
- **Foundational (Phase 2)**: Depends on Setup - BLOCKS all user stories
- **User Stories 1-4 (Phases 3-6)**: All depend on Foundational - can run in parallel after Phase 2
- **User Story 5 (Phase 7)**: Depends on US1-US4 - sequential after they complete
- **Polish (Phase 8)**: Depends on all user stories

### User Story Dependencies

- **US1 (Declared Contract)**: Independent after Foundational
- **US2 (Execution Attempt)**: Independent after Foundational (links to US1 via `declared_contract_hash`)
- **US3 (Artifact Ledger)**: Independent after Foundational
- **US4 (Validation Gate)**: Independent after Foundational (reads US2 manifest + US3 ledger)
- **US5 (Pipeline Integration)**: Requires US1-US4 complete

### Within Each User Story

- Contract tests → Models → Factories/Serialization → Integration helpers
- Tests MUST fail before implementation

### Parallel Opportunities

```
After Phase 2 (Foundational) complete:
  ├─ US1: T019-T026 (models) can run in parallel
  ├─ US2: T030-T035 (models) can run in parallel
  ├─ US3: T039-T044 (models + ledger) can run in parallel
  └─ US4: T050-T056 (validation) can run in parallel

All contract tests (T016-T018, T027-T029, T036-T038, T045-T049) can run in parallel
All quickstart scenarios (T065-T071) can run in parallel
```

---

## MVP Scope

**Stop after Phase 6 (US1-US4 complete)**: Each foundational story independently testable
- `declared-contract.json` + `run-manifest.json` emitted
- `artifact_ledger.tsv` appendable with row/role identity
- `validation_gate` CLI rejects mutations/duplicates/missing
- No `prism.py` integration needed for MVP validation

---

## Task Summary

| Phase | Tasks | Stories |
|-------|-------|---------|
| 1: Setup | 5 | — |
| 2: Foundational | 10 | — |
| 3: US1 | 10 | US1 |
| 4: US2 | 9 | US2 |
| 5: US3 | 9 | US3 |
| 6: US4 | 11 | US4 |
| 7: US5 | 6 | US5 |
| 8: Polish | 11 | — |
| **Total** | **71** | **5** |

**Parallelizable (P-marked)**: 33 tasks (46%)
