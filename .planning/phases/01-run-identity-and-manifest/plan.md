# Implementation Plan: Run Identity and Manifest

**Branch**: `01-run-identity-and-manifest` | **Date**: 2026-07-29 | **Spec**: `.planning/phases/01-run-identity-and-manifest/spec.md`

**Input**: Feature specification from `.planning/phases/01-run-identity-and-manifest/spec.md`

## Summary

Define an immutable, inspectable identity and provenance contract for every PRISM run. Phase 1 establishes the **declared contract** (scientific identity), **execution attempt** (operational identity), **artifact ledger** (row/role/path-keyed TSV), and **validation gate** (pre-consumption integrity check). This wraps the existing pipeline (`prism.py`) without changing stable scientific defaults (NACCESS + TMalign + external Rosetta).

## Technical Context

**Language/Version**: Python 3.11 (host `gtalign_env` interpreter)

**Primary Dependencies**: Biopython, NumPy, pandas, DockQ 2.1.3, optional FreeSASA 2.2.1, optional PyRosetta, external Rosetta 2022.42, TMalign, GTalign, NACCESS, MultiProt, FiberDock

**Storage**: File-based (JSON manifests, TSV artifact ledgers, SHA256 hashes) in isolated `tmp/agent/` run directories

**Testing**: pytest (unit), contract/regression tests in `tests/`

**Target Platform**: Linux HPC (Slurm ai partition), login nodes for orchestration only

**Project Type**: CLI pipeline / scientific workflow orchestration

**Performance Goals**: Run identity/manifest overhead <2% of total pipeline time; artifact hashing parallelizable per-stage

**Constraints**: 
- Preserve dirty worktree user changes (record exact Git revision + diff fingerprint)
- Raw datasets, validated benchmark outputs, legacy tools, project memory = read-only
- Slurm ai QoS: max 8 running / 50 submitted jobs — use internal parallelization within 1 job
- All stages already emit structured outputs; Phase 1 adds unified identity layer

**Scale/Scope**: 1 pipeline entry point (`prism.py`), 6 stage modules, ~14 template assets, benchmark cohort up to ~200 pairs

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Scientific Validity | ✅ PASS | Phase 1 preserves row-level benchmark identity, chain roles, complete mappings, explicit missing/failed states |
| II. Compatibility | ✅ PASS | Stable NACCESS+TMalign+external-Rosetta defaults unchanged; Phase 1 wraps, doesn't replace |
| III. HPC Execution | ✅ PASS | Uses existing Slurm patterns; run identity captured for both local and batch runs |
| IV. Reproducibility | ✅ PASS | Isolated work dirs, frozen inventories, deterministic manifests, SHA256 hashes — all aligned |
| V. Safety | ✅ PASS | Raw data/validated artifacts/legacy tools/project memory protected; dirty worktree preserved |

**No violations** — Phase 1 design fits within existing constraints without new architectural patterns.

## Project Structure

### Documentation (this feature)

```text
.planning/phases/01-run-identity-and-manifest/
├── plan.md                      # This file
├── spec.md                      # User stories (created)
├── research.md                  # Phase 0 output (01-RESEARCH.md)
├── data-model.md                # Phase 1 output
├── quickstart.md                # Phase 1 output
├── contracts/                   # Phase 1 output
│   ├── declared-contract.schema.json
│   ├── run-manifest.schema.json
│   ├── artifact-ledger.schema.json
│   └── validation-gate.schema.json
└── tasks.md                     # Phase 2 output (created)
```

### Source Code (repository root)

```text
# Existing structure (Phase 1 extends, doesn't restructure)
src/
├── run_identity.py          # NEW: run_id, declared contract, execution attempt, context capture
├── artifact_ledger.py       # NEW: append-only TSV ledger writer/reader, symlink handling
├── validation_gate.py       # NEW: pre-consumption hash/identity validator (library + CLI)
├── alignment.py             # EXISTING: TMalign adapter
├── alignment_gtalign.py     # EXISTING: GTalign adapter  
├── alignment_multiprot.py   # EXISTING: MultiProt adapter
├── transformation.py        # EXISTING: candidate filtering/transformation gates
├── refinement.py            # EXISTING: refinement orchestration (as refinement.py)
├── candidate_audit.py       # EXISTING: JSONL audit pattern (reuse for events)
├── evaluation.py            # EXISTING: DockQ/iRMSD scoring (as evaluation.py)
├── ranking.py               # EXISTING: opt-in ranking
└── provenance.py            # EXISTING: investigation_provenance helpers (extend)

benchmark/scripts/
├── investigation_provenance.py  # EXISTING: deterministic hashing, env capture, Git provenance
└── collect_pipeline_verification_baseline.py  # EXISTING: duplicate rejection patterns

prism.py                     # EXISTING: pipeline entry point (add run identity init)

tests/
├── contract/                # TO CREATE: contract tests
├── integration/             # TO CREATE: integration tests  
└── unit/                    # TO CREATE: unit tests (run_identity, artifact_ledger, validation_gate)

tmp/agent/                   # EXISTING: isolated run roots (created per run)
```

**Structure Decision**: Phase 1 adds 3 new modules (`run_identity.py`, `artifact_ledger.py`, `validation_gate.py`) and extends `prism.py` without restructuring existing code. This follows the "Minimal Fix, Surgical Change" principle.

## Complexity Tracking

> **Fill ONLY if Constitution Check has violations that must be justified**

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| (none) | — | — |
