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
