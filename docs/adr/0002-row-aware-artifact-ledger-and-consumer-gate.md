# ADR-0002: Use a row-aware artifact ledger with explicit consumer gates

**Status:** Accepted
**Date:** 2026-07-29
**Scope:** Phase 1 — Run Identity and Manifest

## Context

PRISM output directories contain scientific artifacts alongside logs, caches, symlinked staging views, and external-tool scratch. Path-only manifests cannot distinguish duplicate benchmark rows or scientific roles, while scanning every file creates unstable and noisy evidence.

Historical output roots predate the Phase 1 contract and must remain preserved. New scoring must not silently consume mutated or row-mismatched bytes, but historical evidence cannot be retroactively treated as contract-backed merely because it has a directory layout.

## Decision

Use a fixed-column TSV artifact ledger keyed by:

```text
(dataset_row_id, scientific_role, run_relative_path)
```

Each expected artifact has an explicit record containing at least path, role, row identity, path kind, link target/resolved path when relevant, existence/status, byte size, and SHA256. Scientific artifacts are selected from an expected inventory declared by the stage/run; discovered extras are recorded as `unclassified` observations rather than promoted to evidence. Logs and caches are not scientific artifacts unless explicitly declared.

Symlinked artifacts hash target bytes while retaining logical link metadata. Ledger rows carry deterministic row digests, and run closure records the canonical whole-ledger digest. The ledger is append-oriented during stages; closeout re-reads the expected inventory and verifies current bytes, duplicate keys, missing records, and row identity.

For scored or benchmark runs, `dataset_row_id` is required in the input manifest. Exploratory runs may derive a deterministic synthetic ID from the normalized raw row, but that identity is labeled non-benchmark and cannot silently satisfy a benchmark scoring contract.

New manifest-aware consumers must pass the artifact/row validator before scoring. Unmanifested historical roots require an explicit `--legacy-unverified` path and remain labeled unverified. There is no implicit directory-presence success.

## Consequences

### Positive

- A score consumer can reject changed bytes, duplicate rows, missing expected artifacts, and role/path collisions before interpretation.
- Scientific evidence remains separate from operational clutter.
- Historical outputs remain available without weakening the new contract.
- Ledger integrity and artifact integrity are independently diagnosable.

### Negative

- Every benchmark/scored input path must carry an explicit durable row ID.
- Stage code or launchers must declare expected artifact inventories.
- The validator and legacy bypass become part of every downstream scoring interface.

## Rejected alternatives

- **Hash every run-root file:** too noisy and unstable because run roots contain logs, caches, symlinks, and tool scratch.
- **Use relative path only:** permits row/role collisions and path-only misjoins.
- **Make validation optional for all runs:** allows new consumers to score mutated or mismatched artifacts by omission.
- **Use an append-only JSONL chain as the only public format:** stronger event semantics but conflicts with the selected inspectable TSV contract; row and whole-ledger digests provide the required first-phase integrity with one public ledger.

## Verification implications

- A fixture must mutate a ledgered artifact and observe a non-zero validator result with a mismatch diagnostic.
- A fixture must introduce duplicate/missing/mismatched row identity and observe rejection before a mock scoring call.
- A symlink fixture must preserve link metadata while hashing target content.
- A legacy root must require an explicit bypass and be labeled unverified.

## Related records

- `.planning/phases/01-run-identity-and-manifest/01-CONTEXT.md`
- `.planning/phases/01-run-identity-and-manifest/01-RESEARCH.md`
- `benchmark/scripts/build_pipeline_verification_baseline.py`
- `benchmark/scripts/collect_pipeline_verification_baseline.py`
