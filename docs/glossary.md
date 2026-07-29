# PRISM-prescript Glossary

Terms below are the Phase 1 vocabulary for run identity and artifact provenance.

| Term | Definition | Invariant |
|---|---|---|
| **Declared contract** | The normalized scientific declaration of inputs, selectors, source/configuration inventory, templates, thresholds, backends, and tool identity. | Its canonical hash is independent of run ID, host, timestamps, Slurm placement, status, and evolving artifact observations. |
| **Contract hash** | SHA256 of the canonical serialized declared contract. | Same declaration yields the same hash; any declared source/config/input/tool change yields a different hash. |
| **Execution attempt** | One invocation of the pipeline under a unique readable `run_id`. | Attempts are never overwritten; retries are linked to prior attempts. |
| **Run ID** | Human-readable unique identifier for one execution attempt. | It identifies an attempt, not the scientific declaration. |
| **Raw selector** | The user-supplied input token or row before canonical normalization. | Retained for audit and intent reconstruction. |
| **Normalized selector** | The canonical receptor/ligand/chain representation consumed by pipeline code. | Used for deterministic execution and joins; never replaces the raw selector in provenance. |
| **Dataset row ID** | Explicit durable identity for a benchmark/scored input row. | Required for benchmark/scored runs; synthetic exploratory IDs are labeled non-benchmark. |
| **Source inventory** | Explicit set of source/config/template files whose bytes affect the declared contract. | Includes declared untracked files; unrelated generated worktree clutter is not implicitly included. |
| **Expected artifact** | Scientific artifact declared by a stage/run as required or intentionally unavailable. | Must receive a present, missing, failed, or unavailable ledger record at closeout. |
| **Artifact observation** | One ledger record connecting a row, scientific role, run-relative path, and observed bytes/status. | Identity key is `(dataset_row_id, role, relative_path)`. |
| **Scientific role** | The meaning of an artifact, such as `source_input`, `template_interface`, `alignment`, `transformed_model`, `native`, or `evaluator_output`. | Role is explicit; path names do not define role implicitly. |
| **Path kind** | Whether the observed path is a regular file, symlink, missing path, or another unsupported kind. | Symlink metadata and target bytes are recorded separately. |
| **Ledger digest** | Digest of the canonical complete artifact ledger at run closure. | Detects rewritten ledger rows in addition to changed artifact bytes. |
| **Run closure** | Final validation of expected artifact inventory, row identities, artifact bytes, and ledger digest. | A consumer cannot treat a manifest-aware run as scoreable without a passing closure. |
| **Legacy-unverified** | Explicit consumer mode for historical roots without the Phase 1 manifest contract. | Never silently inferred; output remains labeled unverified. |
| **Operational retry** | A new attempt with the same declared contract. | New `run_id`, same contract hash, explicit parent link. |
| **Corrected retry** | A new attempt whose declared source/config/input/tool contract changed. | New contract hash plus parent/supersedes link. |

## Domain model

```text
Declared contract ──has──> Contract hash
        │
        └──used by──> Execution attempt ──emits──> Artifact observations
                              │                         │
                              └──closes with──> Ledger digest + validation verdict
```

## Naming rule

Use `contract_hash` for declared identity and `run_id` for operational identity. Do not use `run_id` as a substitute for dataset-row identity, artifact identity, or contract hash.
