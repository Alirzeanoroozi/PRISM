# Roadmap

**Project:** PRISM-prescript Pipeline Reliability and Provenance
**Created:** 2026-07-28
**Phases:** 5
**V1 requirements mapped:** 10/10

## Phase Overview

| # | Phase | Goal | Requirements |
|---|-------|------|--------------|
| 1 | Run Identity and Manifest | Define immutable run, source, configuration, and artifact identities before adding more execution evidence | REL-01, REL-04 |
| 2 | Stage and Candidate Contracts | Make stage and candidate outcomes explicit, terminal, and diagnosable | REL-02, REL-03 |
| 3 | Evaluation and Mapping Audit | Ensure scores join to exact models/natives and preserve interface scope and denominators | REL-05, EVAL-01 |
| 4 | Reproducible HPC Operations | Make isolated reruns and bounded Slurm batches recoverable and reproducible | OPS-01, OPS-02 |
| 5 | Evidence Bundle and Regression Gate | Give researchers a clear completion verdict and protect all hardened contracts | DOC-01, TEST-01 |

## Phases

### Phase 1: Run Identity and Manifest

**Goal:** A run has a durable, immutable identity for its declared inputs,
configuration, runtime, tools, resources, and artifacts.

**Requirements:** REL-01, REL-04
**Depends on:** None

**Success criteria:**

- [ ] A manifest records Git revision, raw and normalized selectors, template
  inventory, thresholds, backend choices, environment/tool versions, command,
  output root, and Slurm metadata when present.
- [ ] Every staged source, intermediate, final model, native, and evaluator
  artifact has a stable record containing path, role, size, and SHA256.
- [ ] A fixture can detect a changed artifact or mismatched row identity before
  downstream scoring consumes it.

**Research needed:** Yes — validate the local schema against existing audit,
lineage, and benchmark manifests before consolidating contracts.

### Phase 2: Stage and Candidate Contracts

**Goal:** Every pipeline stage and candidate has an explicit lifecycle outcome
that explains incomplete or failed work.

**Requirements:** REL-02, REL-03
**Depends on:** Phase 1

**Success criteria:**

- [ ] Stage records distinguish started, completed, failed, timeout, partial,
  skipped, and unavailable outcomes with timestamps and return codes.
- [ ] Candidate records distinguish filter rejection, refinement failure,
  score-gate failure, missing output, timeout, and `not_scoreable`.
- [ ] A simulated failed/partial external-tool run remains visible in the
  ledger and cannot be promoted to scientific success by directory presence.
- [ ] Retry behavior appends or supersedes records without erasing the original
  attempt or changing its identity.

**Research needed:** Yes — inspect current refiner, scorer, and collector
failure contracts and preserve compatible fields where possible.

### Phase 3: Evaluation and Mapping Audit

**Goal:** Researchers can trust that evaluation metrics refer to the intended
benchmark row, chain roles, model/native bytes, and interface scope.

**Requirements:** REL-05, EVAL-01
**Depends on:** Phases 1–2

**Success criteria:**

- [ ] Score joins require `dataset_row_id`, model/native hashes, chain mapping,
  template/orientation metadata, and do not use path-only identity.
- [ ] Evaluation output preserves requested receptor–ligand cross-interface
  metrics, GlobalDockQ separately, raw evaluator JSON, evaluator version, and
  explicit `scored`, `no_native_interface`, `not_scoreable`, and failed states.
- [ ] An audit reports the full denominator and identifies every excluded,
  missing, failed, or non-bijective record.
- [ ] Regression fixtures catch swapped chains, duplicate rows, hash mismatch,
  and invalid mapping cases.

**Research needed:** Yes — confirm current DockQ/iRMSD contracts and the
repository’s canonical benchmark role files before changing joins.

### Phase 4: Reproducible HPC Operations

**Goal:** A researcher can rerun a declared case and scale independent tasks
without losing task identity or failure evidence.

**Requirements:** OPS-01, OPS-02
**Depends on:** Phases 1–3

**Success criteria:**

- [ ] An emitted launcher/environment record recreates an isolated smoke run
  with frozen inputs/templates and records the run root.
- [ ] Every Slurm task has a unique directory, command, job/array/task ID,
  resources, start/end state, return code, and output inventory.
- [ ] Bounded array or internal-worker execution respects project CPU/QOS
  constraints and does not share unsafe fixed-name workspaces.
- [ ] A canceled or failed task remains auditable and does not become complete
  solely because sibling tasks succeeded.

**Research needed:** Yes — validate current cluster partitions, QOS behavior,
filesystem semantics, and external-tool execution context on compute nodes.

### Phase 5: Evidence Bundle and Regression Gate

**Goal:** A structural researcher can determine whether a run is complete and
scientifically interpretable from one linked evidence bundle.

**Requirements:** DOC-01, TEST-01
**Depends on:** Phases 1–4

**Success criteria:**

- [ ] A Markdown/TSV/JSON bundle links summary conclusions to manifests,
  ledgers, audit tables, raw outputs, hashes, and failure details.
- [ ] The bundle states a clear completion verdict and separates operational
  completion from biological/scoring interpretation.
- [ ] Focused regression tests cover malformed outputs, timeouts, missing
  interfaces, non-bijective mappings, hash mismatches, and partial refinement.
- [ ] The stable NACCESS + TMalign + external-Rosetta defaults and opt-in
  ranking behavior remain unchanged in the regression suite.

**Research needed:** Yes — test report usability with retained project evidence
and define the smallest stable report schema before adding dashboards.

## Requirement Coverage

| Requirement | Phase |
|-------------|-------|
| REL-01 | 1 |
| REL-02 | 2 |
| REL-03 | 2 |
| REL-04 | 1 |
| REL-05 | 3 |
| EVAL-01 | 3 |
| OPS-01 | 4 |
| OPS-02 | 4 |
| DOC-01 | 5 |
| TEST-01 | 5 |

All 10 v1 requirements map to exactly one phase. V2 requirements are deferred
until the phase contracts and evidence schemas are validated.

---
*Roadmap for: PRISM-prescript Pipeline Reliability and Provenance*
*Created: 2026-07-28*
