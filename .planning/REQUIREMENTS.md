# Requirements

**Project:** PRISM-prescript Pipeline Reliability and Provenance
**Defined:** 2026-07-28
**Scope:** Reliability milestone for structural researchers using the current Python pipeline

## V1 Requirements

### Run and provenance contracts

- [ ] **REL-01:** A researcher can obtain an immutable run manifest containing
  the Git revision, input rows and selectors, template inventory, thresholds,
  tool versions/paths, environment, command, output root, and Slurm resources.
- [ ] **REL-02:** A researcher can inspect per-stage status records with
  started, completed, failed, timeout, partial, skipped, and unavailable
  outcomes, including timestamps, return codes, and failure reasons.
- [ ] **REL-03:** A researcher can inspect per-candidate outcomes, including
  filter rejection, refinement failure, score-gate failure, missing output,
  timeout, and explicit `not_scoreable` states.
- [ ] **REL-04:** A researcher can trace each staged input, intermediate, final
  model, native structure, and evaluator artifact to a SHA256 hash and durable
  benchmark row identity.
- [ ] **REL-05:** A researcher can verify that model/native score joins preserve
  `dataset_row_id`, raw and normalized chain roles, mappings, template, and
  orientation rather than relying on paths alone.

### Evaluation and evidence

- [ ] **EVAL-01:** A researcher can receive an evaluation audit that preserves
  requested receptor–ligand cross-interface metrics, GlobalDockQ separately,
  evaluator versions, raw JSON, chain mappings, explicit score scope, and the
  denominator used for summary metrics.
- [ ] **DOC-01:** A researcher can read one evidence bundle linking human-readable
  conclusions to machine-readable ledgers, raw outputs, hashes, and a clear
  completion verdict.

### Operations and validation

- [ ] **OPS-01:** A researcher can reproduce a run from an emitted command and
  environment record in an isolated work directory using frozen inputs and
  templates.
- [ ] **OPS-02:** A researcher can run bounded HPC batches with unique task
  directories, recorded Slurm job/array IDs and resources, and recoverable
  per-task exit records.
- [ ] **TEST-01:** The project has regression fixtures for malformed outputs,
  timeouts, missing native interfaces, non-bijective mappings, hash mismatches,
  and partial refinement outputs.

## V2 Requirements

- **REUSE-01:** A researcher can safely reuse content-addressed artifacts when
  input, tool, configuration, and contract hashes match exactly.
- **QUERY-01:** A researcher can query which models and results used a given
  input, tool, threshold, template, or evaluator without scanning raw logs.
- **BATCH-01:** A researcher can inspect a stable batch index or dashboard built
  from the validated status and evidence schemas.
- **WORKFLOW-01:** The project provides an isolated Snakemake or Nextflow
  orchestration prototype with behavior and provenance equivalent to the
  current Python path.
- **RANK-01:** A researcher can use a learned candidate reranker supported by
  independent held-out native-complex evaluation and coverage reporting.

## Out of Scope

- **SCOPE-01:** Changing the canonical NACCESS + TMalign + external-Rosetta
  defaults. These remain the stable reference path during the milestone.
- **SCOPE-02:** Enabling deterministic or learned ranking by default. Ranking
  remains opt-in until independent quality evidence supports a change.
- **SCOPE-03:** Silent fallback between aligners or refiners. Backend choice and
  any recovery path must remain explicit in provenance records.
- **SCOPE-04:** Reducing multichain evaluation to GlobalDockQ alone. Requested
  interface metrics and mapping evidence remain required.
- **SCOPE-05:** Claiming historical current/legacy equivalence without a paired
  validation that freezes inputs, templates, runtime, and evaluator contracts.
- **SCOPE-06:** Replacing raw benchmark inputs, curated native references,
  validated outputs, project memory, or the retained legacy tree as routine
  implementation work.

## Requirement Dependencies

```text
REL-01 run manifest
    └──> REL-02 stage status ──> REL-03 candidate status
              └──> REL-04 artifact hashes ──> REL-05 identity-safe joins
                                                   └──> EVAL-01 evaluation audit
                                                               └──> DOC-01 evidence bundle

REL-01 + REL-02 + REL-03 ──> OPS-02 bounded batch recovery
REL-02 + REL-04 + EVAL-01 ──> TEST-01 contract regression fixtures
```

## Traceability Rules

- Every v1 requirement must have at least one focused test or validation case.
- Scientific completion requires both operational terminal status and artifact/
  evaluator contracts; Slurm completion alone is insufficient.
- Missing, failed, partial, unavailable, and `not_scoreable` states remain
  explicit records and are never silently dropped from audit denominators.
- Any future comparison or ranking work must preserve the stable defaults and
  use frozen, row-identity-safe evidence.

---
*Requirements for: PRISM-prescript Pipeline Reliability and Provenance*
*Defined: 2026-07-28*
