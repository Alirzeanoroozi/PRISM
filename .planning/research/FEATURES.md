# Feature Research

**Domain:** Researcher-facing reliability and provenance for protein–protein docking pipelines
**Researched:** 2026-07-28
**Confidence:** HIGH for provenance/research-software expectations; MEDIUM for feature priority because the user group is inferred from the project goal rather than a formal interview

## Table Stakes

Features users expect by default. Missing these makes a scientific result difficult to trust or reproduce.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Run manifest | Researchers need the exact code, parameters, input selectors, tool paths, versions, and resources used | MEDIUM | Write one immutable manifest per isolated run; include Git revision, environment, Slurm IDs, and command lines. NIH guidance explicitly calls for software versions and platform-specific documentation ([NIH](https://datascience.nih.gov/tools-and-analytics/best-practices-for-sharing-research-software-faq)). |
| Stage and candidate status ledger | A surviving output directory is not enough to establish scientific success | HIGH | Record started/completed/failed/skipped/timeout/partial/not-scoreable states, return codes, reasons, and output paths. |
| Input and output identity | Chain roles, template IDs, orientations, and model/native hashes must be auditable | HIGH | Preserve `dataset_row_id`, raw selectors, normalized IDs, mapping direction, and SHA256 values. |
| Complete evaluation report | Researchers need both metrics and the denominator used to compute them | HIGH | Report score scope, interface mappings, missing/failed rows, DockQ/iRMSD versions, and explicit exclusions. DockQ supports interface-level values and JSON output ([DockQ](https://github.com/wallnerlab/DockQ)). |
| Reproducible rerun command | A result should be executable by another project collaborator | MEDIUM | Emit a shell-safe launcher or machine-readable command/environment record; use isolated workspaces and frozen assets. |
| Regression tests for failure contracts | Scientific pipelines fail at tool boundaries, not only in pure functions | MEDIUM | Test timeout, malformed output, missing native interface, non-bijective mapping, hash mismatch, and partial refinement paths. |

## Differentiators

Features that set a research pipeline apart by reducing interpretation risk rather than merely adding more docking methods.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Researcher-readable evidence bundle | A single report explains what completed, what failed, and whether comparisons are valid | HIGH | Generate Markdown/TSV/JSON artifacts with links to raw outputs, hashes, and stage summaries; aligns with DFG emphasis on documentation and verifiability ([DFG](https://www.dfg.de/en/basics-topics/digital-topics/research-software/principles)). |
| Contract-aware cross-method comparison | Prevents invalid claims caused by unequal inputs, tools, or evaluators | HIGH | Freeze source/template/evaluator manifests, compare stage attrition, and label observational vs causal evidence. |
| Content-addressed or hash-aware reuse | Avoids rerunning expensive alignment/refinement when exact inputs and contracts match | HIGH | Start with explicit hash indexes and safe reuse reports; consider workflow-engine persistence only after shared-filesystem testing. |
| Bounded, failure-resilient batch execution | Makes large HPC runs finish with recoverable per-task evidence | HIGH | Use Slurm arrays or one internally parallelized job with bounded workers; Slurm officially exposes task IDs and array concurrency controls ([Slurm](https://slurm.schedmd.com/job_array.html)). |
| Provenance-aware researcher queries | Lets users answer “which models used this input/tool/threshold?” without scanning logs | MEDIUM | Add stable identifiers and tabular indexes before introducing a database or dashboard. |

## Anti-Features

Features that seem attractive but create scientific or maintenance problems.

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|----------------|-------------|
| Automatic “success” from scheduler completion | Convenient batch summaries | A job can exit zero with missing/partial biological outputs; scheduler state is not scientific validity | Require stage contracts, output integrity, and evaluator/audit completion. |
| One global score per multichain complex | Easy ranking table | It can hide the requested receptor–ligand interface and confuse global with cross-interface quality | Preserve interface-scoped scores and GlobalDockQ separately. |
| Silent fallback between aligners/refiners | Improves apparent completion rate | It changes the method being evaluated and makes comparisons irreproducible | Make backend choice explicit and record fallback as a distinct status, or fail closed. |
| Default learned reranking | Promises better candidates | Current evidence does not establish generalization; it can reduce candidate coverage invisibly | Keep deterministic ranking opt-in and require independent held-out evaluation. |
| Dashboard before stable data contracts | Makes results look accessible | A UI can fossilize ambiguous statuses and encourage interpretation of invalid aggregates | Stabilize manifests, ledgers, and report schemas first. |

## Feature Dependencies

```text
Frozen input/tool manifest
    └──requires──> Stage-status ledger
                         └──requires──> Candidate/output identity hashes
                                               └──requires──> Evaluation audit

Stage-status ledger ──enhances──> Evidence bundle
Evaluation audit ──enhances──> Cross-method comparison

Workflow-engine migration ──conflicts──> Unfrozen current pipeline contracts
```

### Dependency Notes

- **Stage status requires frozen inputs and commands:** otherwise a failure
  record cannot identify what was actually attempted.
- **Evaluation audit requires output identity:** labels and scores must join to
  the exact model bytes and benchmark row.
- **Evidence bundles enhance comparison:** they expose denominator, attrition,
  and method-contract differences before researchers make quality claims.
- **Workflow migration conflicts with unfrozen contracts:** changing execution
  infrastructure and scientific semantics at once makes failures ambiguous.

## MVP Definition

### Launch With (v1)

- [ ] Immutable run manifest — establishes the exact execution context.
- [ ] Per-stage and per-candidate status ledger — makes incomplete work explicit.
- [ ] Hash- and row-validated evaluation audit — prevents incorrect score joins.
- [ ] Focused regression fixtures — protects failure and output contracts.
- [ ] Researcher-readable evidence bundle — provides a usable completion verdict.

### Add After Validation (v1.x)

- [ ] Safe content-addressed reuse — add after hash contracts are stable and
  rerun equivalence is demonstrated.
- [ ] Slurm batch dashboard/index — add after status schemas remain stable across
  several benchmark runs.

### Future Consideration (v2+)

- [ ] Snakemake or Nextflow orchestration — evaluate only after the current
  Python path has equivalent contracts and migration fixtures.
- [ ] Learned candidate reranking — requires independent native complexes and
  coverage-aware quality evaluation.
- [ ] Public researcher portal — out of scope until local evidence artifacts are
  stable and access/licensing requirements are understood.

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Run manifest | HIGH | MEDIUM | P1 |
| Stage/candidate status ledger | HIGH | HIGH | P1 |
| Hash-validated evaluation audit | HIGH | HIGH | P1 |
| Evidence bundle | HIGH | MEDIUM | P1 |
| Failure-contract regression tests | HIGH | MEDIUM | P1 |
| Bounded batch recovery | HIGH | HIGH | P1 |
| Content-addressed reuse | MEDIUM | HIGH | P2 |
| Workflow-engine migration | MEDIUM | HIGH | P3 |
| Learned reranking | MEDIUM | HIGH | P3 |

**Priority key:** P1 must have for the reliability milestone; P2 should follow
validated contracts; P3 requires independent evidence or a separate milestone.

---
*Feature research for: researcher-facing reliable docking pipelines*
*Researched: 2026-07-28*
