# Research Summary

**Domain:** Researcher-facing reliability and provenance for protein–protein docking pipelines
**Researched:** 2026-07-28
**Confidence:** HIGH for the core recommendations; MEDIUM for future workflow-engine migration

## Executive Summary

The ecosystem supports reproducible, scalable scientific workflows through
explicit software environments, workflow DAGs, scheduler metadata, reports,
and provenance. For PRISM-prescript, the highest-value move is incremental:
harden the current Python stage contracts and evidence schemas first, then
evaluate Snakemake or Nextflow as an orchestration layer. A new engine cannot
repair ambiguous chain identity, incomplete output semantics, or mismatched
evaluation contracts.

## Recommended Stack

### Primary Technologies

| Technology | Version | Role |
|------------|---------|------|
| Python | 3.11.x | Maintain current `prism.py` and benchmark adapters |
| Biopython/NumPy/pandas | Existing pinned recipe | Structural parsing, coordinates, and durable tables |
| DockQ | 2.1.3 | Scoped native-vs-model evaluation with JSON preservation |
| Slurm | Site-managed | Bounded HPC execution, arrays, dependencies, and resource metadata |
| Git + run manifests | Repository version | Code, configuration, input/output hashes, and evidence provenance |
| Snakemake | 9.23.1 docs line | Future Python-native DAG/provenance prototype, not an immediate rewrite |

### Key Stack Decisions

- **Keep the current Python stack:** it is operationally validated and already
  exposes the needed stage boundaries.
- **Use manifest-first provenance:** NIH and DFG guidance emphasize versioning,
  metadata, documentation, and traceability for research software ([NIH](https://datascience.nih.gov/tools-and-analytics/best-practices-for-sharing-research-software-faq), [DFG](https://www.dfg.de/en/basics-topics/digital-topics/research-software/principles)).
- **Use scoped DockQ artifacts:** the official tool emits interface-level data
  and JSON, so preserve those before deriving summaries ([DockQ](https://github.com/wallnerlab/DockQ)).
- **Treat Snakemake/Nextflow as alternatives, not prerequisites:** both provide
  useful scaling/provenance features ([Snakemake](https://snakemake.readthedocs.io/en/stable/), [Nextflow](https://nextflow.io/)).

## Table Stakes Features

Features that must be in v1 — users expect these by default:

- [ ] Immutable run manifest with code, inputs, templates, parameters, tools,
  environment, and scheduler resources.
- [ ] Per-stage and per-candidate ledger with success, failure, timeout,
  partial, unavailable, and `not_scoreable` states.
- [ ] Row-, chain-, mapping-, and hash-validated artifact identity.
- [ ] Scope-preserving DockQ/iRMSD evaluation audit with explicit denominator.
- [ ] Reproducible rerun command and isolated work directory.
- [ ] Focused regression fixtures for external-tool and evaluator failure paths.
- [ ] Researcher-readable evidence bundle linking summary, machine tables, raw
  outputs, and completion verdict.

## Key Architecture Decisions

### System Shape

Keep the synchronous `prism.py` pipeline and its explicit external-tool
adapters. Add a manifest/provenance layer around stage boundaries, then add
machine-readable completion and evidence reporting. Prototype a workflow engine
only after equivalent contracts exist.

### Critical Boundaries

| Boundary | What It Separates | Why It Matters |
|----------|-------------------|----------------|
| Run manifest / stage adapters | Declared experiment / actual execution | Makes commands, inputs, tools, and thresholds auditable. |
| Stage artifacts / status ledger | Files / meanings and failure states | Prevents surviving files from being mistaken for valid results. |
| Candidate models / evaluator | Generated pose / native comparison | Ensures exact model/native hashes and chain mappings are scored. |
| Cross-interface scores / global summaries | Requested biology / aggregate convenience | Prevents multichain aggregates from hiding the target interface. |
| Current pipeline / legacy compatibility | Maintained behavior / historical evidence | Avoids unsupported equivalence claims. |

### Recommended Build Order

1. Run manifest and durable identity schema.
2. Shared stage/candidate status and subprocess metadata envelope.
3. Artifact inventory and hash ledger.
4. Refiner/evaluator output contracts and mapping audit.
5. Researcher evidence bundle and completion verdict.
6. Bounded Slurm execution improvements.
7. Isolated Snakemake/Nextflow migration prototype, if still justified.

## Top Pitfalls

The most dangerous mistakes for this domain, ranked by severity:

| # | Pitfall | Severity | Prevention |
|---|---------|----------|------------|
| 1 | Scheduler completion treated as scientific completion | CRITICAL | Require terminal stage, output-integrity, evaluator, and audit contracts. |
| 2 | Joining by PDB/path instead of row and chain identity | CRITICAL | Require `dataset_row_id`, selectors, mappings, and model/native hashes. |
| 3 | Unpaired method comparison | CRITICAL | Freeze source, template, thresholds, backend, evaluator, and candidate budget. |
| 4 | Global-only scoring | HIGH | Preserve requested interface rows, mappings, raw JSON, and GlobalDockQ separately. |
| 5 | Opaque external execution | HIGH | Capture executable versions, commands, return codes, stderr, timeouts, and hashes. |

## Primary Recommendation

Prioritize a manifest-first reliability milestone around the existing Python
pipeline: every stage and candidate should produce an explicit, hash-linked,
researcher-readable record before any new aligner, refiner, workflow engine, or
learned ranking claim is introduced. This directly supports the project’s core
value—trustworthy scientific interpretation—while keeping stable defaults and
legacy boundaries intact.

## Confidence Assessment

- **HIGH:** Official documentation supports workflow reports/provenance,
  scheduler task metadata, RCSB structured APIs, DockQ interface/JSON output,
  and research-software traceability requirements.
- **MEDIUM:** The precise best workflow engine for this repository depends on
  shared-filesystem behavior, external binary packaging, and operator needs.
- **LOW:** No new low-confidence ecosystem claim is required for the v1 scope;
  migration decisions should be tested locally rather than inferred from tools’
  marketing capabilities.

## Gaps

- No formal researcher interviews or acceptance thresholds were available;
  feature priority is grounded in the confirmed project goal and repository
  evidence.
- The cluster’s exact Snakemake/Nextflow/container support and filesystem
  persistence behavior were not tested during this research step.
- External binary licensing and redistribution constraints require local
  confirmation before packaging any environment or container.

---
*Research summary for: reproducible protein–protein docking pipelines*
*Researched: 2026-07-28*
*Sources: STACK.md, FEATURES.md, ARCHITECTURE.md, PITFALLS.md*
