# PRISM-prescript Pipeline Reliability and Provenance

## What This Is

PRISM-prescript is a Python pipeline for protein–protein docking: it downloads
and materializes chain-qualified structures, prepares template interfaces,
aligns structures, filters and transforms candidates, refines them, and can
score final poses against native complexes. This project improves the existing
pipeline for structural researchers who need docking results that are
reproducible, interpretable, and explicit about incomplete or failed work.

## Core Value

Researchers must be able to trust what a PRISM result means, how it was
produced, and whether the evidence is complete enough to compare scientifically.

## Requirements

### Validated

- ✓ A current Python 3 orchestration entry point runs input, template,
  surface, alignment, transformation, refinement, and optional comparison
  stages — existing `prism.py` pipeline.
- ✓ Multiple explicit aligner and refiner backends are available, including
  TMalign, GTalign, MultiProt, external Rosetta, PyRosetta, and FiberDock —
  current source adapters.
- ✓ Candidate audits, stage-status records, benchmark manifests, hashes, and
  DockQ/iRMSD evaluation contracts exist as reproducibility foundations —
  current `src/` and `benchmark/scripts/` tooling.
- ✓ Focused unit and contract tests cover pipeline helpers, transformations,
  ranking, scoring, provenance, and benchmark preparation — current `tests/`.

### Active

- [ ] Make every scientifically relevant stage and candidate outcome
  observable enough to distinguish success, failure, timeout, partial output,
  and explicit `not_scoreable` status.
- [ ] Preserve complete provenance for inputs, templates, commands, runtime,
  scheduler resources, transformations, refinements, scores, mappings, and
  output hashes in reproducible isolated runs.
- [ ] Provide researcher-facing validation and comparison artifacts that make
  completeness and failure causes clear before interpreting DockQ, iRMSD, or
  ranking results.
- [ ] Add focused regression coverage for newly hardened output contracts and
  failure paths without changing the stable production workflow implicitly.

### Out of Scope

- Changing the canonical NACCESS + TMalign + external-Rosetta defaults — these
  remain the stable reference path while reliability work is evaluated.
- Enabling experimental candidate ranking or a learned reranker by default —
  ranking utility requires independent evidence beyond mechanical validation.
- Claiming current MultiProt/FiberDock output is historically equivalent to the
  legacy Python 2 pipeline — that requires a separately matched validation.
- Replacing raw benchmark inputs, curated native references, validated outputs,
  or project memory as part of routine implementation — those are preserved
  evidence and remain read-only by default.

## Context

- The repository is a brownfield Python 3/HPC scientific codebase with checked-
  in external tools, Slurm launchers, benchmark collectors, and a retained
  legacy compatibility tree.
- The primary operational environment is the host’s `gtalign_env` interpreter;
  Rosetta 2022.42, GTalign, NACCESS, DockQ, FreeSASA, and optional PyRosetta
  have separate runtime contracts.
- The codebase already records substantial validation history and project
  memory. Current evidence distinguishes mechanical pipeline health from
  biological quality and treats incomplete or observational comparisons
  conservatively.
- Relative paths, module-level state, external binaries, fixed tool filenames,
  and mutable downloaded inputs are known sources of reproducibility risk.

## Constraints

- **Scientific validity**: preserve row-level benchmark identity, chain roles,
  complete mappings, explicit missing/failed states, and separate cross-
  interface from global metrics — otherwise comparisons can be misleading.
- **Compatibility**: retain stable defaults and explicit opt-in experimental
  backends — existing results and downstream workflows depend on them.
- **HPC execution**: heavy alignment, refinement, scoring, and batch loops run
  through Slurm with recorded job/resource metadata — login nodes are not a
  substitute for compute execution.
- **Reproducibility**: use isolated work directories, frozen source/template
  inventories, deterministic manifests where possible, and hashes — live
  downloads and external binaries otherwise make reruns ambiguous.
- **Safety**: preserve the dirty worktree’s user changes and treat raw data,
  validated artifacts, legacy tools, and project memory as protected inputs.

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Improve reliability and provenance before changing scientific defaults | Researchers need trustworthy evidence before quality claims can be interpreted | — Pending |
| Keep the stable NACCESS + TMalign + external-Rosetta path unchanged | It is the documented operational baseline | ✓ Good |
| Keep legacy MultiProt/FiberDock evidence separate from current-pipeline claims | Existing inputs, runtimes, and evaluator contracts are not fully paired | ✓ Good |
| Treat ranking as opt-in until independent quality evidence exists | Current validation proves mechanics and load reduction, not general quality gain | ✓ Good |

---
*Last updated: 2026-07-28 after brownfield project initialization questioning*
