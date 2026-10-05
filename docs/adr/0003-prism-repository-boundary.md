# ADR-0003: PRISM repository boundary and integration ownership

Date: 2026-09-12
Status: accepted for the current integration cycle

## Decision

`PRISM-prescript` remains the maintained PRISM pipeline and the authority for
provenance, benchmark preparation, scoring, validation, and scientific
comparison evidence.

The separate `PRISM` repository is an experimental integration branch. It may
contribute reviewed optional backends, ranking adapters, orientation-safe
output naming, and CLI compatibility behavior, but it is not merged wholesale
into `PRISM-prescript`.

The legacy tree at
`working_version/Multiprot-new/prism-fiberdock-cli/` remains a reference and
compatibility arm only. It is not the implementation owner for the maintained
pipeline.

## Canonical interfaces

- `PRISM-prescript/prism.py` owns the maintained CLI and pipeline defaults.
- Stable defaults remain NACCESS + TMalign + external Rosetta.
- Prescript refinement remains enabled by default; `--no-refine` is the
  explicit opt-out.
- Ranking remains opt-in and must occur after candidate generation and before
  refinement.
- Run-scoped, orientation-safe paths are required for transformed and refined
  candidates.
- Backend adapters must expose explicit input/output contracts, terminal
  status, subprocess return code, failure reason, and score metadata where
  applicable.
- Immutable declared contract identity is separate from execution attempts,
  artifact observations, and run closure/validation.

## Consequences

- Features from `PRISM` must be ported through the prescript contracts and
  focused regression tests before adoption.
- The reverted modular orchestration refactor is not revived as part of this
  cycle; it requires a separate design and migration plan.
- Benchmark and scientific claims must use prescript manifests, ledgers,
  hashes, and evaluator contracts as their evidence source.
