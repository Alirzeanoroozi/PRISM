# Add an opt-in comparable alignment contract

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`,
`Decision Log`, and `Outcomes & Retrospective` current while implementing.

## Purpose / Big Picture

Provide a reproducible comparison mode in which TMalign, USalign, GTalign, and
MultiProt use the same transformation acceptance gates: matched-residue count,
interface coverage, orientation policy, protocol-filter mode, and clash policy.
Native aligner score contracts remain available by default. The comparison mode
does not reinterpret MultiProt's RMSD-derived proxy as a TM-score.

## Progress

- [x] Inspect the active PRISM checkout, existing plan, source contracts, and focused tests.
- [x] Add the explicit comparison contract to configuration and CLI.
- [x] Apply the contract consistently in transformation and stepwise replay gates.
- [x] Add focused regression tests and documentation.
- [x] Run focused tests, compile checks, and diff/provenance review.

## Surprises & Discoveries

- The active checkout is intentionally dirty with prior user evidence and source changes; only the transformation configuration, transformation gates, CLI, tests, and this plan may be touched.
- Native TMalign/USalign/GTalign score gates are not semantically comparable to the MultiProt RMSD proxy. A common match/coverage mode is therefore safer than forcing a fake common TM-score.

## Decision Log

- Decision: Keep `native` as the default and add `common_match_coverage` as an explicit opt-in mode.
  Rationale: Preserve existing production behavior while enabling a fair threshold-controlled method comparison.
  Date/Author: 2026-09-28 / Codex.
- Decision: In comparable mode use the common count/coverage thresholds and inclusive coverage comparison for every aligner; do not apply `tm_score_threshold`.
  Rationale: MultiProt does not provide a TMalign-compatible TM-score. TM-scores remain reported for post-hoc analysis.
  Date/Author: 2026-09-28 / Codex.

## Outcomes & Retrospective

Implemented and verified. Native behavior remains the default; the new mode is
run-selectable and has no scheduler or production-compute side effect. Focused
tests pass, syntax/whitespace checks pass, and CLI help exposes the new mode.
Full production reruns and independent scientific review remain outstanding.

## Context and Orientation

Primary active checkout: `/scratch/rshadi25/GitHub/PRISM-prescript`.

Relevant files are `src/transformation_config.py`, `src/transformation.py`,
`src/stepwise_analysis.py`, `prism.py`, `tests/test_pipeline_configuration.py`,
and `tests/test_transformation_thresholds.py`.

## Plan of Work

Add a validated `alignment_gate_mode` field with `native` and
`common_match_coverage` values. Resolve it from `PRISM_ALIGNMENT_GATE_MODE`
and expose `--alignment-gate-mode`. In common mode, both native and TM-family
records must satisfy the same minimum match count and size-adjusted coverage,
with the same inclusive comparator; TM-score is recorded but not used as a
pass/fail gate. Include the selected mode in serialized threshold evidence.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, update configuration and CLI while preserving existing defaults. Completed.
2. Update transformation and stepwise replay logic so production/native and comparable decisions are explicit and consistent. Completed.
3. Add tests for CLI/environment resolution, equal gates across TMalign and MultiProt, boundary inclusivity, and unchanged native behavior. Completed.
4. Run focused pytest, compileall, shell/diff checks, and inspect the final diff. Completed.

## Validation and Acceptance

- Default `native` behavior still rejects a TMalign record with TM-score below its threshold and preserves the MultiProt native gate.
- `common_match_coverage` accepts/rejects TMalign and MultiProt identically for equal count/coverage inputs regardless of MultiProt's proxy score.
- Boundary coverage is handled identically and is recorded in the audit threshold dictionary.
- No scheduler command or production computation is run by this change.
- `pytest -q tests/test_pipeline_configuration.py tests/test_transformation_thresholds.py tests/test_transformation.py tests/test_transformation_audit.py tests/test_stepwise_analysis.py tests/test_transformation_filter_mode.py tests/test_prism_cli_parity.py tests/test_alignment_parsing.py`: 50 passed.
- `python -m compileall -q src prism.py`: passed.
- `git diff --check` on touched paths: passed.
- `python prism.py --help`: exposes `--alignment-gate-mode {native,common_match_coverage}`.

## Idempotence and Recovery

The mode is run-selectable and writes no generated results by itself. Rerunning
tests is safe. Rollback is deleting only the new mode/CLI/test hunks and the
plan; unrelated dirty files must remain untouched.

## Artifacts and Notes

No raw data, transformed models, refinement outputs, or scheduler state is
modified. Comparison runs must set `PRISM_ALIGNMENT_GATE_MODE=common_match_coverage`
and record the resolved threshold dictionary in their manifest/audit output.

## Interfaces and Dependencies

Python dataclass configuration, PRISM CLI argument parsing, transformation gate
functions, stepwise gate replay, pytest, and the existing PRISM alignment JSON
contracts. No new runtime dependency is required.
