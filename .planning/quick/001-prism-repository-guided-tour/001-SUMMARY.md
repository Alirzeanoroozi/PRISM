# Quick Task 001 Summary

**Task:** Create a presenter's guided tour of PRISM-prescript commands, code locations, tools, and ongoing functionality
**Completed:** 2026-07-31

## What was done

Created an evidence-backed onboarding guide that maps the current CLI and
environment controls to pipeline functions, external tools, inputs, outputs,
execution boundaries, and maturity status. Added short and full presentation
routes, a source-level trace of one receptor–ligand pair, safe demonstration
commands, Socratic checkpoints, recovery prompts, a presenter run-of-show, and
a participant worksheet.

## Files changed

- `docs/PRISM_REPOSITORY_GUIDED_TOUR.md`: presentation-ready repository tour.
- `.planning/quick/001-prism-repository-guided-tour/001-PLAN.md`: scoped quick-task plan.
- `.planning/quick/001-prism-repository-guided-tour/001-SUMMARY.md`: execution and verification record.
- `.planning/STATE.md`: working copy updated with quick-task status, but intentionally not staged because it contained pre-existing user changes.

## Verification

- `prism.py --help`: passed; current flags and explicit boolean-value behavior confirmed.
- All source, launcher, and test paths cited by the guide exist.
- Named stage functions were cross-checked with `rg`.
- `git diff --check`: passed for the guide and plan.
- Focused pytest command: **30 passed, 1 failed**.
- Existing failure: `tests/test_transformation_thresholds.py::test_multiprot_uses_native_gate_not_tmalign_proxy` expects a MultiProt record with `tm_score = 0.0` to pass, while current `src/transformation.py:alignment_score_passes()` requires true TM-score `>= 0.3`. The guide records this source/test drift; no pipeline code or test was changed.

## Commits

- `ec9f571d013` — verified command-to-code tour.
- `86bc911ed9e` — interactive presenter run-of-show and participant worksheet.

## Deviations

- `.planning/STATE.md` is not included in the final quick-task commit because staging it would include unrelated pre-existing user changes. Its working copy retains the quick-task update.
- The plan's two aspects were initially authored together in the guide; the second atomic commit adds the dedicated run-of-show and worksheet as the independently reviewable presentation layer.

## Notes for downstream

- Resolve the MultiProt gate source/test disagreement before presenting the focused suite as green.
- Reconcile the stale DockQ interpreter command in `docs/STABLE_PIPELINE.md` in a separately routed task.
- Phase 1 provenance remains ongoing: the core exists, but fail-closed CLI behavior and complete `prism.py` integration are not yet accepted.
