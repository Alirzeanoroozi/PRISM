# Quick Task 002 Summary

**Task:** Add a verified test case for the PRISM repository showcase
**Completed:** 2026-07-31

## What was done

Documented the retained `5zngA / 4eylA` query and `1a0cCD` template as an
8–10-minute showcase. The fixture exposes two transformed orientations, two
successful PRODIGY scores, complete score evidence, and a top-1 decision that
forwards `o1`; it also separates mechanical selection from any biological
quality claim. The main repository tour now links directly to the case.

## Files changed

- `docs/PRISM_SHOWCASE_TEST_CASE.md`: standalone showcase contract and presenter script.
- `docs/PRISM_REPOSITORY_GUIDED_TOUR.md`: link to the showcase.
- `.planning/quick/002-prism-showcase-test-case/002-PLAN.md`: scoped plan.
- `.planning/quick/002-prism-showcase-test-case/002-SUMMARY.md`: completion record.
- `.planning/STATE.md`: working copy updated but intentionally not staged because it contains unrelated prior changes.

## Verification

- Read-only fixture assertion: passed.
- Both candidates have `status=scored` and `return_code=0`.
- Expected affinities: `o1=-65.827`, `o2=-65.274` kcal/mol.
- Expected selection: two candidates before ranking, one (`o1`) after top-1.
- Referenced combined PDB, stdout, and stderr files exist.
- `tests/test_prodigy_ranker.py`: 4 passed in 2.19 seconds.
- Relative guide link resolves.
- `git diff --check`: passed.

## Commits

- `acc9c359320` — standalone verified showcase case.
- `058217cad27` — link from the main repository tour.

## Deviations

- `.planning/STATE.md` is excluded from the commit to avoid absorbing unrelated user changes.

## Notes for downstream

- The showcase intentionally uses retained evidence and does not regenerate the full docking pipeline.
- A clean full-pipeline reproduction requires a separately staged isolated template workspace and Slurm execution contract.
