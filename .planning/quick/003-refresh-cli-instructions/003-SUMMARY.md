# Quick Task 003 Summary

**Task:** Verify the CLI parity update and refresh repository instructions
**Completed:** 2026-07-31

## What was done

Verified local branch `feature/prism-cli-parity`, compatibility commit
`20ebf4c3cc3`, the live `prism.py --help`, parser implementation, backend
wiring, accepted project-memory decision, and parity tests. Added a canonical
CLI reference and synchronized the repository tour and showcase with new
aliases, bare flags, backend controls, and refinement behavior. Per user
instruction, no revert instructions or history-changing commands were added or
executed.

## Files changed

- `docs/PRISM_CLI_REFERENCE.md`: canonical parser/runtime command reference.
- `docs/PRISM_REPOSITORY_GUIDED_TOUR.md`: updated command-to-code map and boolean guidance.
- `docs/PRISM_SHOWCASE_TEST_CASE.md`: linked reference and updated full-run command template.
- `.planning/quick/003-refresh-cli-instructions/003-PLAN.md`: scoped plan.
- `.planning/quick/003-refresh-cli-instructions/003-SUMMARY.md`: completion record.
- `.planning/STATE.md`: working copy updated but intentionally not staged because it contains unrelated prior changes.

## Verification

- Live `prism.py --help`: inspected successfully.
- Parser-only matrix: 21 documented command cases passed.
- `tests/test_prism_cli_parity.py tests/test_prism_cli.py`: 7 passed in 1.92 seconds.
- Relative CLI-reference links resolve.
- Stale boolean guidance search completed.
- `git diff --check`: passed.
- No `git revert` or `git reset` instruction appears in the new reference.

## Important finding

- `--inputs_csv`/`--inputs-csv` is accepted and passed to `pdb_downloader()`,
  but `src/transformation.py:transformer()` still reads module-level
  `PRISM_INPUTS_CSV`/`inputs.csv`. The documentation marks this as partially
  wired and gives a same-absolute-path CLI+environment workaround. This should
  receive a separately routed code fix and integration regression.

## Commits

- `d544dc2e414` — verified CLI compatibility reference.
- `09c70f9b06a` — synchronized tour and showcase.

## Deviations

- `.planning/STATE.md` is excluded from the commit to avoid absorbing unrelated user changes.

## Notes for downstream

- Add an end-to-end parser/main test proving a custom CSV reaches download and transformation before describing the option as fully wired.
- The local feature branch has no configured upstream; local metadata cannot establish remote publication state.
