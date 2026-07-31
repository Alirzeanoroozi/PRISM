---
wave: 1
depends_on: []
files_modified:
  - docs/PRISM_CLI_REFERENCE.md
  - docs/PRISM_REPOSITORY_GUIDED_TOUR.md
  - docs/PRISM_SHOWCASE_TEST_CASE.md
autonomous: true
single_layer_justified: true
objective: "Presenters and users have one verified CLI reference for the feature/prism-cli-parity branch, and all tour examples use the current parser contract."
---

# Quick Task 003: Refresh CLI Instructions

<objective>
Verify commit `20ebf4c3cc3` and the current `build_parser()` contract, then document the new aliases, bare boolean flags, backend paths, explicit refinement controls, defaults, caveats, and combined example. Exclude Git revert instructions and do not perform history-changing operations.
</objective>

## Tasks

<task id="003-01">
<title>Create the verified CLI reference</title>
<files>
- docs/PRISM_CLI_REFERENCE.md
</files>
<action>
Create a canonical command reference covering working directory/interpreter, input CSV, template generation/limit, surface backends, aligners and backend controls, refinement modes, ranking, DockQ comparison, a combined command, HPC boundary, and read-only version identification. Use exact aliases/defaults from `prism.py:build_parser()` and behavior wiring from `main()`. Explicitly flag any parser-supported option whose end-to-end propagation is incomplete. Omit all revert/reset instructions.
</action>
<verify>
Run parser-only checks for defaults, aliases, bare/explicit booleans, no-refine, backend paths, PyRosetta initialization options, ranking, comparison, and the combined example. Run `tests/test_prism_cli_parity.py` and relevant CLI tests. Confirm no `git revert` or `git reset` text appears.
</verify>
<done>
Every documented command parses under the validated interpreter and maturity caveats distinguish parser compatibility from end-to-end execution validation.
</done>
</task>

<task id="003-02">
<title>Synchronize the tour and showcase</title>
<files>
- docs/PRISM_REPOSITORY_GUIDED_TOUR.md
- docs/PRISM_SHOWCASE_TEST_CASE.md
</files>
<action>
Link both documents to `PRISM_CLI_REFERENCE.md`; replace the obsolete statement that bare `--compare` is rejected; document bare template/comparison/DockQ flags, input CSV and backend runtime controls, explicit `--no-refine`, and the refinement-on default. Keep scientific maturity and Slurm caveats unchanged.
</action>
<verify>
Confirm all relative links resolve, search the three files for stale boolean guidance, and run `git diff --check`.
</verify>
<done>
The tour, showcase, and reference agree with the current parser and no longer teach superseded CLI behavior.
</done>
</task>
