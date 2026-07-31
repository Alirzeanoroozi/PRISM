---
wave: 1
depends_on: []
files_modified:
  - docs/PRISM_SHOWCASE_TEST_CASE.md
  - docs/PRISM_REPOSITORY_GUIDED_TOUR.md
autonomous: true
single_layer_justified: true
objective: "A presenter can demonstrate one verified PRISM candidate-ranking case with exact evidence, expected observations, and safe execution boundaries."
---

# Quick Task 002: PRISM Showcase Test Case

<objective>
Document the retained `5zngA / 4eylA` and `1a0cCD` two-orientation PRODIGY case as a reproducible presentation fixture. The case must distinguish evidence inspection and selector replay from full pipeline regeneration, and it must not imply that top-1 selection proves docking-quality improvement.
</objective>

## Tasks

<task id="002-01">
<title>Create the verified showcase case</title>
<files>
- docs/PRISM_SHOWCASE_TEST_CASE.md
</files>
<action>
Create a presenter script containing the scientific setup, retained evidence root, preflight checks, expected hashes/scores/statuses, step-by-step inspection commands, optional lightweight selector replay, assertions, Socratic prompts, failure interpretation, cleanup policy, and Slurm-only extension. Ground every expected value in `tmp/agent/20260730-prodigy-ranking-paired-test/summary-corrected.json` and its score records.
</action>
<verify>
Validate the retained summary schema and expected values with a read-only Python assertion; verify both score JSON files, combined PDBs, stdout, and stderr files exist; run the smallest relevant PRODIGY ranking unit tests without invoking the external scorer.
</verify>
<done>
The document lets a presenter show two candidate orientations, explain why the more negative affinity wins, verify that top-1 forwards only orientation `o1`, and state what the result does and does not prove.
</done>
</task>

<task id="002-02">
<title>Link the showcase from the repository tour</title>
<files>
- docs/PRISM_REPOSITORY_GUIDED_TOUR.md
</files>
<action>
Add a concise showcase callout near the tour-selection section that links to the dedicated case and recommends it for the retained-evidence demonstration. Preserve all existing tour content.
</action>
<verify>
Confirm the relative Markdown link resolves to the showcase document and `git diff --check` passes.
</verify>
<done>
The main guided tour points presenters directly to the verified showcase case.
</done>
</task>
