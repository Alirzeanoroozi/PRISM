---
wave: 1
depends_on: []
files_modified:
  - docs/PRISM_REPOSITORY_GUIDED_TOUR.md
autonomous: true
single_layer_justified: true
objective: "A presenter can lead a verified, junior-friendly tour from PRISM commands through pipeline code, external tools, outputs, and current development boundaries."
---

# Quick Task 001: PRISM Repository Guided Tour

<objective>
Create one presentation-ready repository tour grounded in the current CLI, source modules, operational documentation, tests, and project memory. The guide must connect each demonstrated command to its implementation and outputs while separating lightweight inspection from Slurm-only compute and stable behavior from experimental or incomplete capabilities.
</objective>

## Tasks

<task id="001-01">
<title>Build the verified command-to-code tour map</title>
<files>
- docs/PRISM_REPOSITORY_GUIDED_TOUR.md
</files>
<action>
Inspect `prism.py`, the directly imported `src/` stage modules, stable smoke and validation launchers, focused tests, `docs/STABLE_PIPELINE.md`, and Phase 1 provenance files. Create a concise architecture and command map that records command/purpose, CLI flag or environment control, implementing function/module, external tool, principal input/output, execution safety, and maturity status. Use exact repository-relative paths and symbols. Mark facts that are operationally validated separately from experimental, diagnostic, legacy, or ongoing work.
</action>
<verify>
Cross-check every named CLI flag with `prism.py --help` or its parser definition, every symbol with `rg`, and every operational command with the referenced launcher/documentation. Confirm the guide contains no claim that directory presence alone proves pipeline success.
</verify>
<done>
The guide contains a traceable command-to-code map covering the main pipeline, backend selection, ranking/comparison, smoke validation, evidence inspection, and Phase 1 provenance work.
</done>
</task>

<task id="001-02">
<title>Turn the map into an interactive presentation script</title>
<files>
- docs/PRISM_REPOSITORY_GUIDED_TOUR.md
</files>
<action>
Add a timed presenter agenda, prerequisites, safe live-demo sequence, optional precomputed-output demonstrations, Socratic questions, expected observations, transitions, caveats, recovery notes, and a closing learning recap. Include a short and a full route so the same document supports a brief overview or a detailed onboarding session. Do not prescribe heavy pipeline execution on a login node.
</action>
<verify>
Read the guide end-to-end as a presenter and confirm each segment states what to run or show, where its implementation lives, what the audience should notice, and whether it is safe locally, requires Slurm, or should use retained evidence.
</verify>
<done>
A presenter can conduct the session without reconstructing commands or code locations, and a junior participant is prompted to explain the data flow and validation evidence rather than copy commands blindly.
</done>
</task>
