---
wave: 1
depends_on: []
files_modified:
  - tools/build_pipeline_function_inventory.py
  - docs/PRISM_FUNCTION_INVENTORY.md
  - docs/PRISM_CLI_REFERENCE.md
  - docs/PRISM_REPOSITORY_GUIDED_TOUR.md
autonomous: true
single_layer_justified: true
objective: "The current PRISM pipeline has a generated, source-checked function/file inventory linked from its presenter and CLI instructions."
---

# Quick Task 004: Current Pipeline Function Inventory

<objective>
Cover the maintained current pipeline surface (`prism.py` plus `src/`) without pretending benchmark and legacy trees are part of the same runtime. Generate the inventory from Python AST data so function and file drift is detectable and the documentation uses the repository's existing script structure.
</objective>

## Task

<task id="004-01">
<title>Generate and link the maintained function inventory</title>
<files>
- tools/build_pipeline_function_inventory.py
- docs/PRISM_FUNCTION_INVENTORY.md
- docs/PRISM_CLI_REFERENCE.md
- docs/PRISM_REPOSITORY_GUIDED_TOUR.md
</files>
<action>
Add an argparse/Path-based standard-library tool that scans `prism.py` and every `src/**/*.py`, records classes, methods, functions, private helpers, signatures, doc summaries, and line numbers, and writes a generated Markdown file with the current file structure. Link the inventory from the CLI reference and repository tour and document how to regenerate it.
</action>
<verify>
Run the generator, `py_compile` it, compare generated module and callable counts against an independent AST scan, run `git diff --check`, and run the focused CLI/ranking tests.
</verify>
<done>
The generated inventory covers all 41 maintained modules and 256 classes/functions, its generator is reusable after source changes, and presenter documentation links to it without expanding scope into benchmark or legacy trees.
</done>
</task>
