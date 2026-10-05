## VALAR ACTIVE RUN CONTEXT

Run ID: run-a2ac256d43744ae084719726833b3ef6
Run directory: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6
Project evidence namespace: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6
Project report: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6/report.md
Project logs: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6/logs
Evidence index: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6/evidence/index.json
Manifest: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6/manifest.json
Checkpoint: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6/checkpoint.md
Handoff: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6/handoff.md
Recover these exact active-run files first. Do not scan or resume a sibling run directory.

## VALAR continuation after walltime checkpoint

Recover this exact run, its initial prompt, checkpoint, handoff, and isolated worktree before acting:
- run directory: /home/rshadi25/valar-agent-framework/evidence/prism-prescript/run-a2ac256d43744ae084719726833b3ef6
- isolated worktree: /scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-a2ac256d43744ae084719726833b3ef6
- prior Slurm attempt: 1656846, walltime checkpoint after partial implementation

Continue the same bounded PRISM-prescript repair. Preserve all existing edits; do not reset, discard, overwrite, merge, push, or touch the canonical checkout.

The previous attempt only partially modified src/alignment.py and src/alignment_multiprot.py and created src/alignment_usalign.py. Finish and verify the actual bounded objective:
1. Wire USalign into the real prism.py CLI/stage flow, including an explicit executable/path contract and fail-closed behavior when USalign is unavailable. Do not download or build USalign and do not silently fall back to another aligner.
2. Keep the MultiProt Seccomp/32-bit diagnosis explicit and fail-closed without PRISM_MULTIPROT_FORCE; ensure diagnostics are reproducible and do not claim runtime compatibility.
3. Add or repair the smallest focused deterministic tests for parser/runtime/CLI behavior. Run those tests and static checks in the isolated worktree. Do not run heavy pipelines or broad benchmark jobs.
4. Inspect the resulting diff for import/path/schema errors and preserve stable TMalign/NACCESS/external-Rosetta behavior.
5. Write durable report.md, handoff.md, checkpoint.md, decisions.jsonl, and evidence with OBSERVATION/EVIDENCE/INFERENCE/HYPOTHESIS/UNKNOWN classifications and separate IMPLEMENTED/TESTED/VALIDATED/REVIEWED statuses. If live USalign or MultiProt cannot be executed, record that limitation rather than claiming validation.

The success condition is evidence-backed completion of the two bounded repairs, not merely a process exit or Slurm completion.