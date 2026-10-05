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

## VALAR continuation-02 after checkpointed focused-test diagnosis

Recover the exact run, checkpoint, handoff, report, decisions, evidence, and isolated worktree before acting. Continue the same Copilot session and bounded PRISM-prescript repair; do not start a sibling run.

Prior job 1656877 was walltime-checkpointed after static compilation passed and the focused suite reported 28 passed and 6 failed. The worker then established that the failures were caused by synthetic fixture path/naming mismatches, not yet by a proven production regression:
- the pipeline convention is `{template}_{chain}_int.pdb` using the full template ID;
- alignment output convention is `{query}_{template}_{chain}.json`;
- a manual stub repro confirmed fail-closed missing-USalign behavior and produced `1abcA_1defA_A.json` for a successful pair;
- tests were corrected in the isolated worktree but were not rerun before walltime closeout.

Continue from the preserved isolated worktree:
/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-a2ac256d43744ae084719726833b3ef6

Required bounded actions:
1. Rerun the corrected focused tests for USalign, MultiProt runtime diagnostics, and existing CLI parity; diagnose any remaining failures as fixture, implementation, or environment evidence.
2. Run bounded static checks and inspect the isolated diff for imports, paths, JSON schema, error handling, and preservation of default TMalign/NACCESS/external-Rosetta behavior. Do not run heavy pipelines.
3. Keep USalign fail-closed and explicit: do not download/build/run a live unavailable executable and do not silently fall back to another aligner. Synthetic stub/parser tests are not live-tool validation.
4. Keep MultiProt Seccomp/32-bit incompatibility explicit and fail-closed without setting PRISM_MULTIPROT_FORCE; do not claim live compatibility.
5. Write durable report.md, handoff.md, checkpoint.md, decisions.jsonl, and evidence artifacts with OBSERVATION/EVIDENCE/INFERENCE/HYPOTHESIS/UNKNOWN classifications and separate IMPLEMENTED/TESTED/VALIDATED/REVIEWED states. Record exact commands/results and remaining uncertainty.
6. Preserve all existing project changes and never touch the canonical checkout, reset, delete, merge, promote, or push.

Do not claim completion from model confidence, process exit 0, or Slurm COMPLETED alone. The run is complete only when the focused evidence and durable closeout satisfy the bounded objective.