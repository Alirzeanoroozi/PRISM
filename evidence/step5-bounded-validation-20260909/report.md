# Step 5 report — bounded validation

## Verdict

The proposed transformation, USalign, and PRODIGY work is `TESTED` and bounded portions are `VALIDATED`, but the final Step 5 gate is not complete: the fresh deterministic PRISM smoke is `UNKNOWN` because Slurm was unavailable. No canonical project file, raw input, or canonical output was changed.

## Status separation

### IMPLEMENTED

USalign and PRODIGY changes exist only in isolated worktrees:

- USalign: `/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-1e6873e93a4f4082a236f7348218ecd5`
- PRODIGY: `/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-4f0e2099ca254fde8dbd32b71a96e740`

Neither candidate was copied into or merged with the canonical checkout.

### TESTED

- Canonical transformation, CLI, selector, and PRODIGY tests: `22 passed in 2.03s` with `/home/rshadi25/.conda/envs/gtalign_env/bin` on PATH.
- Isolated USalign tests: `20 passed in 1.78s`.
- Isolated PRODIGY tests: `9 passed in 0.44s`.
- Step 5 Slurm script syntax: `bash -n` passed.

An unqualified canonical test invocation initially had one PATH-resolution failure for the command name `prodigy`; the environment-corrected rerun passed 22/22. This dependency is recorded rather than hidden.

### VALIDATED

Retained runtime evidence provides bounded quantitative validation:

- Job `1656966`: USalign wrote a PRISM-compatible success record with 214 matches and TM-score 0.98359.
- Job `1656993`: PRODIGY state sequence `available -> executed`, affinity `-65.827`, one selected candidate.
- Job `1656992`: PRODIGY state sequence `available -> failed` on no contacts; both candidates remained selected.
- Job `1400311`: matched current-tree ranking smoke had 4 transformed pairs; baseline selected/refined 4/2 and ranked top-1 selected/refined 1/1, both return code 0.

These results validate compatibility, state handling, and load reduction within bounded fixtures. They do not validate biological ranking quality, full-panel equivalence, or speedup.

### REVIEWED

- Canonical defaults remain `tmalign`, `rank=false`, `rank_method=baseline`, `top_k=5`, `refine=true`.
- Canonical Git status remains 373 entries with SHA-256 `f63a1ccc230321db5ab801ad1f0d7db3ef1e11f497b804eb546dae842b1a8c2e`.
- No raw data or canonical output changes are attributable to Step 5.
- Pre-existing `git diff --check` failures remain in dirty benchmark/source files; Step 5 did not modify them.

## Source and stage review

The maintained runtime path is `prism.py` -> the direct TMalign/GTalign/
MultiProt provider branch -> `src.transformation.transformer` -> optional
candidate selection -> refinement. `src/structural_aligner.py` is not on that
path. USalign and the PRODIGY adapter remain isolated candidates. The current
MultiProt score is a Kabsch-RMSD proxy, not a common TM-score contract, and the
USalign candidate record still lacks an explicit score-contract field.

The retained no-drop ledger quantifies 76,248 potential sides, 2,997 written
alignment records, 560 transformation attempts, 489 clash rejections, and 71
final passes. It records the 73,251 missing/unwritten remainder without
assigning an unsupported provider-failure cause; the 1,877 residual is not
called structural orphanage. External Rosetta output counts remain incomplete
because the current refiner does not persist per-candidate return codes and
score-gate reasons. A native DockQ quality comparison was not run.

The tool matrix now retains all planned Tier 1–4 arms: 31 rows total, with
three bounded completed rows, eight retained/not-comparable rows, and twenty
scheduler-blocked rows. No unrun combination is silently omitted.

The machine-readable matrix and ledgers are linked from
`evidence/tool-matrix.json`, `evidence/tool-matrix.tsv`,
`evidence/alignment-stage-ledger.json`, `evidence/refiner-stage-ledger.json`,
`evidence/matched-candidate-ledger.json`, and
`evidence/no-drop-ledger.json`. Executable provenance and the restricted
MultiProt probe are in
`evidence/tool-preflight-20260909T021142+0300.json`.

The critical-region review of these evidence changes is recorded in
`evidence/change-review.md`.

The requirement-by-requirement completion audit is recorded in
`evidence/completion-audit.md`; it leaves the overall task incomplete.

The established daily-organizer mapping is
`projects/research/prism_refactoring`; its required `26-09-09-02.md` audit
entry was not written because that checkout is outside the writable boundary
and already contains unrelated dirty files. The synchronization blocker is
recorded in `evidence/daily-organizer-sync-blocker.md`; no commit or push was
performed.

## Scheduler blocker

Fresh job `1657005` was submitted with canonical TMalign, isolated USalign, and isolated PRODIGY arms, all staged under `/scratch/tmp`. It has produced no visible output. Repeated `scontrol ping` reconciliations through `2026-09-09T02:11:42+03:00` reported both `headnode01` and `headnode02` Slurm controllers down, and `sacct` could not contact accounting. The job submission is therefore `UNKNOWN`; it is not counted as a run, pass, or failure, and no duplicate was submitted.

## Residual risks and next action

The remaining required evidence is a fresh current-tree one-pair/one-template `--no-refine` smoke after scheduler recovery, followed by inspection of alignment JSON counts, transformation outputs, candidate selection/state logs, and no-drop accounting. Native DockQ comparison and full template-panel evaluation remain separate scientific work. Preserve all current candidates and evidence until that rerun is reconciled.
