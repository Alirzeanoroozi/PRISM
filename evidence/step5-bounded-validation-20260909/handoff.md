# Step 5 handoff — bounded validation

## Current disposition

`PARTIAL_BLOCKED_BY_SCHEDULER`; keep the goal active. Deterministic tests and prior bounded runtime evidence pass, but the fresh Step 5 pipeline smoke has not executed.

## Evidence package

- [final-validation.json](evidence/final-validation.json)
- [error-inventory.json](evidence/error-inventory.json)
- [current-pipeline-map.md](evidence/current-pipeline-map.md)
- [report.md](report.md)
- [checkpoint.md](checkpoint.md)
- [decisions.jsonl](decisions.jsonl)
- Step 5 submission script: `slurm/step5_validation.sbatch`
- Step 5 submitted job: `1657005`

## Verified now

- Canonical focused tests: 22/22.
- Isolated USalign tests: 20/20.
- Isolated PRODIGY tests: 9/9.
- Retained Slurm probes: USalign success, PRODIGY success, PRODIGY failure preservation, and current-tree ranking load reduction. Scheduler rechecks through `2026-09-09T02:51:38+03:00` still report both controllers down; the latest accounting query also failed. The newest check is `evidence/scheduler-recheck-20260909T025138+0300.md`.
- Canonical status count/hash unchanged: 373 / `f63a1ccc230321db5ab801ad1f0d7db3ef1e11f497b804eb546dae842b1a8c2e`.
- Source/runtime map and machine-readable tool, alignment, refiner, matched-
  candidate, and no-drop ledgers are now under `evidence/`; they preserve
  unknown causes and explicitly distinguish provider score semantics.
- The mapped daily-organizer checkout is present but outside the writable
  boundary and has unrelated dirty files; the synchronization blocker and
  required target filename are recorded in
  `evidence/daily-organizer-sync-blocker.md`.
- Evidence-change critical-region review: `evidence/change-review.md`.
- Completion-gate audit: `evidence/completion-audit.md`; requirements remain
  partial, unknown, or externally blocked, so the goal stays active.

## Required next action

After `scontrol ping` and accounting report a healthy controller, reconcile whether job `1657005` ever ran. If no durable artifacts exist, submit exactly one fresh copy of `slurm/step5_validation.sbatch`. Inspect all three staged arms and record alignment JSON count, transformed PDB count, selected count, PRODIGY states, return codes, and failure/no-drop disposition. Do not merge or promote.

## Do not infer

Do not infer pipeline success from the job ID, `sbatch` exit status, model confidence, or prior retained jobs. Do not infer native ranking quality from the 4-to-1 load reduction. Do not pool MultiProt's Kabsch-RMSD proxy with TMalign, GTalign, or USalign TM-scores, and do not treat the disconnected `src/structural_aligner.py` as the maintained runtime path.
