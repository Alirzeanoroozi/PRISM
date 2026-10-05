# Step 5 checkpoint — bounded validation

Recorded: 2026-09-09 (Europe/Istanbul)

## Gate state

`PARTIAL_BLOCKED_BY_SCHEDULER`.

The deterministic test and retained bounded-runtime gates pass. The fresh Step 5 PRISM smoke is `UNKNOWN`, not failed or successful, because repeated reconciliations through `2026-09-09T02:51:38+03:00` found both Slurm controllers down after job `1657005` submission and no durable job output appeared. The latest check from `ai01` also could not contact accounting; no replacement was submitted.

## Completed

- Canonical focused suite: 22 passed with `gtalign_env/bin` on PATH (2.03s).
- Isolated USalign suite: 20 passed (1.78s).
- Isolated PRODIGY suite: 9 passed (0.44s).
- Canonical defaults reviewed: TMalign, ranking disabled, baseline method, top-k 5, refinement enabled.
- Retained current-tree paired smoke reviewed: 2 audit records, 4 transformed pairs; baseline selected/refined 4/2, ranked top-1 selected/refined 1/1; both return code 0.
- Retained USalign Slurm probe reviewed: job 1656966, 214 matches, TM-score 0.98359.
- Retained PRODIGY probes reviewed: success job 1656993 (`available`, `executed`, affinity -65.827) and failure-preservation job 1656992 (`available`, `failed`, selected count 2).
- Canonical Git status remains 373 entries with unchanged hash `f63a1ccc230321db5ab801ad1f0d7db3ef1e11f497b804eb546dae842b1a8c2e`.
- Source review confirms `prism.py` directly calls the maintained TMalign,
  GTalign, and MultiProt providers; `src/structural_aligner.py` is disconnected
  and USalign remains isolated. Tool provenance is recorded in
  `evidence/tool-preflight-20260909T021142+0300.json`.
- Machine-readable tool, alignment, refiner, matched-candidate, and no-drop
  ledgers were added under `evidence/`. The retained no-drop accounting records
  73,251 missing/unwritten alignment records without inventing their cause and
  preserves the 1,877 residual as unclassified rather than structural orphanage.
- The tool matrix explicitly retains all planned Tier 1–4 arms: 31 rows with
  completed, not-comparable, or scheduler-blocked status; no planned arm was
  silently omitted.

## Not completed

- Fresh one-pair/one-template PRISM smoke: job 1657005 remains unreconciled and has produced no visible artifacts.
- Fresh current-tree USalign and PRODIGY integration runs.
- Native DockQ comparison, full template-panel rerun, and promotion review.
- GPU/CPU/multithread GTalign comparison and complete per-candidate Rosetta
  return-code/score-gate accounting.

## Safety decision

Do not merge, promote, reset, delete, or change production defaults. Preserve job `1657005` as a scheduler-state-unknown submission and rerun the same bounded script only after Slurm controller/accounting health is restored and the job identity is reconciled. The latest check is recorded in `evidence/scheduler-recheck-20260909T025138+0300.md`.
