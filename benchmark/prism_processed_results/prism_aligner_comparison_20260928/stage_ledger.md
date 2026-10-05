# PRISM aligner comparison stage ledger

Date: 2026-09-28  
Selected project: `/scratch/rshadi25/GitHub/PRISM-prescript`  
Primary historical evidence: `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919/compact_results_20260927`  
Panel lane: historical/working 19,855 templates, 257 BM5.5 cases  
Requested manifest lanes are frozen in `exact_panel_evidence/package_manifest.json`:
19,948 checked; 19,062 calculated; 19,058 materialized; 19,855 historical.

Status vocabulary: `VALID_REUSABLE`, `VALID_BUT_DIFFERENT_PANEL`, `PARTIAL`,
`INVALID_CONTRACT`, `MISSING`, `NOT_REQUIRED`.

| stage | TMalign | MultiProt | USalign | GTalign |
|---|---|---|---|---|
| alignment | `VALID_REUSABLE`: 257/257 cases; 20,410,940 observed records; compact validation | `VALID_REUSABLE`: 257/257 cases; 20,410,940 observed records; compact validation | `VALID_REUSABLE` for bounded pilot: job 1708928 completed default/`-fast` at 1/2/4/8/16 workers, each 1,892/1,892 records; exact current-panel jobs 1659360/1659361 remain alignment-only | `VALID_BUT_DIFFERENT_PANEL`: historical full GPU search, 257 cases and 724,544 parsed records; exact current-panel alignment-only jobs 1659356/1659358 are retained, no downstream quality claim |
| transformation | `VALID_REUSABLE`: 18,171 retained candidate rows/groups in compact candidate audit | `VALID_REUSABLE`: 70,104 retained candidate rows/groups in compact candidate audit | `SUBMITTED`: common replay/compaction array 1709007 queued after USalign provider 1708992; corrected transformed-DockQ array 1709046 queued after compaction | `VALID_REUSABLE` for historical run: 15,440 transformed groups / 30,880 receptor+ligand parts, subject to panel label |
| filtering | `VALID_REUSABLE`: 3,068 generated; 15,103 clash rejected; explicit threshold failures | `VALID_REUSABLE`: 66,827 generated; 3,277 clash rejected; explicit threshold failures; old RMSD proxy not TM-score | `SUBMITTED`: replay writes explicit generated/rejection summaries before deleting raw alignment JSON | `VALID_REUSABLE` in historical run; common-contract reconciliation pending |
| ranking | `MISSING` for this matched lane; historical run used rank false | `MISSING` for this matched lane; historical run used rank false | `NOT_REQUIRED` until alignment/transformation exists | `VALID_BUT_DIFFERENT_PANEL`: rank false historical source; top-k analysis still required |
| refinement | `PARTIAL`: selected candidates are in active common external-Rosetta run | `PARTIAL`: selected candidates are in active common external-Rosetta run | `NOT_REQUIRED` until upstream exists | `SUBMITTED`: 14,470 retained GTalign models; common worker array 1709167, execution/scientific validation pending |
| DockQ | `PARTIAL`: corrected evaluator exists; compact final GlobalDockQ/cross-interface reaggregation pending | `PARTIAL`: corrected evaluator exists; compact final GlobalDockQ/cross-interface reaggregation pending | `NOT_REQUIRED` until valid transformed candidates exist | `VALID_REUSABLE` for historical transformed lane: 15,440 rows, 14,300 scored and 1,140 score_failed across 216 cases; corrected GlobalDockQ/cross-interface fields retained, refined lane pending |
| timing | `VALID_REUSABLE` for historical alignment/stage timings; not a clean speed ranking; exact 946/19058 alignment-only timings also retained | `VALID_REUSABLE` for historical alignment/stage timings; candidate-count confounded | `VALID_REUSABLE` for bounded alignment pilot: default_w16 105.337 records/s, `-fast`_w16 110.677 records/s, zero failures; `-fast` changed 1,385/1,892 mappings/scores and 1,306 transforms, so production selection is default with 16 workers | `PARTIAL`: historical GPU timing available; exact current-panel alignment-only timing retained; matched reconciliation pending |

## Evidence counts

- Expected BM5.5 cases: 257.
- Historical compact source: 19,855 templates; 79,420 expected records/case;
  20,410,940 expected records/aligner.
- TMalign compact candidate audit: 257 audit files, 10,205,470 audit records,
  18,171 retained rows, 3,068 generated, 15,103 clash rejected.
- MultiProt compact candidate audit: 257 audit files, 10,205,470 audit records,
  70,104 retained rows, 66,827 generated, 3,277 clash rejected.
- Historical candidate identity overlap (template, query pair, orientation,
  chain pair) is compacted for all 257 cases under
  `candidate_overlap/historical_tmalign_multiprot_gtalign/`: TMalign versus
  MultiProt has 82 shared candidates and mean case-wise Jaccard `0.0014676`;
  GTalign versus TMalign has 210 shared and mean Jaccard `0.0532620`; GTalign
  versus MultiProt has 19 shared and mean Jaccard `0.0003044`. These are
  candidate-space observations only; TMalign has candidates in 195 cases,
  MultiProt in 257, and GTalign in 216.
- USalign historical attempt: 1 failed case, 256 not started, 79,420 observed
  records in the failed/partial lane, zero complete transformed candidates.
- Existing bounded 946-template arm: 10 rigid cases, 946 templates, 37,840
  alignment records; TMalign completed with 37,840 successes and 3 transformed
  pairs; MultiProt completed with 2,621 successes and 132 transformed pairs;
  the USalign-labelled arm completed 37,840 records but produced zero
  transformed pairs and its recorded command/progress is not sufficient to
  establish a clean USalign worker comparison.
- GTalign historical source: 724,544 parsed alignment records and 15,440
  transformed groups; corrected transformed scoring has 14,300 scored and
  1,140 score_failed rows across 216 cases. This is not silently pooled with
  the exact 19,058 lane or treated as a refined result. The compact score CSV,
  summary, input manifest, and hashes are retained under
  `gtalign_historical_19855_transformed/`.
- GTalign common-refinement preparation retained 14,470 materialized models
  and recorded 970 rejected rows (missing retained model/native consumer).
  Array 1709167 (`0-144%8`, 100 candidates/shard) runs the same
  FiberDock/external-Rosetta/corrected-DockQ worker contract as the active
  TMalign/MultiProt refinement. The prior artifact poll found 91 checkpoint
  records: 86 top-level completed, 4 running, and 1 failed; external Rosetta
  has 43 completed and 43 explicit `no_model` outcomes. This is execution
  evidence, not scientific validation.
  Read-only aggregation job 1709178 (test-only 1709177) is dependency-linked
  after the refinement array and will require all 14,470 checkpoints before
  producing cleanup eligibility.
- Exact current-panel alignment-only evidence: checked-prefix 946 is a
  one-query, 1,892-interface-side run; materialized calculated panel is a
  one-query, 38,116-interface-side run. TMalign, USalign, and GTalign have
  completed records for these lanes, but transformation/filtering/ranking/
  refinement/DockQ were not run on this exact panel.
- Corrected USalign pilot job 1708928 passed 10 configurations (default and
  `-fast`, workers 1/2/4/8/16), each with 1,892/1,892 execution successes and
  zero failures. The fastest measured configuration was `-fast`/16 workers,
  but it changed 1,385/1,892 mappings and normalized scores relative to
  default/16; therefore the scientifically selected production configuration
  is default USalign with 16 workers. Compact pilot evidence and cleanup
  manifest are retained under `usalign_pilot_946/`.
- Historical-panel USalign production was submitted as resumable VALAR array
  `1708992` (`1-26%4`) after dry-run `1708991`, using default USalign and 16
  workers with refinement disabled. At ledger update time, four tasks were
  running and the remainder pending; this is `SUBMITTED`, not execution or
  scientific validation. Batch completion markers and compact counts are
  required before any production result is accepted.
  Latest read-only poll found all four provider tasks still active, with
  approximately 660,648 / 754,490 / 595,650 / 753,537 raw JSON files in
  batches 1--4 and no compact completion markers. The raw record trees are
  therefore not eligible for deletion.
- A dependent postprocess array `1709007` (dry-run `1709006`) is queued
  `afterany:1708992`. It will replay only the existing common transformation
  with candidate-audit output, compact per-case dual TM-score/candidate
  summaries, and delete raw alignment JSON only after those summaries
  validate.
- A dependent corrected transformed-DockQ array `1709046` (dry-run `1709044`)
  is queued `afterany:1709007`. It consumes compact generated-candidate rows,
  uses the existing bijective GlobalDockQ/cross-interface evaluator, writes
  per-candidate checkpoints and TSVs, and deletes only per-candidate combined
  scoring scratch after those compact records are durable. Transformed USalign
  halves are intentionally retained for the pending common-refinement consumer;
  score failures retain them for audit as well. Cleanup is gated on common
  refinement and paired DockQ validation.
- A dependent read-only USalign batch aggregation job `1709192` (dry-run
  `1709191`) is queued `afterany:1709046`. It will merge only validated
  per-batch candidate/DockQ/interface tables into one compact aggregate and
  will not delete or copy raw alignment/transformation files. The tested
  `prepare_usalign_refinement_manifest.py` step is gated on that aggregate
  reporting complete and on all referenced transformed halves remaining
  materialized.
- A fail-closed USalign refinement-preparation job `1709233` (dry-run
  `1709232`) is queued `afterany:1709192`. It will create the common-refinement
  selected/rejected manifest only when the batch aggregation reports
  `validated_compacted`; it does not submit an array or delete structures.
- Latest GTalign common-refinement artifact poll at 2026-09-28T05:54:41+03:00
  contains 92 checkpoint records: 87 top-level completed, 4 running, and 1
  failed. The failed row is
  `medium_1wq1_045 / 1de4AC / o1`, rejected during input normalization because
  native ligand `G*` requires two chain segments while the retained assembled
  model exposes one predicted chain `D`; it is retained as an explicit
  failure, not reconstructed with an unsupported chain assumption. FiberDock
  has 90 completed stages and 86 DockQ scores; external Rosetta has 43
  completed and 43 explicit `no_model` outcomes.
  These are execution facts, not quality results.
- The previously observed KUACC TMalign/MultiProt common-refinement task
  `3132295_874` resolved to child `3140177` with exit `0:0`; controller
  `3132296` also completed with exit `0:0` and scheduled continuation arrays
  `3140365` and `3140401`--`3140407`. The remote run currently exposes 16,651
  checkpoint files but no aggregate package, so it remains execution-pending;
  no KUACC cleanup or job intervention was performed.

## Contract and cleanup gates

- TM scores: retain query/Structure_1 and reference/Structure_2 independently;
  common score is explicitly reference-normalized, never `max()`.
- Compact USalign alignment summaries now retain explicit unavailable/truncated
  status counts but calculate numerical TM-score and match-count means only
  from complete successful records; the contract delta is recorded in
  `usalign_production_19855/source_snapshot_after_alignment_summary_contract_fix.json`
  and was tested with the provider still active. No provider rerun is needed.
- Final matched aggregation was strengthened before execution: ranking output
  now distinguishes no-ranking/all, deterministic baseline, PRODIGY when
  retained affinities exist, top-1/top-3/top-5, and diagnostic oracle; it also
  emits transformed/refined coverage, cross-interface quality summaries, and a
  method-by-split summary table.
  Focused tests pass 3/3, with provenance in
  `usalign_production_19855/source_snapshot_after_final_ranking_contract_fix.json`.
- DockQ: retain GlobalDockQ, requested cross-interface DockQ, interface scope,
  Fnat/iRMSD/LRMSD, mapping, status, evaluator/version, hashes, elapsed time.
- USalign cleanup: `cleanup_usalign_run.py` is dry-run/apply tested and
  requires validated batch aggregation, refinement handoff, common-refinement
  aggregate eligibility, retained-artifact hashes, and final-package
  validation. It has not been applied.
- No intermediate deletion is approved by this ledger until all downstream
  consumers finish and the corresponding compact package, provenance, hashes,
  counts, failures, and commands validate.  The active large-refinement root
  and its user-owned jobs remain untouched.

## Next acceptance gate

Latest scheduler poll at 2026-09-28T06:12:29+03:00 found USalign task
1708992_1 validated at exit 0:0 with `available=10 unavailable=0` and empty
stderr; tasks 3, 5, 6, and 7 remain active, later tasks are pending, and no
compact marker exists. All post-provider stages remain dependency-held. GTalign
now has 96 checkpoint files (91 completed, 4 running, 1 explicit failure). This is live
execution evidence only; no cleanup is eligible.

The requested template-panel counts and bounded USalign pilot are validated;
the resumable historical-panel production batches are submitted with default
USalign and 16 workers. Validate provider execution, common
transformation/filtering, and corrected transformed DockQ before any
refinement interpretation. The identical
pair correctness gate is complete: TMalign and USalign both parsed 28 aligned
residues, identical mappings/transforms, RMSD 0.00, and both normalized
TM-scores 1.0; common transformation outputs were byte-identical.  Production
USalign and final EDA remain blocked until the worker sweep and full-panel
production evidence pass.

## Live checkpoint — 2026-09-28T06:22:15+03:00

The corrected GTalign checkpoint probe found 97 files: 92 completed, 4
running, and 1 explicit failure. The failure remains the recorded
`medium_1wq1_045 / 1de4AC / o1` input-normalization/cardinality mismatch.
USalign provider tasks `1708992_3`, `1708992_5`, `1708992_6`, and `1708992_7`
remain running; all dependent compact, DockQ, aggregation, and handoff jobs
remain pending. No cleanup is eligible.

## Live checkpoint — 2026-09-28T06:25:50+03:00

The authorized KUACC common-refinement continuation is progressing: array
`3140365` has a running task `_919`, and array `3140401` has running tasks
through `_786/_787`; later shards remain pending under `AssocMaxJobsLimit`.
The large-tree remote checkpoint count was not updated because its read-only
probe timed out. No cleanup or job intervention was performed.

## Live checkpoint — 2026-09-28T06:27:19+03:00

GTalign array `1709167` now exposes 99 checkpoint files: 94 completed, 4
running, and 1 explicit failure. KUACC continuation `3140401` reached running
task `_862`, while `3140365` remains active at `_919`. USalign provider and
all dependent local stages remain active or dependency-held; no cleanup is
eligible.

## Live checkpoint — 2026-09-28T06:28:26+03:00

GTalign array `1709167` now exposes 100 checkpoint files: 95 completed, 4
running, and 1 explicit failure. KUACC continuation `3140401` advanced to
task `_900`, with `3140365` still active at `_919`. USalign and all dependent
local stages remain active or dependency-held; no cleanup is eligible.

## Live checkpoint — 2026-09-28T06:29:06+03:00

KUACC continuation `3140401` advanced to running task `_919`, while
`3140365` remains active at `_919`. GTalign remains at 100 checkpoint files
(95 completed, 4 running, 1 explicit failure); USalign and all dependent local
stages remain active or held. No cleanup is eligible.

## Live checkpoint — 2026-09-28T06:30:32+03:00

KUACC continuation `3140401` advanced to running task `_964`, while `3140365`
remains active at `_919`. GTalign remains at 100 checkpoint files (95
completed, 4 running, 1 explicit failure). No local compact, corrected-DockQ,
aggregation, or cleanup marker exists.

## Live checkpoint — 2026-09-28T06:32:05+03:00

KUACC continuation `3140401` reached tasks `_997`–`_999`, and shard `3140402`
opened tasks `0`–`10`; later tasks remain association-limit pending. Local
USalign and GTalign workers remain running, with downstream stages pending on
their declared dependencies.

## Live checkpoint — 2026-09-28T06:31:24+03:00

KUACC continuation `3140401` advanced to running task `_988`, while `3140365`
remains active at `_919`. Local USalign and GTalign workers remain running and
all declared downstream jobs remain pending. A broad marker scan timed out and
was not treated as evidence of completion.

## Live checkpoint — 2026-09-28T06:34:29+03:00

KUACC refinement shard `3140402` advanced to running task `_93`; tasks
`_94–999` remain association-limit pending. Local USalign/GTalign workers
remain active with downstream stages held; GTalign remains at 95 completed, 4
running, 1 explicit failure.

## Live checkpoint — 2026-09-28T06:32:41+03:00

KUACC refinement shard `3140402` advanced to running tasks through `_31`; its
remaining tasks are association-limit pending. Local USalign/GTalign workers
remain running and downstream stages remain dependency-held. GTalign remains
at 100 checkpoints (95 completed, 4 running, 1 explicit failure).

## Live checkpoint — 2026-09-28T06:35:30+03:00

KUACC refinement shard `3140402` advanced to running task `_122`; tasks
`_123–999` remain association-limit pending. Local provider/refinement arrays
remain active and GTalign remains at 95 completed, 4 running, 1 explicit
failure. No dependency or cleanup gate has opened.

## Live checkpoint — 2026-09-28T06:33:43+03:00

KUACC refinement shard `3140402` advanced to running task `_65`; tasks
`_66–999` remain association-limit pending. Local USalign/GTalign workers
remain active and downstream jobs remain held; GTalign remains at 95 completed,
4 running, 1 explicit failure.

## Live checkpoint — 2026-09-28T06:46:49+03:00

KUACC refinement shard `3140402` advanced to running task `_511`; tasks
`_512–999` remain association-limit pending. GTalign remains at 96 completed,
4 running, 1 explicit failure; USalign task `1708992_8` remains active and
downstream local stages remain held.

## Live checkpoint — 2026-09-28T06:36:08+03:00

KUACC refinement shard `3140402` advanced to running task `_146`; tasks
`_147–999` remain association-limit pending. Local USalign/GTalign workers
remain active; GTalign remains at 95 completed, 4 running, 1 explicit failure.

## Live checkpoint — 2026-09-28T06:45:40+03:00

KUACC refinement shard `3140402` advanced to running task `_471`; tasks
`_472–999` remain association-limit pending. GTalign remains at 101
checkpoints (96 completed, 4 running, 1 explicit failure). USalign task
`1708992_8` remains active and downstream local stages remain held.

## Live checkpoint — 2026-09-28T06:37:59+03:00

KUACC refinement shard `3140402` advanced to running task `_208`; tasks
`_209–999` remain association-limit pending. Local USalign/GTalign workers
remain active; GTalign remains at 95 completed, 4 running, 1 explicit failure.
No downstream gate has opened.

## Live checkpoint — 2026-09-28T06:38:43+03:00

USalign provider task `1708992_8` is now running while tasks `9–26` remain
array-limit pending. KUACC refinement shard `3140402` advanced to task `_234`.
GTalign remains at 95 completed, 4 running, 1 explicit failure; local
downstream stages remain dependency-held.

## Live checkpoint — 2026-09-28T07:06:46+03:00

Authoritative scheduler state remains nonterminal: GTalign workers `0–3`,
USalign tasks `1708992_5–8`, and common-refinement tasks `1709167_0–3` are
running. Their dependent compact/DockQ/aggregate stages remain held. KUACC
shard `3140403` has progressed through task `_172`; tasks `_173–999` remain
association-limited while late shard-18 tasks continue draining. No newly
validated compact consumer or cleanup gate has opened.

## Live checkpoint — 2026-09-28T07:04:19+03:00

Read-only GTalign checkpoint audit reports 103 completed, 4 running, and 1
explicit input-normalization failure. All completed records contain the
expected input-normalization, FiberDock, both DockQ, and external-Rosetta
stage keys; no completed record has a missing expected stage key. KUACC shard
`3140403` and VALAR USalign/common-refinement workers remain active, so
downstream compaction, DockQ, and aggregation gates remain held.

## Live checkpoint — 2026-09-28T07:05:09+03:00

GTalign completed-record stage audit: 103/103 passed input normalization and
FiberDock; 51/103 produced an external-Rosetta model and scored with DockQ;
52/103 are explicitly `no_model`/`not_run_no_model` because the refiner
produced no model. These remain explicit valid-unscored states, not zeros.
Active KUACC/VALAR jobs and dependency-held downstream stages are unchanged;
no cleanup is authorized until compact consumers validate.

## Live checkpoint — 2026-09-28T07:05:56+03:00

GTalign now has 104 completed, 4 running, and 1 explicit input-normalization
failure. Refreshed audit: no completed record is missing an expected stage key;
external Rosetta/DockQ are `52 scored` and `52 no_model`/
`not_run_no_model`. KUACC shard `3140403` is active while late shard-18 tasks
drain; VALAR USalign/common-refinement workers remain active and all downstream
gates stay dependency-held.

## Live checkpoint — 2026-09-28T07:03:06+03:00

KUACC refinement shard `3140403` advanced through running task `_54`; tasks
`_55–999` remain association-limit pending. Shards `3140404–3140407` and
controller `3140408` remain held. GTalign advanced to 102 completed, 4 running,
1 explicit input-normalization failure. USalign tasks `1708992_5–8` remain
active, tasks `9–26` remain array-limit pending, and downstream stages remain
dependency-held.

## Live checkpoint — 2026-09-28T06:39:34+03:00

USalign task `1708992_8` remains running while tasks `9–26` remain array-limit
pending. KUACC shard `3140402` advanced to `_262`; later tasks remain
association-limited. GTalign remains at 95 completed, 4 running, 1 explicit
failure.

## Live checkpoint — 2026-09-28T06:44:45+03:00

GTalign common refinement now has 101 checkpoint files: 96 completed, 4
running, and 1 explicit failure. KUACC refinement shard `3140402` advanced to
task `_440`; tasks `_441–999` remain association-limit pending. Local
downstream stages remain dependency-held.

## Live checkpoint — 2026-09-28T06:44:02+03:00

KUACC refinement shard `3140402` advanced to running task `_415`; tasks
`_416–999` remain association-limit pending. USalign task `1708992_8` and
GTalign workers remain active; GTalign remains at 95 completed, 4 running, 1
explicit failure.

## Live checkpoint — 2026-09-28T06:43:21+03:00

KUACC shard `3140402` advanced to running task `_390`; tasks `_391–999`
remain association-limit pending. USalign task `1708992_8` and GTalign
workers remain active; GTalign remains at 95 completed, 4 running, 1 explicit
failure.

## Live checkpoint — 2026-09-28T06:40:31+03:00

USalign task `1708992_8` remains active; tasks `9–26` remain array-limit
pending. KUACC shard `3140402` advanced to `_288`; later tasks remain
association-limited. GTalign remains at 95 completed, 4 running, 1 explicit
failure.

## Live checkpoint — 2026-09-28T06:42:01+03:00

KUACC shard `3140402` advanced to running task `_341`; tasks `_342–999`
remain association-limit pending. USalign task `1708992_8` and all four
GTalign workers remain active; GTalign remains at 95 completed, 4 running, 1
explicit failure.

## Live checkpoint — 2026-09-28T06:42:43+03:00

KUACC shard `3140402` advanced to running task `_369`; tasks `_370–999`
remain association-limit pending. USalign task `1708992_8` and all GTalign
workers remain active; GTalign remains at 95 completed, 4 running, 1 explicit
failure.

## Live checkpoint — 2026-09-28T06:58:49+03:00

KUACC refinement shard `3140402` advanced through running task `_924`; tasks
`_925–999` remain association-limit pending. Shards `3140403–3140407` and
controller `3140408` remain dependency/association-limited. GTalign advanced to
98 completed, 4 running, 1 explicit input-normalization failure. USalign task
`1708992_8` remains active, tasks `9–26` remain array-limit pending, and all
downstream USalign stages remain dependency-held.

## Live checkpoint — 2026-09-28T07:01:38+03:00

KUACC refinement shard `3140402` is draining its final tasks while shard
`3140403` has started (tasks `_0–4` running; `_5–999` association-limited).
Shards `3140404–3140407` and controller `3140408` remain held. GTalign advanced
to 101 completed, 4 running, 1 explicit input-normalization failure. USalign
tasks `1708992_5–8` remain active, tasks `9–26` remain array-limit pending, and
downstream stages remain dependency-held.

## Live checkpoint — 2026-09-28T07:07:25+03:00

KUACC shard `3140403` advanced to task `_197`; tasks `_198–999` remain
association-limited, with late shard-18 tasks still running. GTalign remains
104 completed, 4 running, 1 explicit failure. VALAR USalign tasks
`1708992_5–8` and common-refinement tasks `1709167_0–3` remain active; all
dependent stages remain pending by dependency.

## Live checkpoint — 2026-09-28T07:08:30+03:00

GTalign advanced to 105 completed, 4 running, 1 explicit failure. KUACC
refinement shard `3140403` advanced to task `_225`; tasks `_226–999` remain
association-limited, with late shard-18 tasks still running. VALAR USalign
tasks `1708992_5–8` and common-refinement tasks `1709167_0–3` remain active;
no downstream compact or DockQ output has appeared.

## Live checkpoint — 2026-09-28T07:09:02+03:00

KUACC shard `3140403` advanced to task `_240`; task `_211` is completing and
tasks `_241–999` remain association-limited. Late shard-18 tasks remain active.
GTalign remains 105 completed, 4 running, 1 explicit failure. VALAR USalign
tasks `1708992_5–8` and common-refinement tasks `1709167_0–3` remain active;
no downstream marker or validated compact output has appeared.

## Live checkpoint — 2026-09-28T07:09:34+03:00

KUACC refinement shard `3140403` advanced to task `_258`; tasks `_259–999`
remain association-limited, with late shard-18 task `_578` still active.
GTalign remains 105 completed, 4 running, 1 explicit failure. VALAR USalign
and common-refinement workers remain active; the new checkpoint file is not a
terminal downstream gate.

## Live checkpoint — 2026-09-28T07:10:06+03:00

KUACC refinement shard `3140403` advanced to task `_274`; tasks `_275–999`
remain association-limited, while shard-18 task `_578` remains active. GTalign
remains 105 completed, 4 running, 1 explicit failure. VALAR USalign and
common-refinement workers remain active; no dependency has opened.

## Live checkpoint — 2026-09-28T07:10:38+03:00

KUACC refinement shard `3140403` advanced to task `_292`; tasks `_293–999`
remain association-limited, while shard-18 task `_578` remains active. GTalign
remains 105 completed, 4 running, 1 explicit failure. VALAR USalign and
common-refinement workers remain active; no newly validated downstream artifact
or dependency transition is present.

## Live checkpoint — 2026-09-28T07:11:36+03:00

KUACC refinement shard `3140403` advanced to task `_319`; task `_282` is
completing and tasks `_320–999` remain association-limited. Shard-18 task
`_578` remains active. GTalign remains 105 completed, 4 running, 1 explicit
failure. VALAR USalign and common-refinement workers remain active; no
downstream dependency has opened.

## Live checkpoint — 2026-09-28T07:12:07+03:00

KUACC refinement shard `3140403` advanced to task `_334`; task `_312` is
completing and tasks `_335–999` remain association-limited. Shard-18 task
`_578` remains active. GTalign remains 105 completed, 4 running, 1 explicit
failure. VALAR USalign and common-refinement workers remain active at
approximately four hours; no downstream marker or validated compact artifact
exists.

## Live checkpoint — 2026-09-28T07:12:39+03:00

KUACC refinement shard `3140403` advanced to task `_353`; tasks `_354–999`
remain association-limited, while shard-18 task `_578` remains active. GTalign
remains 105 completed, 4 running, 1 explicit failure. VALAR USalign and
common-refinement workers remain active; no dependency or compact-result
transition has occurred.

## Live checkpoint — 2026-09-28T07:13:37+03:00

GTalign advanced to 106 completed, 4 running, 1 explicit failure. KUACC
refinement shard `3140403` advanced to task `_385`; tasks `_387–999` remain
association-limited, with long-running shard-18 task `_578` active. The new
GTalign checkpoint is within the active common-refinement run; no downstream
aggregation or cleanup gate has opened.

## Live checkpoint — 2026-09-28T07:14:10+03:00

KUACC refinement shard `3140403` advanced to task `_402`; tasks `_403–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. No downstream compact, DockQ, aggregation, or cleanup artifact has
appeared.

## Live checkpoint — 2026-09-28T07:16:26+03:00

KUACC refinement shard `3140403` advanced to task `_476`; tasks `_477–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR USalign/common-refinement workers remain active; no downstream
validation artifact has appeared.

## Live checkpoint — 2026-09-28T07:14:57+03:00

KUACC refinement shard `3140403` advanced to task `_430`; tasks `_431–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR USalign and common-refinement jobs remain active or
dependency-held; no downstream artifact has opened.

## Live checkpoint — 2026-09-28T07:15:53+03:00

KUACC refinement shard `3140403` advanced to task `_460`; tasks `_461–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no validated compact consumer has
appeared.

## Live checkpoint — 2026-09-28T07:17:18+03:00

KUACC refinement shard `3140403` advanced to task `_502`; tasks `_503–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no dependency or validated downstream
artifact has opened.

## Live checkpoint — 2026-09-28T07:21:17+03:00

KUACC refinement shard `3140403` advanced to task `_624`; tasks `_625–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no downstream dependency or compact
result has opened.

## Live checkpoint — 2026-09-28T07:22:24+03:00

KUACC refinement shard `3140403` advanced to task `_656`; tasks `_657–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no terminal dependency or compact output
has appeared.

## Live checkpoint — 2026-09-28T07:30:00+03:00

KUACC refinement shard `3140403` advanced to task `_886`; tasks `_887–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR jobs remain active or dependency-held; no downstream gate has
opened.

## Live checkpoint — 2026-09-28T07:48:56+03:00

KUACC shard-20 advanced through task `_460`; tasks `_461–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:58:34+03:00

GTalign remains `116 completed`, `4 running`, and `2 explicit failures`; the
failure audit is unchanged and zero completed records are missing required
stages. KUACC shard-20 task `_739` is active; tasks `_740–999` remain
association-limited. VALAR workers remain active and dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:51:41+03:00

GTalign checkpoint audit now reports `113 completed`, `4 running`, and `1
explicit failure`; zero completed records are missing required stages. KUACC
shard-20 advanced through task `_538`; tasks `_539–999` remain
association-limited. VALAR USalign and common-refinement workers remain active,
with dependent jobs held.

## Live checkpoint — 2026-09-28T07:52:18+03:00

KUACC shard-20 task `_523` reached terminal `COMPLETED` status; tasks through
`_558` remain active and `_559–999` are association-limited. GTalign remains
`113 completed`, `4 running`, and `1 explicit failure`; zero completed records
are missing required stages. VALAR workers remain active, with dependent jobs
still held.

## Live checkpoint — 2026-09-28T07:51:08+03:00

KUACC shard-20 advanced through task `_522`; tasks `_523–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:50:34+03:00

KUACC shard-20 advanced through task `_508`; tasks `_509–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:50:04+03:00

KUACC shard-20 advanced through task `_494`; tasks `_495–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:49:30+03:00

KUACC shard-20 advanced through task `_477`; tasks `_478–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:48:25+03:00

KUACC shard-20 advanced through task `_447`; tasks `_448–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; no dependent stage has
opened.

## Live checkpoint — 2026-09-28T07:47:54+03:00

KUACC shard-20 advanced through task `_431`; tasks `_433–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; no dependent stage has
opened.

## Live checkpoint — 2026-09-28T07:47:19+03:00

KUACC shard-20 advanced through task `_412`; tasks `_413–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; no USalign or refinement
dependency is terminally validated yet.

## Live checkpoint — 2026-09-28T07:46:44+03:00

KUACC shard-20 advanced through task `_395`; tasks `_396–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; no dependent stage has
opened.

## Live checkpoint — 2026-09-28T07:31:23+03:00

KUACC refinement shard `3140403` advanced to task `_925`; tasks `_926–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no terminal dependency has opened.

## Live checkpoint — 2026-09-28T07:32:11+03:00

GTalign advanced to 107 completed, 4 running, 1 explicit failure. Contract
audit: no completed record is missing an expected stage key; external
Rosetta/DockQ are `55 scored` and `52 no_model`/`not_run_no_model`. KUACC
refinement shard `3140403` advanced to task `_949`; tasks `_950–999` remain
association-limited, while shard-18 task `_578` and shard-16 task `_919` remain
active.

## Live checkpoint — 2026-09-28T07:30:47+03:00

KUACC refinement shard `3140403` advanced to task `_907`; task `_860` is
completing and tasks `_908–999` remain association-limited. Shard-18 task
`_578` and shard-16 task `_919` remain active. GTalign remains 106 completed,
4 running, 1 explicit failure. VALAR jobs remain active; no downstream gate
has opened.

## Live checkpoint — 2026-09-28T07:23:00+03:00

KUACC refinement shard `3140403` advanced to task `_672`; tasks `_673–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no dependency or validated compact output
has opened.

## Live checkpoint — 2026-09-28T07:23:33+03:00

KUACC refinement shard `3140403` advanced to task `_688`; tasks `_690–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR jobs remain active or dependency-held; no terminal dependency
or compact result has appeared.

## Live checkpoint — 2026-09-28T07:24:05+03:00

KUACC refinement shard `3140403` advanced to task `_707`; tasks `_708–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no terminal dependency or downstream
compact artifact has appeared.

## Live checkpoint — 2026-09-28T07:29:09+03:00

KUACC refinement shard `3140403` advanced to task `_862`; task `_844` is
completing and tasks `_863–999` remain association-limited. Shard-18 task
`_578` and shard-16 task `_919` remain active. GTalign remains 106 completed,
4 running, 1 explicit failure. VALAR workers remain active; no downstream gate
has opened.

## Live checkpoint — 2026-09-28T07:24:41+03:00

KUACC refinement shard `3140403` advanced to task `_728`; tasks `_729–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR stages remain active or dependency-held; no downstream artifact
has appeared.

## Live checkpoint — 2026-09-28T07:25:16+03:00

KUACC refinement shard `3140403` advanced to task `_745`; tasks `_746–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no terminal dependency or downstream
artifact has appeared.

## Live checkpoint — 2026-09-28T07:26:26+03:00

KUACC refinement shard `3140403` advanced to task `_784`; tasks `_785–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active or dependency-held; no downstream
artifact has opened.

## Live checkpoint — 2026-09-28T07:26:59+03:00

KUACC refinement shard `3140403` advanced to task `_802`; tasks `_803–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no dependency or downstream artifact has
opened.

## Live checkpoint — 2026-09-28T07:25:53+03:00

KUACC refinement shard `3140403` advanced to task `_766`; tasks `_767–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no terminal dependency or downstream
artifact has appeared.

## Live checkpoint — 2026-09-28T07:21:49+03:00

KUACC refinement shard `3140403` advanced to task `_636`; tasks `_637–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR jobs remain active or dependency-held; no validated downstream
output has appeared.

## Live checkpoint — 2026-09-28T07:20:09+03:00

KUACC refinement shard `3140403` advanced to task `_586`; tasks `_587–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no terminal dependency or new compact
output has appeared.

## Live checkpoint — 2026-09-28T07:20:45+03:00

KUACC refinement shard `3140403` advanced to task `_602`; tasks `_603–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no downstream dependency or compact
output has appeared.

## Live checkpoint — 2026-09-28T07:18:00+03:00

KUACC refinement shard `3140403` advanced to task `_526`; tasks `_527–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no dependency has opened.

## Live checkpoint — 2026-09-28T07:19:32+03:00

KUACC refinement shard `3140403` advanced to task `_571`; task `_552` is
completing and tasks `_572–999` remain association-limited. Shard-18 task
`_578` and shard-16 task `_919` remain active. GTalign remains 106 completed,
4 running, 1 explicit failure. VALAR workers remain active and downstream
dependencies remain held.

## Live checkpoint — 2026-09-28T07:18:59+03:00

KUACC refinement shard `3140403` advanced to task `_552`; tasks `_555–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no terminal dependency or compact result
has appeared.

## Live checkpoint — 2026-09-28T07:27:32+03:00

KUACC refinement shard `3140403` advanced to task `_816`; task `_792` is
completing and tasks `_818–999` remain association-limited. Shard-18 task
`_578` and shard-16 task `_919` remain active. GTalign remains 106 completed,
4 running, 1 explicit failure. VALAR jobs remain active or dependency-held; no
downstream artifact has opened.

## Live checkpoint — 2026-09-28T07:28:17+03:00

KUACC refinement shard `3140403` advanced to task `_839`; tasks `_840–999`
remain association-limited, while shard-18 task `_578` and shard-16 task
`_919` remain active. GTalign remains 106 completed, 4 running, 1 explicit
failure. VALAR workers remain active; no terminal dependency or downstream
compact artifact has appeared.

## Live checkpoint — 2026-09-28T07:32:59+03:00

KUACC refinement shard `3140403` advanced through task `_973`; tasks `_974–999`
remain association-limited. Shard-18 task `_578` and shard-16 task `_919` remain
active. GTalign remains 107 completed, 4 running, 1 explicit failure. VALAR
USalign and common-refinement workers remain active; downstream compact, DockQ,
and aggregate stages remain dependency-held.

## Live checkpoint — 2026-09-28T07:35:09+03:00

KUACC refinement shard `3140403` advanced through task `_987`; shard-18 task
`_578` and shard-16 task `_919` remain active, and shard-20 has started running
tasks. Later shard tasks remain association-limited. GTalign remains 107
completed, 4 running, 1 explicit failure. VALAR USalign and common-refinement
workers remain active; no downstream compact, DockQ, or aggregate artifact is
yet validated.

## Live checkpoint — 2026-09-28T07:35:55+03:00

GTalign checkpoint audit now reports `111 completed`, `4 running`, and `1
explicit failure` across 116 JSON records. KUACC shard-20 refinement tasks are
running; shard-19, shard-18, and shard-16 still have active work. VALAR USalign
and common-refinement workers remain active, with downstream stages
dependency-held.

## Live checkpoint — 2026-09-28T07:37:02+03:00

KUACC shard-20 advanced through task `_95`; tasks `_96–999` remain
association-limited. Shards 16, 18, and 19 still have active refinement work.
GTalign remains `111 completed`, `4 running`, and `1 explicit failure`. VALAR
USalign and common-refinement workers remain live; no downstream compact,
DockQ, or aggregate artifact is validated yet.

## Live checkpoint — 2026-09-28T07:37:55+03:00

KUACC shard-20 advanced through task `_123`; tasks `_124–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains
`111 completed`, `4 running`, and `1 explicit failure`. VALAR USalign and
common-refinement workers remain active, and downstream stages remain
dependency-held.

## Live checkpoint — 2026-09-28T07:38:49+03:00

GTalign contract audit: `112 completed`, `4 running`, and `1 explicit failure`;
all 112 completed records contain `fiberdock`, `external_rosetta`,
`dockq_fiberdock`, and `dockq_rosetta`. External Rosetta is 57 completed and
scored, with 55 explicit `no_model`/`not_run_no_model` states. KUACC shard-20
advanced through task `_153`; tasks `_154–999` remain association-limited. VALAR
USalign and common-refinement workers remain active; downstream stages remain
dependency-held.

## Live checkpoint — 2026-09-28T07:40:27+03:00

KUACC shard-20 advanced through task `_203`; tasks `_204–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`; the completed-record audit
reports zero missing required stages. VALAR USalign and common-refinement
workers remain active, with downstream stages held.

## Live checkpoint — 2026-09-28T07:39:53+03:00

KUACC shard-20 advanced through task `_187`; tasks `_188–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with the completed-record
stage contract valid. VALAR USalign and common-refinement workers remain active;
downstream stages remain dependency-held.

## Live checkpoint — 2026-09-28T07:41:01+03:00

KUACC shard-20 advanced through task `_219`; tasks `_220–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`; the completed-record audit
still reports zero missing required stages. VALAR workers remain active and no
downstream dependency has opened.

## Live checkpoint — 2026-09-28T07:41:50+03:00

KUACC shard-20 advanced through task `_247`; tasks `_248–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`; the completed-record audit
remains zero-missing-stage. VALAR USalign and common-refinement workers remain
active, with dependent jobs pending.

## Live checkpoint — 2026-09-28T07:42:23+03:00

KUACC shard-20 advanced through task `_264`; tasks `_265–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`; its completed-record audit
remains zero-missing-stage. VALAR workers remain active and no downstream output
is yet eligible for validation.

## Live checkpoint — 2026-09-28T07:42:55+03:00

KUACC shard-20 advanced through task `_279`; tasks `_280–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; no terminal downstream
batch is available yet.

## Live checkpoint — 2026-09-28T07:43:29+03:00

KUACC shard-20 advanced through task `_296`; tasks `_297–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; no compact or aggregate
output is yet validated.

## Live checkpoint — 2026-09-28T07:44:19+03:00

KUACC shard-20 advanced through task `_318`; tasks `_320–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; no dependent stage has
opened.

## Live checkpoint — 2026-09-28T07:44:51+03:00

KUACC shard-20 advanced through task `_340`; tasks `_341–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; all dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:46:11+03:00

KUACC shard-20 advanced through task `_377`; tasks `_378–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; no dependent stage has
opened.

## Live checkpoint — 2026-09-28T07:45:39+03:00

KUACC shard-20 advanced through task `_360`; tasks `_361–999` remain
association-limited. Shards 16, 18, and 19 remain active. GTalign remains `112
completed`, `4 running`, and `1 explicit failure`, with zero completed records
missing required stages. VALAR workers remain active; all dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:52:54+03:00

KUACC shard-20 task `_581` is active; tasks `_582–999` remain
association-limited. Earlier shard-20 task `_523` is terminally completed.
GTalign remains `113 completed`, `4 running`, and `1 explicit failure`, with
zero completed records missing required stages. VALAR USalign and
common-refinement workers remain active; dependent jobs remain held.

## Live checkpoint — 2026-09-28T07:54:10+03:00

GTalign audit now reports `114 completed`, `4 running`, and `2 explicit
failures`; all 114 completed records contain the required refinement/DockQ
stages. The second failure is the same documented `medium_1wq1_045`
input-normalization contract failure (`gtalign_source_*_R.pdb` has one chain
`D`, expected two). It remains explicit and unscored; no synthetic chain was
introduced. KUACC shard-20 has active tasks `_621–622`; `_623–999` remain
association-limited. VALAR workers remain active and dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:55:33+03:00

GTalign remains `114 completed`, `4 running`, and `2 explicit failures`; the
two failures are both the documented `medium_1wq1_045` input-normalization
case, and zero completed records are missing required stages. KUACC shard-20
task `_657` is active; tasks `_658–999` remain association-limited. VALAR
USalign and common-refinement workers remain active, with dependent jobs held.

## Live checkpoint — 2026-09-28T07:56:19+03:00

GTalign audit now reports `115 completed`, `4 running`, and `2 explicit
failures`; the two failures remain the documented `medium_1wq1_045`
input-normalization records, and zero completed records are missing required
stages. KUACC shard-20 task `_677` is active; tasks `_678–999` remain
association-limited. VALAR workers remain active and dependent jobs remain
held.

## Live checkpoint — 2026-09-28T07:57:16+03:00

GTalign remains `115 completed`, `4 running`, and `2 explicit failures`; the
same two `medium_1wq1_045` input-normalization failures are retained, with
zero completed records missing required stages. KUACC shard-20 task `_707` is
active; tasks `_708–999` remain association-limited. VALAR workers remain
active and dependent jobs remain held.

## Live checkpoint — 2026-09-28T07:57:54+03:00

GTalign audit now reports `116 completed`, `4 running`, and `2 explicit
failures`; only the two documented `medium_1wq1_045` normalization failures
remain, and zero completed records are missing required stages. KUACC shard-20
task `_725` is active; tasks `_727–999` remain association-limited. VALAR
workers remain active and dependent jobs remain held.

## Live checkpoint — 2026-09-28T07:59:12+03:00

GTalign audit now reports `117 completed`, `4 running`, and `2 explicit
failures`; the two failures remain the documented `medium_1wq1_045`
normalization records, and zero completed records are missing required stages.
KUACC shard-20 task `_761` is active; tasks `_762–999` remain
association-limited. VALAR workers remain active and dependent jobs remain
held.

## Live checkpoint — 2026-09-28T08:01:15+03:00

GTalign audit now reports `118 completed`, `4 running`, and `2 explicit
failures`; the two failures remain the documented `medium_1wq1_045`
normalization records, and zero completed records are missing required stages.
KUACC shard-20 task `_799` is active; tasks `_800–999` remain
association-limited. VALAR workers remain active and dependent jobs remain
held.

## Live checkpoint — 2026-09-28T08:02:54+03:00

GTalign audit now reports `119 completed`, `4 running`, and `2 explicit
failures`; the two failures remain the documented `medium_1wq1_045`
normalization records, and zero completed records are missing required stages.
KUACC shard-20 task `_854` is active; tasks `_855–999` remain
association-limited. VALAR workers remain active and dependent jobs remain
held.

## Live checkpoint — 2026-09-28T08:03:55+03:00

GTalign remains at `119 completed`, `4 running`, and `2 explicit failures`;
all completed records contain the required refinement/DockQ stages, and the
two failures remain the documented `medium_1wq1_045` normalization records.
KUACC shard-20 task `_898` is active; tasks `_899–999` remain
association-limited. VALAR workers remain active and dependent jobs remain
held.

## Live checkpoint — 2026-09-28T08:05:08+03:00

GTalign remains `119 completed`, `4 running`, and `2 explicit failures`; the
latest completed checkpoint is fresh, all completed records contain the
required refinement/DockQ stages, and the two failures remain the documented
`medium_1wq1_045` normalization records. USalign production batches 1–4 have
validated completion markers; batches 5–8 are active and their pipeline logs
are still advancing. KUACC shard-20 task `_898` is active; tasks `_899–999`
remain association-limited. Downstream aggregation/refinement jobs remain
correctly dependency-held.

## Live checkpoint — 2026-09-28T08:08:17+03:00

GTalign advanced to `120 completed`, with `4 running` and `2 explicit
failures`; completed-record stage validation remains clean, and the two
failures remain the documented `medium_1wq1_045` normalization records.
USalign batches 1–4 remain validated complete; batches 5–8 remain active with
growing pipeline logs. KUACC refinement advanced to shard-21 task `_30`; later
tasks and shards remain association-limited, with the controller still
dependency-held.

## Live checkpoint — 2026-09-28T08:09:43+03:00

GTalign advanced to `121 completed`, `4 running`, and `2 explicit failures`;
zero completed records are missing required refinement/DockQ stages. USalign
batches 1–4 have validated completion markers; batches 5–8 remain active
without exit markers, with pipeline logs increasing. KUACC shard-21 reached
active task `_73`; remaining tasks/shards are still association-limited and
the controller remains dependency-held.

## Live checkpoint — 2026-09-28T08:11:01+03:00

GTalign remains `121 completed`, `4 running`, and `2 explicit failures`; no
completed record has lost required refinement/DockQ stages. USalign batches 1–4
remain validated complete; batches 5–8 remain active and their logs grew since
the previous poll. KUACC shard-21 advanced to active task `_116`; later
tasks/shards remain association-limited and the controller remains
dependency-held.

## Live checkpoint — 2026-09-28T08:12:06+03:00

GTalign remains `121 completed`, `4 running`, and `2 explicit failures`;
active worker logs are still changing and no completed checkpoint is missing
required stages. USalign batches 1–4 remain validated complete; batches 5–8
remain active without exit markers. KUACC shard-21 has progressed through task
`_142` (task `_119` is completing); later tasks/shards remain
association-limited and the controller remains dependency-held.

## Live checkpoint — 2026-09-28T08:12:53+03:00

GTalign remains `121 completed`, `4 running`, and `2 explicit failures`; the
four active worker logs remain live and completed-stage validation is clean.
USalign batches 1–4 remain complete and validated; batches 5–8 remain active
without exit markers. KUACC shard-21 advanced through active task `_178`;
later tasks/shards remain association-limited and the controller remains
dependency-held.

## Live checkpoint — 2026-09-28T08:16:53+03:00

GTalign audit now reports `124 completed`, `4 running`, and `3 explicit
failures`; completed records remain stage-complete. The new failure is
`medium_1ijk_021`, where source groups have 2+1 chains but the native
receptor/ligand groups have 1+2 chains. The selected-candidate audit contains
268 rows with this cross-partition cardinality pattern. The current worker
incorrectly requires per-side source/native cardinalities to match; the total
chain count is still bijective for this class. The two `medium_1wq1_045`
failures remain true total-chain-count mismatches and are retained explicitly.
USalign batches 1–4 remain validated complete and batches 5–8 remain active;
KUACC refinement remains active and no downstream compact marker exists.

## Live checkpoint — 2026-09-28T08:20:22+03:00

GTalign has `124 completed`, `4 running`, and `3 explicit failures`; the newly
observed `medium_1ijk_021` failure confirmed the cross-partition cardinality
defect. The retained repair wrapper passes a real-candidate smoke test
(`ABC:ABC` model/native mapping) and still rejects the known `1wq1`
total-chain mismatch. The active array remains on the original worker for
provenance; the repair lane will be submitted only after the original array
has finished and its failed-row set is final. USalign batches 1–4 remain
validated complete and batches 5–8 remain active. KUACC refinement remains
active and no downstream compact marker exists; cleanup remains ineligible.

## Live checkpoint — 2026-09-28T08:21:32+03:00

GTalign has `132 completed`, `4 running`, and `4 explicit failures`; the
failures are one cross-partition `medium_1ijk_021` record and three
`medium_1wq1_045` total-chain mismatches. No completed record is missing a
required refinement/DockQ stage. The corrected repair wrapper remains
smoke-tested but is not submitted while the original array is active. USalign
batches 1–4 remain validated complete; batches 5–8 remain active. KUACC
refinement remains active and dependency-held downstream stages remain
unchanged; cleanup is still ineligible.

## Live checkpoint — 2026-09-28T08:24:41+03:00

GTalign has `134 completed`, `4 running`, and `4 explicit failures`; the four
failures remain one repairable cross-partition record and three explicit
`1wq1` total-chain mismatches. The compact 268-row repair manifest is prepared
and hashed at `gtalign_common_refinement_19855/repair/manifest.json`; it is
not submitted while array `1709167` remains active. USalign batches 1–4 remain
validated complete and batches 5–8 remain active. All downstream
aggregate/compact jobs remain dependency-held; cleanup is still ineligible.

## Live checkpoint — 2026-09-28T08:25:31+03:00

GTalign has `135 completed`, `4 running`, and `4 explicit failures`; all
completed checkpoints still contain the required refinement/DockQ stages. The
original array `1709167` remains active, so the hashed 268-row repair manifest
remains prepared but unsubmitted. USalign batches 1–4 remain validated
complete; batches 5–8 remain active. KUACC shard-21 has advanced through task
`_543`; downstream jobs remain dependency-held and cleanup remains ineligible.

## Live checkpoint — 2026-09-28T08:26:36+03:00

GTalign has `136 completed`, `4 running`, and `4 explicit failures`; all
completed checkpoints remain complete-stage valid. Array `1709167` is still
active, so the 268-row repair manifest remains unsubmitted. USalign batches 1–4
remain validated complete; batches 5–8 remain active. KUACC shard-21 advanced
through task `_581`; downstream compact jobs remain dependency-held and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:29:25+03:00

GTalign has `137 completed`, `4 running`, and `4 explicit failures`; all
completed checkpoints remain complete-stage valid. Array `1709167` is still
active, so the hashed 268-row repair manifest remains unsubmitted. USalign
batches 1–4 remain validated complete; batches 5–8 remain active. KUACC
refinement remains active and shard-21 has advanced through task `_616`;
downstream compact/aggregate jobs remain dependency-held and cleanup remains
ineligible.

## Live checkpoint — 2026-09-28T08:30:11+03:00

GTalign has `139 completed`, `4 running`, and `4 explicit failures`; all
completed checkpoints remain complete-stage valid. Array `1709167` is still
active, so the hashed 268-row repair manifest remains unsubmitted. USalign
batches 1–4 remain validated complete; batches 5–8 remain active. KUACC
refinement remains active and shard-21 has advanced through task `_688`;
downstream compact/aggregate jobs remain dependency-held and cleanup remains
ineligible.

## Live checkpoint — 2026-09-28T08:32:19+03:00

GTalign has `141 completed`, `4 running`, and `4 explicit failures`; all
completed checkpoints contain the required input-normalization, refinement, and
DockQ stages. Array `1709167` remains active, so the hashed 268-row
cross-partition repair manifest remains unsubmitted. USalign batches 1–4 remain
validated complete; batches 5–8 remain active. KUACC refinement and downstream
dependency-held jobs remain active/held; cleanup remains ineligible.

## Live checkpoint — 2026-09-28T08:34:11+03:00

GTalign has `145 completed`, `3 running`, and `5 explicit failures`; all
completed checkpoints contain the required input-normalization, refinement, and
DockQ stages. Failures are one repairable equal-total cross-partition
`medium_1ijk_021` candidate and four explicit `medium_1wq1_045` total-chain
mismatches. Array `1709167` remains active, so the hashed 268-row repair
manifest remains unsubmitted. USalign batches 1–4 remain validated complete;
batches 5–8 remain active. KUACC refinement/downstream jobs remain active or
dependency-held; cleanup remains ineligible.

## Live checkpoint — 2026-09-28T08:35:13+03:00

GTalign has `146 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the required nested-stage audit. The failure
classification is unchanged: one repairable equal-total cross-partition
candidate and four explicit `medium_1wq1_045` total-chain mismatches. Array
`1709167` remains active, so repair submission remains held. USalign batches 1–4
remain validated complete; batches 5–8 remain active. KUACC refinement and
downstream jobs remain active or dependency-held; cleanup remains ineligible.

## Live checkpoint — 2026-09-28T08:36:59+03:00

GTalign remains at `146 completed`, `4 running`, and `5 explicit failures`, with
all completed checkpoints passing the nested-stage audit. The bounded
three-shard cross-partition repair is submitted as `1709749`, dependency-held
by `afterany:1709167:1709178` so it cannot overlap the original owner or first
aggregate. Submission provenance is retained in
`gtalign_common_refinement_19855/repair/submission.json`. USalign batches 1–4
remain validated, batches 5–8 remain active, and cleanup remains ineligible.

## Live checkpoint — 2026-09-28T08:38:03+03:00

GTalign has `150 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749` is
still dependency-held by `1709167` and `1709178`, with no repair output yet.
The original refinement array remains active. USalign batches 1–4 remain
validated, batches 5–8 remain active, and cleanup remains ineligible.

## Live checkpoint — 2026-09-28T08:38:53+03:00

GTalign has `151 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:39:30+03:00

GTalign has `152 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:40:07+03:00

GTalign has `154 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:40:46+03:00

GTalign has `155 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:41:52+03:00

GTalign has `157 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:42:30+03:00

GTalign has `158 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:52:48+03:00

GTalign has `159 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:53:54+03:00

GTalign has `160 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T08:54:58+03:00

GTalign has `161 completed`, `4 running`, and `5 explicit failures`; all
completed checkpoints pass the nested-stage audit. Repair array `1709749`
remains dependency-held by `1709167` and `1709178`, with no repair output yet.
USalign batches 1–4 remain validated, batches 5–8 remain active, and cleanup
remains ineligible.

## Live checkpoint — 2026-09-28T09:00:05+03:00

GTalign has `166 completed`, `4 running`, and `5 explicit failures`; all 166
completed checkpoint records pass the corrected nested-stage audit. Repair array
`1709749` remains dependency-held by `1709167` and `1709178`, with no repair
output yet. USalign batches 1–4 remain validated, batches 5–8 remain active,
and cleanup remains ineligible.

## Live checkpoint — 2026-09-28T09:01:36+03:00

GTalign has `167 completed`, `4 running`, and `5 explicit failures`; all 167
completed checkpoint records pass the corrected nested-stage audit. Repair array
`1709749` remains dependency-held by `1709167` and `1709178`, with no repair
output yet. USalign batches 1–4 remain validated, batches 5–8 remain active,
and cleanup remains ineligible.

## Live checkpoint — 2026-09-28T09:03:13+03:00

GTalign has `169 completed`, `4 running`, and `5 explicit failures`; all 169
completed checkpoint records pass the corrected nested-stage audit. Repair array
`1709749` remains dependency-held by `1709167` and `1709178`, with no repair
output yet. USalign batches 1–4 remain validated, batches 5–8 remain active,
and cleanup remains ineligible.

## Live checkpoint — 2026-09-28T09:04:01+03:00

GTalign has `170 completed`, `4 running`, and `5 explicit failures`; all 170
completed checkpoint records pass the corrected nested-stage audit. Repair array
`1709749` remains dependency-held by `1709167` and `1709178`, with no repair
output yet. USalign batches 1–4 remain validated, batches 5–8 remain active,
and cleanup remains ineligible.
