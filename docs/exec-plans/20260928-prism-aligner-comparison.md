# PRISM aligner comparison — 2026-09-28 execution plan

## Purpose / Big Picture

Complete an evidence-backed comparison of TMalign, MultiProt, USalign, and
GTalign within the existing PRISM pipeline.  Alignment is the interchangeable
provider; transformation, filtering, ranking, refinement, and corrected DockQ
remain common wherever the contract permits.  The durable scientific result is
a compact, reproducible package, not a copy of the KUACC working tree.

The selected project is `/scratch/rshadi25/GitHub/PRISM-prescript`.  Existing
PRISM execution evidence is read from the run-scoped roots below; no sibling
source tree is modified.

## Progress

- Recovered project memory, framework status, plans, and the existing Graphify
  graph; no graph rebuild is needed.
- TMalign and MultiProt have validated compact alignment/candidate ledgers for
  257 cases on the historical 19,855-template panel.
- GTalign has a reusable historical full search/transformation result, but its
  panel and downstream scoring require reconciliation before matched claims.
- The USalign parser contract now preserves `tm_score_query` and
  `tm_score_ref`, records `tm_score_contract`, and uses Structure_2/reference
  normalization for the shared `tm_score` field.  The previous `max(score1,
  score2)` collapse is no longer used in this parser.
- Focused parser and contract tests pass: 15 tests.
- Identical-pair TMalign/USalign smoke passes: both parsers report 28 matches,
  equal mappings/transforms, RMSD 0.00, both normalized TM-scores 1.0, and
  byte-identical common transformation output.
- The existing 946-template/10-case run was retrieved as compact pilot
  evidence.  It has useful timing/counts but does not record the required
  1/2/4/8/16 worker sweep and its USalign-labelled command/progress needs
  provenance repair before scientific use.
- Exact current-panel evidence has now been copied into the dated compact
  package: 19,948 checked, 19,062 calculated, 19,058 materialized, and the
  historical 19,855 hash/panel lane. Exact alignment-only jobs and their
  limitations are recorded in the stage ledger.
- A reusable USalign provider branch now routes through the existing PRISM
  `align()` → transformation/filtering pipeline, with explicit `USalign`
  provenance, `-fast` support, and run-scoped `alignment_usalign` output.
- The first live worker-sweep wrapper attempt (1708915) was invalidated after
  detecting serial dispatch; its reason is retained and raw partial records
  were cleaned. The corrected bounded sweep job 1708928 completed all 10
  configurations with zero execution failures. `-fast` changed 1,385/1,892
  mappings/scores relative to default, so default/16 is the selected
  production configuration.
- A live KUACC external-Rosetta refinement run remains active through
  continuation arrays. Its jobs and run-scoped artifacts are not cancelled or
  cleaned until their downstream compact results are validated.
- Existing GTalign full-BM55 transformed scoring is now reconciled as a
  reusable historical-panel baseline: 15,440 rows, 14,300 scored, and 1,140
  score failures over 216 cases, with GlobalDockQ kept distinct from the
  diagnostic best/interface score. It is not pooled with the exact-panel
  alignment-only lane or treated as refined.
- The common transformed scorer's retention gate and the USalign compactor's
  case-identity join have test coverage. A read-only corrected refinement
  aggregator is staged in the active KUACC run; it computes paired deltas
  only when both compact global scores exist and does not delete files.
- A resumable historical-panel USalign provider array 1708992 is running on
  VALAR `kutem` with default USalign and 16 workers; its compaction array
  1709007 is dependency-queued.  A corrected transformed-DockQ array 1709046
  is dependency-queued after compaction and uses the existing bijective scorer.
- The transformed scorer retention contract was corrected before execution:
  job 1709046 retains scored transformed halves for the pending common
  refinement consumer and deletes only combined/scoring scratch. The focused
  continuation suite passes 34 tests, with the retention fixture included.
- A compact `aggregate_matched_comparison.py` reducer and two regression tests
  now emit candidate, case-level, and top-k diagnostic tables from compact
  inputs; missing DockQ remains empty and oracle top-k is explicitly labelled
  diagnostic.
- GTalign common-refinement array 1709167 is executing the shared
  FiberDock/external-Rosetta/corrected-DockQ worker; the latest artifact poll
  found 92 checkpoint records: 87 top-level completed, 4 running, and 1
  explicit input-normalization failure. Its read-only aggregate 1709178 is
  dependency-gated and no cleanup is eligible.
- USalign batch aggregation is now a tested, read-only stage.  Job 1709192
  (test-only 1709191) is queued after corrected transformed DockQ 1709046;
  it will validate all 26 compact batch packages before producing the single
  downstream input for common refinement.  The tested USalign refinement
  manifest generator records selected/rejected rows without copying files.
- Added `aggregate_final_matched_comparison.py` with a regression test. It
  joins exact candidate identities across transformed and refined tables,
  emits method-by-case/candidate/ranking/paired-delta tables, and reports a
  normal-approximation 95% CI only for exact transformed/refined pairs.
- Added the dependency-gated, fail-closed USalign refinement handoff job
  1709233 (dry-run 1709232). It will prepare, but deliberately not submit,
  the common-refinement manifest after validated 26-batch compaction.
- Added `aggregate_timing_resources.py` with a regression test. It will
  compact stage wall time, CPU-core-hours, GPU-hours, and status counts from
  existing timing/refinement ledgers without retaining raw work directories.
- Added the fail-closed `cleanup_usalign_run.py` gate and regression test. It
  requires validated USalign compaction, a validated refinement handoff and
  common-refinement aggregate, and a validated final package before deleting
  only the run-scoped batch/refinement work directories. No cleanup has been
  applied.
- Revalidated the comparison-focused regression set after the cleanup and
  provenance changes: 74 tests passed. The DockQ runtime fixture checks the
  explicit `--n_cpu 1` contract while still verifying the mapping argument.
- A bounded audit found explicit `alignment_unavailable` USalign records. The
  compact alignment reducer now excludes non-success rows from numerical
  means, retains status/failure counts and valid-record counts, and records a
  source-snapshot delta. The affected comparison tests pass 35/35; the active
  provider does not require rerun.
- Before final inputs became available, the matched reducer was strengthened
  to report no-ranking/all, deterministic, PRODIGY-when-available, top-1/3/5,
  and diagnostic-oracle strategies, plus transformed/refined coverage and
  cross-interface quality and method-by-split summary tables. Its focused
  tests pass 3/3; the reducer has not yet consumed production tables.
- Latest execution poll: USalign provider tasks 1--6 remain active or newly
  active with raw batches and dependency stages pending; GTalign is at 92
  checkpoints (87 completed, 4 running, 1 failed). KUACC child 3140177 and controller
  3132296 completed with exit 0:0 and scheduled continuation arrays
  3140365/3140401--3140407; that run has 16,651 checkpoints and no aggregate.

## Surprises & Discoveries

- The durable BM55 source run is a 19,855-template historical/working panel,
  not the 19,948 checked / 19,062 calculated / 19,058 materialized panel
  requested for explicit tracking.  These lanes must not be pooled silently.
- TMalign and MultiProt alignment records are complete, but their candidate
  yields differ sharply: 3,068 and 66,827 generated candidates respectively;
  MultiProt's historical RMSD proxy is not a TM-score.
- The previous USalign full array did not constitute a completed comparison:
  one case failed and 256 were not started.  A bounded correctness test and
  946-template worker pilot are therefore required before production.
- Valid GTalign alignment output exists, so GPU search should not be repeated
  merely to regenerate downstream tables.

## Decision Log

1. Use the 19,855-template compact run as a clearly labelled historical lane
   while locating/freezing the 19,948/19,062/19,058 manifest lane.
2. Reuse validated TMalign/MultiProt and GTalign alignment evidence; rerun only
   the first invalid downstream stage.
3. Preserve both TM-score normalizations and use the frozen PRISM orientation
   (`Structure_1=query`, `Structure_2=reference/template`) for the common score.
4. Use corrected DockQ records with separate GlobalDockQ and requested
   cross-interface DockQ; failure and not-scoreable are statuses, never zero.
5. Keep current KUACC jobs untouched.  Submit new CPU work only after live
   eligibility and `sbatch --test-only` validation, and only in resumable,
   case-wise compacting batches.
6. Cleanup is limited to positively identified, run-scoped PRISM intermediates
   after compact result, provenance, checksum, and downstream-consumer checks.

## Outcomes & Retrospective

The outcome must include one case-level row per method and BM5.5 case, useful
candidate-level rows, corrected transformed/refined DockQ records, overlap and
ranking analyses, paired refinement deltas, timing/resource tables, failure
catalogues, hashes, commands, and a README.  Historical/panel differences,
execution facts, interpretations, and unresolved limitations remain separate.

## Context and Orientation

Primary selected checkout:

`/scratch/rshadi25/GitHub/PRISM-prescript`

Historical compact source run:

`/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919`

Compact package:

`/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919/compact_results_20260927`

Active refinement run:

`/scratch/users/rshadi25/valar-remote-runs/prism-prescript-large-refine-20260927`

The historical source manifest records 257 pairs, 19,855 templates, expected
79,420 records per case and 20,410,940 records per aligner.  The compact
validation package is the authority for completed TMalign/MultiProt counts.

Live scheduler evidence must be refreshed immediately before any new job:
VALAR uses `/home/rshadi25/gpu-status.sh`, `squeue`, `sinfo`, and `scontrol`;
KUACC uses explicit BatchMode SSH and `squeue`/`sinfo`/`scontrol` without a
`sacct` dependency.

## Plan of Work

1. Reconcile the compact historical lane and freeze the exact comparison
   manifest/panel lanes, including the three requested template counts.
2. Complete DockQ contract reconciliation and produce a validated compact
   scoring audit for reusable TMalign, MultiProt, and GTalign candidates.
3. Reuse TMalign/MultiProt ledgers; retain explicit failures and rejection
   reasons; continue only missing downstream stages.
4. Validate the USalign parser on identical real input pairs, then benchmark
   946 representative templates at worker counts 1/2/4/8/16 and default/
   `-fast` configurations subject to resources.  Reuse existing 946-template
   counts/timings only as a bounded historical reference, not as the sweep.
5. If the corrected pilot passes, submit measured USalign production batches
   with the common transformation/filtering provider, no-refine alignment
   stage first, resumable batch manifests, and streaming compaction.
6. Reconcile existing GTalign output and skip GPU rerun unless a required
   alignment artifact is genuinely absent or incompatible.
7. Score transformed candidates with the common evaluator before refinement;
   retain transformed structures until all common-refinement consumers finish,
   then run common external-Rosetta refinement only for the selected matched
   set.
8. Produce top-k/ranking, overlap/disagreement, timing, and paired refinement
   analyses; validate and compact each stage before cleanup.
9. Update durable project/framework status and report implementation, testing,
   validation, review, scientific interpretation, and limitations separately.

## Concrete Steps

### Reconciliation and compact evidence

- Retrieve only declared compact files (`validation.json`, candidate validation,
  comparison/case/timing/failure tables, manifests, schema, and cleanup
  manifest) from the KUACC source run.
- Build/update the stage ledger at
  `benchmark/prism_processed_results/prism_aligner_comparison_20260928/`.
- Locate and hash the 19,948/19,062/19,058 manifest evidence; if absent,
  record the missing source rather than infer counts from the 19,855 lane.

### Correctness and execution

- Run identical-pair TMalign/USalign checks for mapping, aligned length, RMSD,
  matrix, both TM scores, and downstream transformation.
- Inspect live resource/QOS/storage state and use `sbatch --test-only` before
  any bounded pilot or production submission.
- Monitor job state plus persistent logs and expected artifacts.  Submission or
  disappearance alone never establishes success.

### Stage-wise retention

For every stage: process, validate counts/metrics, write a compact result,
record command/version/hash/job metadata, then delete only explicitly listed
disposable run-scoped files whose consumers have finished.  Preserve failed
statuses and reasons before deleting failed working artifacts.

## Validation and Acceptance

- Parser tests cover Structure_1/Structure_2 and reject incomplete USalign
  contracts.
- Every retained case has an explicit status: scored, valid_unscored,
  scored_cross_only, not_scoreable, score_failed, or the appropriate alignment
  failure/no-output state.
- TMalign/MultiProt/GTalign reusable artifacts have provenance and panel labels;
  USalign production is not accepted without the correctness smoke and pilot.
- Final tables pass row/count/hash/schema checks and distinguish transformed
  from refined and paired from unpaired observations.
- No final claim is based solely on pooled mean DockQ or candidate count.
- Retained package is sufficient to reconstruct inputs, parameters, software,
  commands, outputs, failures, timings, and cleanup decisions.

## Idempotence and Recovery

The 2026-09-28T06:12:29+03:00 poll found USalign task 1708992_1 still
validated at exit 0:0 with `available=10 unavailable=0` and empty stderr;
tasks 3, 5, 6, and 7 remain active with later tasks pending, and no compact
marker exists. GTalign has 96 checkpoint files (91 completed, 4 running, 1
explicit failure). The
dependency chain remains intact and no cleanup is eligible.

Case manifests and checkpoints are authoritative.  A rerun skips validated
case/stage outputs, reprocesses only missing/invalid descendants, and writes to
an explicit backend-qualified run namespace.  If a parser fails, reparse raw
output where retained; rerun alignments only when raw evidence is unavailable
or invalid.  Active unrelated/protected jobs are never modified.

## Artifacts and Notes

The stage ledger, compact tables, manifests, hashes, cleanup manifest, logs,
and validation reports belong under the dated comparison directory.  Large raw
alignment, transformed PDB, DockQ, and refinement trees remain temporary
run-scoped artifacts and are eligible for cleanup only after the corresponding
compact evidence is validated.

## Interfaces and Dependencies

The shared interfaces are `src/alignment.py`, common transformation/filtering,
`src/eval/dockq.py`, ranking/EDA modules, existing Slurm wrappers, TMalign,
USalign, MultiProt, GTalign, and the verified DockQ 2.1.3 environment.  The
USalign binary is the validated 20241108 executable reached through the
PRISM-compatible wrapper; the external tool repository is not modified.

## Live checkpoint — 2026-09-28T06:22:15+03:00

GTalign common refinement array `1709167` now has 97 checkpoint files (92
completed, 4 running, 1 explicit failure). USalign production tasks
`1708992_3`, `1708992_5`, `1708992_6`, and `1708992_7` remain active; every
downstream dependency is still held. This is execution evidence only: no
compact, corrected-DockQ, final aggregation, or cleanup gate has passed.

## Live checkpoint — 2026-09-28T06:25:50+03:00

KUACC common refinement is making scheduler-visible progress: continuation
array `3140365` reached running task `_919`, and continuation `3140401` has
running tasks `_786/_787`. Other continuation shards remain association-limit
pending. A large-tree checkpoint count probe timed out and was not treated as
terminal evidence.

## Live checkpoint — 2026-09-28T06:27:19+03:00

GTalign common refinement has 99 checkpoint files (94 completed, 4 running,
1 explicit failure). KUACC continuation `3140401` advanced to task `_862`,
with continuation `3140365` still active at `_919`. No local dependency has
opened and no cleanup gate has passed.

## Live checkpoint — 2026-09-28T06:28:26+03:00

GTalign common refinement has 100 checkpoint files (95 completed, 4 running,
1 explicit failure). KUACC continuation `3140401` reached `_900`, while
`3140365` remains active at `_919`. No local dependency has opened.

## Live checkpoint — 2026-09-28T06:29:06+03:00

KUACC continuation `3140401` advanced to `_919` and continuation `3140365`
remains active at `_919`. GTalign remains at 100 checkpoint files (95
completed, 4 running, 1 explicit failure); no local dependency has opened.

## Live checkpoint — 2026-09-28T06:30:32+03:00

KUACC continuation `3140401` reached `_964`, while `3140365` remains active
at `_919`. GTalign remains at 100 checkpoints (95 completed, 4 running, 1
explicit failure); no local compact, DockQ, aggregation, or cleanup marker
exists.

## Live checkpoint — 2026-09-28T06:31:24+03:00

KUACC continuation `3140401` reached `_988`, while `3140365` remains active at
`_919`. Local USalign and GTalign workers remain active; all downstream jobs
remain dependency-held. The broad marker scan exceeded the bounded probe and
was not interpreted as terminal evidence.

## Live checkpoint — 2026-09-28T06:32:05+03:00

KUACC continuation `3140401` reached `_997`–`_999` and shard `3140402` began
tasks `0`–`10`; remaining tasks are association-limit pending. Local USalign
and GTalign workers remain active and downstream stages remain dependency-held.

## Live checkpoint — 2026-09-28T06:34:29+03:00

KUACC shard `3140402` reached running task `_93`; tasks `_94–999` remain
association-limit pending. Local USalign/GTalign workers remain active and
downstream stages remain held; GTalign remains at 95 completed, 4 running, and
1 explicit failure.

## Live checkpoint — 2026-09-28T06:42:43+03:00

KUACC shard `3140402` reached running task `_369`; tasks `_370–999`
remain association-limit pending. USalign task `1708992_8` and all GTalign
workers remain active; GTalign remains at 95 completed, 4 running, and 1
explicit failure.

## Live checkpoint — 2026-09-28T06:46:49+03:00

KUACC shard `3140402` reached running task `_511`; tasks `_512–999`
remain association-limit pending. GTalign remains at 96 completed, 4 running,
and 1 explicit failure; USalign task `1708992_8` remains active and
downstream local stages remain held.

## Live checkpoint — 2026-09-28T06:44:45+03:00

GTalign common refinement has 101 checkpoints (96 completed, 4 running, 1
explicit failure). KUACC shard `3140402` reached `_440`; tasks `_441–999`
remain association-limit pending. Local downstream stages remain held.

## Live checkpoint — 2026-09-28T06:44:02+03:00

KUACC shard `3140402` reached running task `_415`; tasks `_416–999`
remain association-limit pending. USalign task `1708992_8` and GTalign
workers remain active; GTalign remains at 95 completed, 4 running, and 1
explicit failure.

## Live checkpoint — 2026-09-28T06:43:21+03:00

KUACC shard `3140402` reached running task `_390`; tasks `_391–999`
remain association-limit pending. USalign task `1708992_8` and GTalign
workers remain active; GTalign remains at 95 completed, 4 running, and 1
explicit failure.

## Live checkpoint — 2026-09-28T06:35:30+03:00

KUACC shard `3140402` reached running task `_122`; tasks `_123–999` remain
association-limit pending. Local provider/refinement arrays remain active;
GTalign remains at 95 completed, 4 running, and 1 explicit failure. No
dependency or cleanup gate has opened.

## Live checkpoint — 2026-09-28T06:32:41+03:00

KUACC shard `3140402` reached running tasks through `_31`; tasks `_32–999`
remain association-limit pending. Local USalign/GTalign workers remain active
and downstream stages remain dependency-held; GTalign remains at 95 completed,
4 running, and 1 explicit failure.

## Live checkpoint — 2026-09-28T06:33:43+03:00

KUACC shard `3140402` reached running task `_65`; tasks `_66–999` remain
association-limit pending. Local USalign/GTalign workers remain active and
downstream stages remain held; GTalign remains at 95 completed, 4 running, and
1 explicit failure.

## Live checkpoint — 2026-09-28T06:36:08+03:00

KUACC shard `3140402` reached running task `_146`; tasks `_147–999` remain
association-limit pending. Local USalign/GTalign workers remain active;
GTalign remains at 95 completed, 4 running, and 1 explicit failure.

## Live checkpoint — 2026-09-28T06:37:59+03:00

KUACC shard `3140402` reached running task `_208`; tasks `_209–999` remain
association-limit pending. Local USalign/GTalign workers remain active;
GTalign remains at 95 completed, 4 running, and 1 explicit failure. No
downstream gate has opened.

## Live checkpoint — 2026-09-28T06:38:43+03:00

USalign provider task `1708992_8` started; tasks `9–26` remain array-limit
pending. KUACC shard `3140402` reached `_234`. GTalign remains at 95 completed,
4 running, and 1 explicit failure; local downstream stages remain held.

## Live checkpoint — 2026-09-28T06:39:34+03:00

USalign task `1708992_8` remains active and tasks `9–26` remain array-limit
pending. KUACC shard `3140402` reached `_262`; later tasks remain
association-limited. GTalign remains at 95 completed, 4 running, and 1
explicit failure.

## Live checkpoint — 2026-09-28T06:40:31+03:00

USalign task `1708992_8` remains active; tasks `9–26` remain array-limit
pending. KUACC shard `3140402` reached `_288`; later tasks remain
association-limited. GTalign remains at 95 completed, 4 running, and 1
explicit failure.

## Live checkpoint — 2026-09-28T06:42:01+03:00

KUACC shard `3140402` reached running task `_341`; tasks `_342–999`
remain association-limit pending. USalign task `1708992_8` and all four
GTalign workers remain active; GTalign remains at 95 completed, 4 running, and
1 explicit failure.

## Live checkpoint — 2026-09-28T06:45:40+03:00

KUACC shard `3140402` reached running task `_471`; tasks `_472–999`
remain association-limit pending. GTalign remains at 101 checkpoints (96
completed, 4 running, 1 explicit failure). USalign task `1708992_8` remains
active and downstream local stages remain held.
