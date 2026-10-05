# GlobalDockQ multimer re-aggregation — 2026-09-27

## Scope

This evidence records a read-only re-aggregation of the active KUACC
PRISM-prescript large-refinement run. It does not rerun alignment, refinement,
or DockQ and does not rewrite raw checkpoints or raw DockQ JSON.

- Source run: `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-large-refine-20260927`
- New output root: `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-large-refine-20260927-globaldockq-v3`
- Source manifest expectation: 69,895 candidates
- Snapshot: 8,000 checkpoint rows; the source run was still active
- Corrected table SHA-256: `5f877de746f9223b757ad6ec294efd37702cfe9aaa63080717c23a81ccd2da37`

Earlier `globaldockq-v1` and `globaldockq-v2` namespaces are preserved for
audit; v1 had a concurrent login-session retry, and v2 predates the final
summary-statistics fields. The `v3` namespace was run once after those
discoveries and is authoritative for this snapshot.

## Contract

For each multimer DockQ stage, the primary `*_dockq` column is the bounded
`GlobalDockQ` value from the preserved raw JSON. The interface sum is retained
as `*_dockq_sum` and the checkpoint's former value is retained as
`*_dockq_legacy`. A `*_dockq_contract` field distinguishes valid scores from
missing raw JSON and unscored/invalid results. No interface sum is used as a
fallback primary score.

## Snapshot results

| stage | GlobalDockQ count | mean | median | min | max | primary out of range | interface-sum count | sum mean | sum max | sum > 1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| FiberDock | 7,431 | 0.083766 | 0.008749 | 0.000013 | 0.911140 | 0 | 7,431 | 0.257321 | 13.392731 | 248 |
| Rosetta | 2,315 | 0.137307 | 0.020768 | 0.002134 | 0.856815 | 0 | 2,315 | 0.365458 | 13.214773 | 232 |

The values above 1 occur only in the preserved interface sums, which are
expected for multimers because the official DockQ output sums selected native
interface scores and separately provides `GlobalDockQ` as the normalized
multimer score. They are not valid primary model-quality scores.

Contract counts in this snapshot:

- FiberDock: 7,431 `GlobalDockQ`, 356 `missing_raw_json`, 213 `missing`.
- Rosetta: 2,315 `GlobalDockQ`, 5,472 `missing_raw_json`, 213 `missing`.

The mean of the out-of-range interface sums alone is 2.948805 for FiberDock
(248 rows; median 1.971656) and 1.919420 for Rosetta (232 rows; median
1.426393). These are diagnostic sums, not means of invalid normalized scores.

The source checkpoint snapshot was incomplete: 7,431 `completed`, 356
`completed_with_stage_failures`, 212 `failed`, and 1 `running`. The table must
be regenerated after the source run reaches its terminal reconciliation point
before scientific ranking or promotion claims are made.

## Validation

- Focused aggregation tests: 3 passed.
- Existing native-chain DockQ preflight tests: 5 passed.
- Python compile check and `git diff --check`: passed.
- Active scheduler jobs and source/raw artifacts: unchanged.
