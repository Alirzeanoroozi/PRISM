# PRISM KUACC Compact Results and Retention

## Purpose / Big Picture

Create a small, reproducible evidence package from the completed PRISM BM55
KUACC case ledgers before any storage cleanup. The package must preserve
case-level outcomes, alignment/transformation totals, provenance hashes,
failure accounting, and a reviewable cleanup manifest without copying the
millions of raw alignment JSON files or transformed structures.

This is an evidence-preservation and reporting task. It does not authorize
deletion, cancellation, requeueing, promotion, or scientific production
reruns.

## Progress

- [x] Read the applicable execution, storage-hygiene, KUACC, and data-validation skills.
- [x] Recovered the BM55 submission manifest and durable per-case ledgers.
- [x] Confirmed TMalign and MultiProt have 257 durable case ledgers each; USalign has one failed and 256 not-started cases.
- [x] Confirmed the original aggregators timed out while copying raw artifacts; their status files are stale `running` markers.
- [x] Confirmed the current run has ranking, PRODIGY, and refinement disabled; no DockQ/PyRosetta results are present in this run.
- [x] Inspect candidate-level audit/score artifacts and select the compact-table grain.
- [x] Generate an isolated, idempotent case/stage compact-results package on KUACC.
- [x] Validate completeness, totals, hashes, duplicate keys, and status semantics for the case/stage package.
- [x] Complete streaming candidate-audit aggregation into retained-candidate and rejection-summary tables (KUACC job `3124073`, validated).
- [ ] Produce a dry-run cleanup manifest; obtain explicit authorization before any deletion.
- [ ] Decide whether raw/transformed artifacts must be retained for future DockQ/PyRosetta work.

## Surprises & Discoveries

- The existing `aggregate_cpu_case_batches.py` is not summary-only: it copies
  every raw alignment JSON and transformation artifact into a second output
  tree before writing its final summary. This is the direct operational reason
  the TMalign and MultiProt aggregators hit their one-day walltime.
- Slurm disappearance is not treated as success. Case ledgers and durable
  attempt summaries are the evidence of case completion; scheduler state is
  retained only as submission context.
- The current manifest requests `ranking=false`, `prodigy=false`, and
  `refinement=false`; therefore a compact table from this run cannot claim
  DockQ or PyRosetta evaluation.
- Candidate-audit JSONL is present per selected attempt. It contains one record
  per candidate/orientation and exposes threshold rejection reasons, match
  coverage, TM scores, clashes, contacts, and generated candidates. The
  follow-up collector retains non-threshold-rejected records individually and
  reduces threshold-rejected records to small distribution summaries.
- Remote GPU inventory is incomplete because the configured KUACC GPU-status
  script is absent. This run is CPU-only, so that does not change the present
  summary-only operation.

## Decision Log

1. Write summaries into a new run-scoped namespace under
   `.../compact_results_20260927/`; do not modify original case outputs.
2. Use standard-library CSV/JSON so the collector does not depend on pandas,
   pyarrow, or a particular KUACC environment.
3. Preserve one row per expected `(pipeline, case_index)` and represent
   missing, failed, not-started, and completed states explicitly.
4. Keep raw alignment and transformed artifacts classified as `REVIEW`, not
   `DELETE_CANDIDATE`, until the need for DockQ/PyRosetta is resolved.
5. No deletion is permitted in this phase. Any later deletion requires a
   reviewed manifest and a new explicit user authorization.

## Outcomes & Retrospective

To be completed after the compact package is generated and validated. Record
what was retained, what remains unresolved, and whether the package is enough
to regenerate the requested analyses without raw artifacts.

## Context and Orientation

- Project: `/home/rshadi25/valar-agent-framework`
- Remote host: `rshadi25@login.kuacc.ku.edu.tr`
- Run root: `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919`
- Submission manifest: `full_submission_manifest.json`
- Submitted arrays: TMalign `3121371`, MultiProt `3121373`, USalign `3121375`
- Submitted aggregators: TMalign `3121372`, MultiProt `3121374`, USalign `3121376`
- Expected cases: 257 per aligner; expected alignment records per aligner: 20,410,940
- Expected template manifest hash:
  `4680d3eda8030861a40373cd193b0e8bef7c21a771a90553c4c49814e964b48d`

## Plan of Work

1. Inspect candidate-level artifacts and the schema of representative case
   summaries.
2. Add/use a durable compact collector that reads only manifests, case
   ledgers, attempt summaries, and small provenance files.
3. Run it into the new compact namespace; emit case results, failure catalog,
   artifact inventory, schema, and cleanup manifest.
4. Validate row counts, expected-key coverage, status partitioning, arithmetic
   totals, template hashes, duplicate keys, and source-file hashes.
5. Retrieve the small summary package for review and update project memory or
   status only with evidence-backed claims.
6. Stop before deletion and present the retention choices and risks.

## Concrete Steps

All remote inspection and execution uses explicit non-interactive SSH. No
`sacct` is used. The collector must be idempotent: rerunning it replaces only
files inside the new compact output namespace, never source artifacts.

The collector records command/version/timestamp, input root, manifest hash,
source paths, row grain, and validation results. It must not copy raw JSON or
PDB files. Any candidate-level extraction that would require reading millions
of raw files is reported as a separate follow-up rather than silently done on
the login node.

## Validation and Acceptance

Acceptance requires:

- exactly 257 rows per aligner in the case table;
- unique `(pipeline, case_index)` keys and explicit status counts;
- expected-record totals equal `20,410,940` where the case ledgers are
  complete, with `successful + failed = alignment_records` per case;
- template count/hash and run configuration agree with the submission manifest;
- no completed row lacks its referenced task status and attempt summary;
- USalign's failed/not-started partition remains explicit;
- compact outputs contain no copied raw alignment JSON or transformed PDB;
- validation results and source hashes are written alongside the tables;
- `git diff --check` and relevant local script checks pass.

## Idempotence and Recovery

Use a new timestamped or task-dedicated output namespace. If collection is
interrupted, rerun against the same source root and output namespace after
checking that no source files changed. Preserve the previous validation JSON
if a rerun finds changed inputs, and create a new namespace rather than
overwriting evidence. A failed collector must leave its error and command
metadata in the compact namespace.

## Artifacts and Notes

Planned compact artifacts:

- `case_results.csv` and `case_results.jsonl`
- `failure_catalog.csv`
- `artifact_inventory.csv`
- `cleanup_manifest.csv`
- `schema.json`
- `validation.json`
- `collection_manifest.json`

The cleanup manifest is advisory until explicitly approved. It must classify
each path as `KEEP`, `REVIEW`, `PROTECTED`, or `DELETE_CANDIDATE`, include a
reason and reversibility/risk note, and distinguish exact sizes from estimates
or unmeasured directories.

## Interfaces and Dependencies

Inputs are the existing submission manifest, per-aligner corrected batch
manifests, per-case `task_status.json`, selected attempt `run_summary.json`,
and small scheduler/log metadata. The collector uses Python 3 standard
library only. Future DockQ/PyRosetta processing requires transformed models,
their provenance, and the relevant original alignment/mapping data; the
current compact case table alone is insufficient for that evaluation.
