# PRISM large refinement processing — 2026-09-27

## Bounded objective

Process the existing generated BM55 transformed candidates for the TMalign and
MultiProt arms through the already validated input-normalization, FiberDock,
external Rosetta, and DockQ stages. Do not recompute alignment, run PRODIGY,
modify raw benchmark outputs, or include USalign (zero transformed pairs).

## Frozen scope

- TMalign generated candidates: 3,068.
- MultiProt generated candidates: 66,827.
- Total eligible generated candidates: 69,895.
- Explicitly excluded clash-rejected rows: 15,103 TMalign and 3,277
  MultiProt; these never passed the transformation gate.
- Source: the 2026-09-27 compact BM55 candidate table and case manifests on
  KUACC; all inputs are referenced by absolute path and remain read-only.
- The run root is a new namespace:
  `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-large-refine-20260927`.

## Execution contract

Each candidate has an atomic checkpoint and independent result directory. The
worker records stage status, elapsed time, input/output hashes, scheduler IDs,
refined models, raw DockQ JSON, and explicit failure reasons. A candidate
stage failure does not stop other stages or array elements; an incomplete or
stage-failed checkpoint is resumable and remains visible to reconciliation.

FiberDock uses the shared read-only bundle but a candidate-specific output
prefix, preventing the fixed-prefix collision found in the earlier worker.
Manifest shards are limited to 1,000 rows because KUACC reports
`MaxArraySize=1001`; each Slurm task reads only its shard rather than reparsing
the full manifest.

## Placement and continuation

KUACC `mid` / account `users` is the primary backend. The explicit `mid` QOS
was rejected by Slurm during test-only validation and is not used. Arrays are
submitted in bounded waves so the reported `MaxJobCount=10000` is not exceeded.
The next wave is launched only after the prior wave reaches a terminal state;
if a wave cannot be submitted because of a live per-user limit, unstarted
shards are left in the manifest for backend-qualified VALAR fallback rather
than duplicated.

## Acceptance gates

1. Test-only submission succeeds with the live partition/account/resources.
2. Two-candidate mixed-pipeline pilot proves disjoint FiberDock outputs and
   resumable stage checkpoints.
3. DockQ is executed through the existing KUACC `prism_dockq_20260919`
   environment and returns raw JSON with `GlobalDockQ`.
4. Large aggregation reports one row per manifest index, including failed,
   skipped, and incomplete states; no scientific summary is updated before
   reconciliation is complete.

## DockQ chain contract

The native/reference PDB is the ground truth for chain labels and the
receptor/ligand interface. A one-chain receptor and one-chain ligand is a
valid two-partner DockQ case; a genuinely monomeric reference has no
protein-protein interface and is recorded as `valid_unscored`. Missing,
overlapping, or otherwise inconsistent native chain groups are recorded as
`invalid_mapping`, not silently reduced to a single-chain score. Unexpected
DockQ exceptions are recorded as `scorer_error` with the exception class.

The worker now performs this preflight before canonicalization and scoring,
and the aggregate table preserves the available/missing native chains,
mapping status, scoring scope, and error class. Existing checkpoints are not
rewritten by this code change; affected candidates must be resumed or rerun
in a new controlled attempt after the updated worker is staged.

The BM55 `T_*.csv` audit covers all 257 cases and 69,895 selected candidates:
all benchmark `Complex` mappings agree after stripping trailing annotation
markers. The assembled native PDB inventory contains 233 exact chain sets,
20 cases with extra blank/auxiliary chains but all expected chains present,
and four cases with genuinely missing expected chains (`3aad`, `1oyv`, and the
two `3p57` alternate complexes). The four missing-coordinate cases remain
excluded from DockQ until their native assemblies are recovered.

## Resource policy

Baseline task request: one CPU and 8 GB RAM, with the partition’s one-day
maximum walltime. This matches the successful 40-candidate run and leaves
parallelism to the scheduler. The plan does not request GPUs: FiberDock,
Rosetta, and DockQ are CPU workloads.

## Risks and limits

- Scratch has approximately 17 TB free but is 91% utilized; results are kept
  in the new root and no cleanup is authorized by this plan.
- The current worker is candidate-serial across FiberDock, Rosetta, and DockQ;
  candidate arrays provide parallelism while preserving stage-level resume.
- Scheduler disappearance is not treated as success; checkpoint and artifact
  reconciliation remains authoritative.

## GlobalDockQ multimer re-aggregation — 2026-09-27

The large-run aggregation path was corrected to match the already validated
small-comparison contract. For each DockQ stage, `*_dockq` now means the
bounded `GlobalDockQ` parsed from the preserved raw JSON. The former
checkpoint value is retained as `*_dockq_legacy`, and the official
interface-sum `best_dockq` is retained as `*_dockq_sum`. Missing or invalid
GlobalDockQ is explicitly unscored rather than falling back to a sum.

The authoritative live snapshot was written to the new namespace
`/scratch/users/rshadi25/valar-remote-runs/prism-prescript-large-refine-20260927-globaldockq-v3`.
It contains 8,000 of 69,895 expected checkpoint rows while the source run was
still active. All 7,431 FiberDock and 2,315 Rosetta primary scores parsed in
that snapshot were within [0, 1]. Interface sums exceeded 1 for 248 FiberDock
and 232 Rosetta rows; the out-of-range-sum means were 2.948805 and 1.919420,
respectively. This is expected multimer behavior and is now retained
only as diagnostic provenance. See
`evidence/prism-prescript/globaldockq-reaggregation-20260927/README.md` for
the exact snapshot metrics and hash.
