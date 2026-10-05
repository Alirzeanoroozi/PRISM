# PRISM Refined-Model Evaluation Preparation

## Purpose / Big Picture

Prepare an auditable downstream path from the completed BM55 TMalign and
MultiProt transformed candidates to PyRosetta refinement and then DockQ/iRMSD
scoring. Do not rerun alignment and do not overwrite the completed run.

The user subsequently instructed that neither refinement nor DockQ should be
submitted yet: DockQ must consume refined-model outputs. This plan therefore
ends at validated preparation and records the runtime blocker rather than
submitting speculative jobs.

## Progress

- [x] Confirmed the current BM55 run has transformed but unrefined candidates.
- [x] Found 88,275 retained transformed candidates: 18,171 TMalign and 70,104 MultiProt.
- [x] Inspected existing refinement/scoring logic and identified its aggregate-root assumption.
- [x] Probed the KUACC PyRosetta bundle; import fails because `rosetta.so` requires `GLIBC_2.27` while KUACC is older.
- [x] Confirmed the VALAR login has glibc 2.28 but does not expose the KUACC run filesystem or bundle.
- [x] Added durable downstream manifest, DockQ batch, and PyRosetta preflight scripts; no jobs submitted.
- [x] Characterized query/transformed chain counts across all 257 cases.
- [ ] Stage and run the downstream manifest preparation after the refined-model backend is available.
- [ ] Submit PyRosetta refinement only on a validated compatible backend.
- [ ] Submit DockQ/iRMSD only after refinement completion and output validation.

## Surprises & Discoveries

- Existing `refine_and_score.py` assumes one completed aggregate run root;
  this BM55 run instead has per-case attempt roots because the aggregate jobs
  timed out while copying raw files.
- The PyRosetta bundle is present but is not runnable on KUACC due to glibc,
  and no PyRosetta module is exposed by KUACC's module system.
- Transformed sides are not uniformly single-chain. The 257 completed query
  cases contain 1–6 chains on the receptor side and 1–6 on the ligand side.
- The most common shapes are `(1,1)` in 144 cases and `(2,1)` in 78 cases;
  the remaining cases include `(1,2)`, `(2,2)`, `(2,3)`, `(2,6)`, `(3,1)`,
  `(4,1)`, `(4,2)`, `(4,4)`, and `(6,1)`.
- Source/model chain IDs are often renamed during transformation; validation
  must compare chain counts and use the dataset's native receptor/ligand chain
  mapping, not assume source IDs equal native IDs.
- Eighteen cases have a query-PDB chain-count mismatch against the native
  evaluator contract and must be reviewed/excluded before any DockQ scoring.

## Decision Log

1. Do not submit DockQ against unrefined models because the requested endpoint
   is refined-model DockQ.
2. Do not submit PyRosetta on KUACC after the import probe failed; the adapter
   is fail-closed and must not silently fall back to external Rosetta.
3. Keep downstream preparation separate from the original run and retain
   per-candidate paths, native mapping, chain counts, and input hashes.
4. Treat the 18 chain-count-mismatch cases as review-required, not silently
   repair them by changing native mappings.

## Outcomes & Retrospective

No downstream jobs were submitted in this phase. Preparation scripts are
implemented and locally syntax-checked; runtime/backend validation remains
the gating item.

## Context and Orientation

- Source run: `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919`
- Compact candidates: `compact_results_20260927/output/candidate_generated.csv`
- Candidate rows: 88,275
- Expected source pipelines: TMalign and MultiProt only for refinement
- PyRosetta bundle: `.../pyrosetta_bundle`
- KUACC PyRosetta wrapper: `.../kuacc_submission_20260919_full_bm55/pyrosetta_run_python.sh`
- DockQ environment: `/kuacc/users/rshadi25/.conda/envs/prism_dockq_20260919/bin/python`

## Plan of Work

1. Keep the existing compact case/candidate evidence as the source of truth.
2. Generate a downstream manifest from candidate rows and per-case attempt
   roots, checking transformed pair existence and native mapping.
3. Run a compatible PyRosetta preflight on the selected backend.
4. Submit resumable batched refinement only after preflight passes.
5. Validate refined output counts, hashes, chain composition, timings, and
   no-drop accounting.
6. Submit resumable DockQ/iRMSD batches against refined outputs only.
7. Aggregate final tables and preserve failures/missing/skipped/not-run states.

## Concrete Steps

The prepared manifest builder uses one JSONL batch per 500 candidates, checks
native and transformed inputs, records chain counts, and writes invalid-input
rows separately. The prepared DockQ worker assembles receptor/ligand parts
with native chain IDs, scores DockQ and backbone iRMSD, and appends a
per-candidate checkpoint before producing each batch CSV. The PyRosetta probe
records Python, glibc, bundle, version, and import error and exits nonzero on
failure.

## Validation and Acceptance

Before submission, require:

- PyRosetta import success on the actual compute backend;
- no unresolved chain-count mismatches for submitted candidates;
- exact source manifest and candidate-table hashes recorded;
- all refined output rows linked to one input candidate;
- DockQ only starts after refinement terminal accounting is complete;
- per-stage timings, return codes, hashes, and failure categories are present.

## Idempotence and Recovery

All workers use new downstream roots, per-candidate/batch checkpoints, atomic
summary writes, and skip completed records. A disappeared scheduler job is
`UNKNOWN` until durable status and output evidence are present. No source
alignment or transformed files are deleted by this plan.

## Artifacts and Notes

Prepared scripts:

- `scripts/prism-prescript/kuacc/prepare_downstream_manifest.py`
- `scripts/prism-prescript/kuacc/score_dockq_batch.py`
- `scripts/prism-prescript/kuacc/score_dockq_batch.sbatch`
- `scripts/prism-prescript/kuacc/probe_pyrosetta_runtime.py`

## Interfaces and Dependencies

Dependencies are the validated candidate table, corrected batch manifests,
per-case task ledgers, transformed PDB pairs, native PDBs, PRISM evaluation
modules, DockQ environment, and a compatible PyRosetta runtime. External
Rosetta is not an implicit substitute for PyRosetta.
