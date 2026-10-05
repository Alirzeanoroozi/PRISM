# KUACC PRISM intermediate cleanup evidence

This is a compact execution record for the confirmed cleanup of upstream
alignment artifacts. The remote namespace is
`/scratch/users/rshadi25/valar-remote-runs/cleanup-manifests/20260928-prism-intermediate-dry-run/`.

## Scope

- Frozen plan: `alignment_cleanup_plan.json`
- Planned targets: 517 exact directories
- Target form: `processed/alignment_*` below the BM55 full-run source root
- Preserved: `processed/transformation`, refinement structures/energies,
  DockQ outputs, checkpoints, manifests, scripts, logs, native/template PDBs,
  and the active large-refinement root

## Execution

- User-confirmed deletion scope: upstream alignment directories only.
- Resumable implementation: `scripts/prism-prescript/kuacc/delete_alignment_dirs_fast.py`
- Authoritative result ledger: `deleted_alignment_directories.tsv`
- Next-round skip list: `alignment_cleanup_skiplist.tsv` (materialized by the
  updated script from only `deleted`/`already_absent` ledger records)
- Event ledger: `alignment_cleanup_events.jsonl`
- A duplicate delayed detached launcher was stopped by exact PID; the single
  directly monitored cleanup process was retained.

## Evidence at record time

- 4 target directories deleted and 1 target was already absent; 513 planned
  target directories remained at the last check.
- Representative transformed data remained present.
- The large refinement and GlobalDockQ-v4 roots remained present.
- KUACC refinement jobs remained queued/running and were not modified.

The deletion process was still running at the time this record was written;
final counts must be read from the remote ledger rather than inferred from
this snapshot.
