# PRISM workload cancellation evidence

## Scope

The user authorized cancellation of PRISM workloads on both VALAR and KUACC
before a controlled TMalign comparison. Only scheduler entries whose job name
matched the PRISM workload names were targeted; unrelated jobs were not
cancelled.

## Snapshots

- VALAR snapshot:
  `/scratch/rshadi25/GitHub/PRISM-prescript/benchmark/prism_processed_results/cancellation_snapshots/20260928T105918Z/`
- KUACC snapshot:
  `/scratch/users/rshadi25/valar-remote-runs/cancellation-snapshots/20260928T110120Z/`

Each snapshot contains the pre-cancellation queue, job metadata, targeted IDs,
cancellation passes, and post-cancellation verification. Persisted result
roots were retained.

## Verification

- VALAR remaining PRISM jobs: 0
- KUACC remaining PRISM jobs: 0
- VALAR comparison/refinement result roots: present
- KUACC refinement checkpoints root: present
- KUACC refinement event ledger: present
- The separate alignment-cleanup controller and its exact child were also
  stopped to freeze the input artifacts before the new comparison.

## Scientific comparison note

The next TMalign run can use the same downstream geometry settings and frozen
query/template panel. MultiProt should not be forced through the TMalign
TM-score gate because its native alignment contract is match count plus
coverage; this distinction must remain explicit in the comparison report.
