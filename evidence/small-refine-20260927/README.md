# PRISM-prescript bounded refinement comparison — 2026-09-27

This evidence package records a bounded, high-confidence probe of retained
TMalign and MultiProt transformed models. It is not a complete BM55 benchmark
and does not support promotion or a final pipeline ranking.

## Scope and provenance

- 40 selected candidates: 20 TMalign and 20 MultiProt.
- Selection rule: generated rows with both transformed PDBs and the native PDB
  present; descending minimum of the two side TM-scores, then mean TM-score and
  total match count.
- Every candidate was passed through input chain normalization, FiberDock,
  external Rosetta, and DockQ/iRMSD where a refined model existed.
- No alignment, ranking, PRODIGY, or PyRosetta computation was introduced.
- Original KUACC run root:
  `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-small-refine-20260927`
- Source transformed/native artifacts:
  `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919`
- Full comparison array: KUACC job `3124119`, `1-19,21-39%16`.
- Smoke/retry jobs: `3124109` (initial adapter), `3124112` (chain-normalized
  smoke), and `3124116` (resumable DockQ compatibility retry). The initial
  adapter failure is preserved remotely and is not included as a success claim.

## Terminal accounting

- Input normalization: 40/40 completed.
- FiberDock: 40/40 refined models.
- External Rosetta: 34/40 models; 6/40 explicit `no_model` outcomes.
- FiberDock DockQ: 40/40 scored.
- Rosetta DockQ: 34/40 scored; 6/40 `not_run_no_model`.
- Candidate rows/checkpoints: 40/40; no-drop validation passed.

DockQ values in `comparison_results.csv` are corrected bounded `GlobalDockQ`
values read from preserved raw JSON. The legacy deployed KUACC parser used
`best_dockq` (an interface sum) as `dockq`; those legacy values are retained in
the `*_dockq_legacy` and `*_dockq_sum` columns for audit and are not used for
quality conclusions.

## Compact artifacts

- `selected_candidates.csv`: frozen 40-row input manifest with TM-score
  selection metrics, native chain mappings, and input hashes.
- `selection_manifest.json`: selection provenance and counts.
- `comparison_results.csv`: one row per candidate with per-stage statuses,
  timings, FiberDock energy, Rosetta scores, corrected DockQ/iRMSD, mappings,
  raw JSON paths/hashes, and model paths.
- `comparison_summary.json`: no-drop validation, stage accounting, and summary
  statistics.

The refined PDBs, raw DockQ JSON, checkpoints, and event ledger remain in the
isolated KUACC run root above; they were not copied into the framework
repository.

## Validation

- Local Python compilation and `git diff --check` passed for the selector,
  worker, batch wrapper, and aggregator.
- `pytest -q --ignore=tests/test_valar_interactive_qwen.py` passed. The full
  suite's Qwen launcher test is environment-sensitive here and intermittently
  fails at socket creation with `PermissionError`/`ports occupied`; it is
  unrelated to the comparison scripts.
- Local chain-normalization/canonicalization smoke passed with explicit `AB:AB`
  mapping.
- KUACC worker compilation passed before submission.
- `sbatch --test-only` passed for the smoke and full arrays.
- The final aggregate reports 40 selected rows and 40 checkpoints with
  `no_drop=true`.
- Corrected GlobalDockQ range across scored rows was within `[0, 1]`; refined
  model and raw JSON paths validated on KUACC with zero missing paths.

## Per-pipeline summary

Values below are median (mean) seconds or scores over available rows; the
complete per-candidate values are in `comparison_results.csv`.

| Pipeline | FiberDock time | Rosetta time | FiberDock GlobalDockQ | Rosetta GlobalDockQ | FiberDock energy | Rosetta interaction |
|---|---:|---:|---:|---:|---:|---:|
| TMalign (20) | 22.07 (25.09) s | 77.73 (92.79) s | 0.668 (0.532), n=20 | 0.658 (0.589), n=19 | -6.70 (-1.05) | -20.75 (-19.74), n=19 |
| MultiProt (20) | 14.46 (19.80) s | 59.46 (66.04) s | 0.370 (0.366), n=20 | 0.495 (0.470), n=15 | 0.70 (4.55) | -15.18 (-17.10), n=15 |

DockQ scoring time was small relative to refinement: FiberDock-DockQ median
3.34 s (TMalign) / 2.80 s (MultiProt), and Rosetta-DockQ median 3.11 s /
2.54 s, respectively. These values include only the bounded probe and are not
throughput estimates for the full panel.

## Interpretation boundary

This set was deliberately enriched for high alignment TM-score, so its DockQ
distribution is descriptive and selection-biased. It tests that both
refinement paths and scoring contracts execute on representative retained
models; it does not establish that either refiner improves the full pipeline.
