# Add Multi-Chain TM-align Support and Re-run the Full PRISM Comparison

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` current while the benchmark proceeds.

## Purpose / Big Picture

Enable current PRISM to consume the same chain-qualified and multi-chain target identifiers as legacy PRISM, while retaining the existing single-chain path. Then create one normalized manifest from all rigid, medium, and difficult benchmark rows, execute both `TMalign + Rosetta` and `MultiProt + FiberDock` on identical ten-pair batches, re-run DockQ/iRMSD analysis, and publish a pairwise discrepancy report with stage-specific failures and evidence-backed explanations.

## Progress

- [x] Read project memory and identify the current-vs-legacy comparison scope.
- [x] Add and pass target normalization/materialization tests, including multi-chain and legacy residue-range tokens.
- [x] Preserve all source chains through the current Rosetta combine/contact handoff.
- [x] Add shared full-benchmark manifest and ten-pair batch preparation utility.
- [ ] Add reproducible Slurm launchers for both pipelines and validate their batch inputs.
- [ ] Submit and collect every full-benchmark batch.
- [ ] Re-run benchmark analysis and produce pairwise comparison/failure report.
- [ ] Record environment limitation and final verification in project memory.

## Surprises & Discoveries

- The current downloader rejected every target that was not exactly five characters, while legacy preprocessing already supported a four-character PDB ID followed by multiple chain IDs.
- Benchmark identifiers include underscores and residue-range annotations such as `1IK0_A(10)` and empty-chain tokens such as `3LZT_`; normalization must remove these presentation artifacts before chain materialization.
- Current Rosetta preprocessing inferred one source chain from a filename and would discard additional chains in a multi-chain partner. The corrected path assigns a unique chain group to each partner and preserves every requested chain.
- The dedicated `/scratch/tmp/prism-current-test-py311` Python 3.11 environment reaches interpreter initialization with `init_fs_encoding`/filesystem-codec failure when importing modules; focused checks therefore run in `gtalign_env`.

## Decision Log

- Decision: Normalize target IDs to lowercase four-character PDB ID plus sorted unique alphanumeric chain IDs, matching legacy behavior.
  Rationale: This makes `1ABC_BA`, `1ABC_AB`, and `1abcAB` refer to the same multi-chain input while preserving chain identity.
- Decision: For empty chain suffixes, retain the full downloaded structure rather than inventing a chain filter.
  Rationale: Benchmark rows contain tokens such as `3LZT_`, and legacy PRISM treats them as whole-structure inputs.
- Decision: Keep benchmark rows distinct in the shared manifest even when normalized pairs repeat across sets.
  Rationale: Set membership and source-row provenance are needed for faithful per-benchmark analysis.
- Decision: Report FiberDock and Rosetta energies separately; compare DockQ, iRMSD, stage status, and counts instead of numeric energy magnitudes.
  Rationale: The scoring functions and scales are not calibrated.

## Outcomes & Retrospective

To be completed after batch collection and analysis.

## Context and Orientation

Current input preparation is in `src/pdb_download.py`; surface ASA filtering is in `src/naccess_utils.py` and `src/surface_extract.py`; transforms are written by `src/transformation.py`; Rosetta handoff is in `src/rosetta_refinement.py`. Focused tests are in `benchmark/scripts/test_prism_pipeline_helpers.py`, `tests/test_transformation_thresholds.py`, and `tests/test_comparison_batches.py`. Full benchmark rows are `benchmark/data/T_Rigid.csv`, `benchmark/data/T_medium.csv`, and `benchmark/data/T_difficult.csv`.

## Plan of Work

Use the normalized manifest as the only source of pair inputs. Each batch contains at most ten rows and is copied unchanged into current `inputs.csv` and legacy `pair_list` forms. Current and legacy workspaces remain isolated; only the pair manifest is shared. Jobs run on CPU Slurm nodes, with downloads/setup performed before submission or reused from staged PDB/template trees. Analysis consumes both output trees and writes explicit success, missing-output, scoring-error, and execution-failure rows.

## Concrete Steps

1. From the repository root, run the focused `gtalign_env` tests and the batch-preparation utility with `--batch-size 10`.
2. Validate every generated batch has 1–10 rows and that current/legacy input forms contain identical normalized pairs.
3. Stage current and legacy code/tool assets in isolated run roots under `tmp/agent/<run-id>/`; do not overwrite raw benchmark data or prior results.
4. Verify live hostname, Slurm job context, queue, account, and QOS; submit one job per batch per pipeline with explicit logs.
5. Collect output manifests, run existing DockQ/iRMSD analysis, then join the two pipelines by `pair_id` and normalized PDB pair.
6. Classify divergences by last successful stage: input/download, preprocessing/chain selection, alignment, transform filtering, refinement, output discovery, or metric scoring.

## Validation and Acceptance

The implementation must pass focused tests in `gtalign_env`, retain the existing single-chain helper behavior, preserve every requested chain in a multi-chain combined PDB, and produce a manifest covering all 257 benchmark rows in batches of ten. Full acceptance additionally requires collected outputs or explicit failures for every submitted job, repeated benchmark metrics, and a pairwise report separating observations from likely causes.

## Idempotence and Recovery

Batch preparation is deterministic for a given data directory and output root. Never delete previous runs; rerun into a new timestamped `tmp/agent` root. A failed batch can be resubmitted independently, and analysis must retain failed rows rather than silently dropping them.

## Artifacts and Notes

The shared manifest/batch generator is `benchmark/scripts/prepare_full_comparison_batches.py`. Planned run artifacts belong under `tmp/agent/<run-id>/comparison_batches/`, with Slurm logs below each pipeline’s run root and final reports under `benchmark/comparison_reports/` or an explicitly named run root.

## Interfaces and Dependencies

The current pipeline uses `gtalign_env` for Python checks and the Rosetta 2022.42 module/binaries for refinement. The legacy pipeline requires the available Python 2.7 interpreter, MultiProt, FiberDock, NACCESS, and its existing legacy configuration. DockQ/iRMSD analysis uses the established benchmark scripts and must retain chain-group context for multi-chain rows.
