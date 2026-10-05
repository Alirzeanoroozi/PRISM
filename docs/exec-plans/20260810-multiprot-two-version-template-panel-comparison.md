# Two-version MultiProt template-panel comparison

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` aligned with the actual run state.

## Purpose / Big Picture

Run a controlled comparison of the maintained current MultiProt adapter and the legacy `working_version/Multiprot-new/prism-fiberdock-cli` MultiProt path on two fixed docking pairs, using the largest pre-existing shared native-interface panel (770 IDs) from the two differently formatted template trees. Partition the panel into ten deterministic subsets, preserve per-subset outputs, and compare raw alignment results without conflating current and legacy contracts.

## Progress

- Audited both template lists and native interface formats.
- Selected 770 exact shared template IDs in current-list order; every selected ID has both current `*_int.pdb` sides and legacy `*.int` sides.
- Partitioned the frozen panel into 10 deterministic batches of 77 IDs.
- Frozen pairs: `1cew`–`2ghuD` and `2uwjG`–`2uwjE`, using retained surface-extraction inputs for both implementations.
- Scope is the MultiProt alignment stage; downstream filtering/refinement is intentionally excluded until native alignment outputs are compared.

- [ ] Recover current project, template, input, and runtime state.
- [x] Determine whether the two template lists are fundamentally different and freeze a manifest.
- [x] Prepare ten-way partitions and isolated current/legacy run workspaces.
- [x] Validate runners on representative bounded parser cases before/alongside the full shared panel.
- [x] Run the paired matrix through Slurm and collect explicit statuses.
- [x] Compare outputs and independently validate the evidence.

## Surprises & Discoveries

- The full lists are not fundamentally disjoint: the current list has 19,062 IDs and is a subset of the legacy list's 21,072 IDs; the legacy list adds 2,010 IDs and orders the shared IDs differently.
- The literal first-1,000 panels were not equally runnable with retained native assets, so the controlled panel uses 770 exact shared IDs with both current `*_int.pdb` and legacy `*.int` sides.
- Native interface files are mostly structurally equivalent by CA residue keys (1,524/1,540 sides equal; median common-coordinate RMSD 0.0 A), but 16 sides differ in keys and 27 have nonzero coordinate RMSD.
- Both executables have the same SHA256. On identical query/template/chain keys, current accepted 405/6,160 records while legacy accepted 5,968/6,160; current-only success was zero.
- A representative legacy-only case produced identical MultiProt output and 29 parsed matches in both native formats, while current's additional Kabsch reconstruction returned `inf`; this demonstrates an important adapter-contract difference but does not explain every unavailable record.

## Decision Log

- Decision: Keep current-pipeline and legacy MultiProt outputs in separate run roots and compare only after normalizing explicit metadata.
  Rationale: Project memory requires current/legacy FiberDock/MultiProt evidence to remain provenance-separated and warns against historical-equivalence claims from unmatched inputs or evaluators.
  Date/Author: 2026-08-10 / Codex
- Decision: Use one reproducible Slurm allocation with internal parallelism only if the live resource check supports it; otherwise use a throttled array with isolated task directories.
  Rationale: Ten template subsets are independent, but shared working directories and legacy binaries must not overlap.
  Date/Author: 2026-08-10 / Codex
- Decision: Compare the largest pre-existing shared native-interface panel (770 IDs), not a nominal 1,000-entry panel with silently missing legacy assets.
  Rationale: A paired benchmark must preserve identical template IDs and runnable native assets; missing legacy interfaces would confound list choice with asset availability.
  Date/Author: 2026-08-10 / Codex

## Outcomes & Retrospective

Slurm jobs 1492939 (main panel), 1492941 (legacy batch-8 retry), and parser diagnostics 1492944/1492945 completed with exit code 0. The evidence package contains 6,160 records per implementation, ten 77-template batches, native-format asset comparisons, status/match ledgers, and parser diagnostics. Focused tests passed 9/9. Terra independently returned `PASS WITH CAVEATS`: counts are reliable, but this is an end-to-end adapter comparison, not proof of raw-MultiProt equivalence or complete causal attribution. Future runs should record granular failure reasons, source commits, query hashes, and exact environments, then perform same-invocation orientation controls.

## Context and Orientation

The maintained pipeline is rooted at `prism.py` and `src/`; the legacy comparison root is `working_version/Multiprot-new/prism-fiberdock-cli`. The current project memory identifies the current MultiProt adapter as `src/alignment_multiprot.py` and the legacy helper stack as a compatibility/reference arm. Current and legacy results must not be treated as causal or historically equivalent without matched inputs, staging, and evaluator contracts.

The requested comparison is two pipeline versions using their native template-list/asset formats and two pairs. The shared panel contributes 770 identical template IDs, partitioned into ten 77-template subsets. The pair IDs, template-selection rule, chain roles, executable paths, environment, and output roots are frozen in a manifest before submission.

## Plan of Work

First inventory and compare the two template lists by canonical template/chain IDs, hashes, counts, duplicates, and ordering. Then inspect both MultiProt runners and their output contracts, select two valid chain-qualified pairs already supported by project fixtures, and prepare a dry-run manifest. After a bounded single-subset preflight, submit the independent subset jobs through Slurm with no shared mutable working directories. Aggregate raw alignment/status records before comparing downstream candidate/refinement results; report observations, inferences, and unresolved differences separately.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, verify `hostname`, `SLURM_JOB_ID`, `squeue -u "$USER"`, live partition/QOS availability, repository status, and canonical memory.
2. Locate both template lists, both MultiProt entry points, existing pair inputs, and any established runner scripts. Read only the relevant parsers, path constructors, output writers, and legacy instructions.
3. Create a new run root under `tmp/agent/20260810-multiprot-two-version-template-panel-comparison/` containing `manifest.tsv`, `template_lists/`, `partitions/`, `runs/`, `logs/`, and `aggregate/`. Do not alter raw or existing benchmark outputs.
4. Generate the deterministic 770-template shared selection and ten 77-template partitions, recording source paths, selection order, hashes, overlap, and exclusions.
5. Run the smallest safe preflight for each version/list/pair combination on one partition; stop if the legacy runtime, input contract, or output parser is not valid.
6. Submit only after preflight validation. Use explicit account/partition/QOS and either one internally parallelized CPU job or a throttled array, with one isolated output root per version/list/pair/subset and one status record per task.
7. Aggregate completed, failed, timed-out, and missing subsets without silently dropping any task. Compare raw alignments first, then transformed candidates and downstream artifacts only where the relevant stage completed.
8. Run file-scoped tests and validation scripts; inspect representative outputs and hashes. Do not claim biological quality or historical equivalence from this small panel.

## Validation and Acceptance

Acceptance requires: both lists are characterized quantitatively; exactly 1,000 selected templates and ten 100-template partitions are recorded per list; all 8 matrix conditions have explicit terminal status; no output root is shared across versions or subsets; raw alignment counts and failures are aggregated with status fields; representative output schemas parse; and the final comparison distinguishes observation from inference. A scheduler cancellation, missing legacy runtime, or incomplete subset matrix remains an explicit limitation rather than a successful run.

## Idempotence and Recovery

All generated artifacts use the unique run root and subset-specific paths. Reruns must use a new run ID or explicit missing-task retry directory. Never delete or overwrite raw lists, source templates, prior run evidence, or validated outputs. If a task fails, preserve its command, environment, stderr/stdout, exit code, and status JSON; retry only that task after diagnosing the failure.

## Artifacts and Notes

Planned artifacts: `manifest.tsv`, `template_lists/*.tsv`, `partitions/<list>/part_*.tsv`, per-condition `status.json`, Slurm stdout/stderr, raw alignment outputs, and `aggregate/*.tsv`/`*.json`. The final report will link these paths and record the Git state, interpreter, tool versions, Slurm job IDs, resources, and known limitations.

## Interfaces and Dependencies

Dependencies include Python 3.11 in `/home/rshadi25/.conda/envs/gtalign_env`, the current MultiProt adapter and executable, the legacy MultiProt binary/helper files, Biopython/NumPy, Slurm, and the two template-list formats. Exact executable paths and resource assumptions remain unknown until inspection and must be recorded before submission.
