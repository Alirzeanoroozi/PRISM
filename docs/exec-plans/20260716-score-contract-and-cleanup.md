# Align Benchmark Scoring with PRISM-main and Restore Scoreable Benchmark Output

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` up to date as work proceeds.

## Purpose / Big Picture

Stop the currently invalid Benchmark 5.5 scoring path, establish whether `PRISM-prescript` invokes the same DockQ and iRMSD logic as `/scratch/rshadi25/GitHub/PRISM-main/benchmark`, and prove a small current-pipeline result is scoreable before launching a replacement benchmark. Consolidate disposable root-level logs into a timestamped `tmp/agent` quarantine; preserve source, raw benchmark inputs, validated artifacts, and active-run evidence.

## Progress

- [x] Record active Benchmark 5.5 generation/scoring state and preserve its evidence.
- [x] Compare the exact PRISM-main and PRISM-prescript scoring code, inputs, command-line flags, chain mapping, and native selection.
- [x] Identify the first causal mismatch using a small existing generated model and a canonical positive control.
- [x] Implement the smallest scoring/output-contract repair and focused static coverage.
- [x] Run a controlled one-model scoreability smoke that produces valid DockQ and iRMSD; retain the normal compute smoke as optional confirmation because the scheduler was delayed.
- [x] Submit a replacement benchmark only after the smoke has passed; preserve the current run as invalid-generation evidence.
- [x] Produce a dry-run cleanup manifest, quarantine only confirmed log artifacts, and verify the repository remains reproducible.

## Surprises & Discoveries

- Observation: the active 26-batch current run has generated many Rosetta models but interim strict `--no_align` scoring returned no DockQ or iRMSD values because residue correspondence validation failed.
  Evidence: `tmp/agent/20260715-benchmark55-full/scoring/progress_scores2.csv`.
- Observation: direct external-Rosetta output and copied `structures/` output layouts coexist; model discovery was already expanded to cover both.
  Evidence: `benchmark/scripts/score_comparison_models.py` and current run outputs.
- Observation: `irmsd.py`, `irmsd_backbone.py`, and both `rosetta_output` analyzers are byte-identical to PRISM-main.  The current collector instead adds a strict `--no_align` residue-numbering gate and whole-partner mapping; it is not the PRISM-main benchmark contract.
  Evidence: file hashes and `benchmark/scripts/score_comparison_models.py` versus `PRISM-main/benchmark/scripts/rosetta_output/analyze_prism_rigid_results.py`.
- Observation: the unmodified PRISM-main scorer computed iRMSD 29.599 A for current generated model `2gk2AB_1fgnHL_1tfhA...rosetta_0001.pdb`; local DockQ was blocked only when its multiprocessing manager attempted to create a sandbox-forbidden socket.
  Evidence: `tmp/agent/20260716-score-contract-cleanup/main-score-1ahw.csv` and its captured local DockQ error.
- Observation: the root cleanup manifest contains 1,639 untracked `*.out`, `*.err`, or `*.log` artifacts totalling 885,471 bytes, all classified `review-or-quarantine`.
  Evidence: `tmp/agent/20260716-score-contract-cleanup/cleanup_manifest.tsv`.
- Observation: the PRISM-main analyzer resolves its sibling metric scripts from `rosetta_output/`, where they do not exist. A direct invocation therefore returns metric-subprocess exit code 2 even though the analyzer source is byte-identical to the local copy.
  Evidence: `tmp/agent/20260716-score-contract-cleanup/main-contract-local-serial-output/per_prediction.csv`.
- Observation: a byte-identical disposable analyzer copy with adjacent symlinks to the original PRISM-main metric scripts produces numeric DockQ 0.8912906078 and iRMSD 29.599 A for the staged generated 1AHW example, with no scoring error.
  Evidence: `tmp/agent/20260716-score-contract-cleanup/main-contract-local-exact-output/per_prediction.csv` and `main-contract-local-runner/` SHA256 comparison.

## Decision Log

- Decision: do not treat the active generation run as a valid benchmark result until the scoring contract is reconciled.
  Rationale: null metrics from a strict evaluator are not quality measurements.
  Date/Author: 2026-07-16 / Codex.
- Decision: use the PRISM-main benchmark implementation as the reference requested by the user; do not substitute an ad hoc evaluator.
  Rationale: direct script and invocation comparison is required to diagnose measurement drift.
  Date/Author: 2026-07-16 / Codex.
- Decision: score the replacement benchmark by invoking the unmodified PRISM-main analyzer directly.  Stage only symlinks with the filename grammar it requires, preserving current PDB coordinates and chain IDs.
  Rationale: this is the smallest repair that enforces the requested benchmark contract without rewriting a validated metric implementation.
  Date/Author: 2026-07-16 / Codex.

## Outcomes & Retrospective

Pending score-contract audit and controlled smoke.

The score-contract audit is complete. Superseded generation/scoring jobs `1360392` and `1360454` were cancelled. Initial replacement jobs `1361649`/`1361656` failed immediately because relative input roots were resolved from Slurm's spool directory; this exposed and fixed a launcher path bug. Corrected array `1361684` and dependent canonical scoring job `1361687` use a fresh run root with absolute paths. The current full benchmark has not yet completed, so no aggregate quality report exists. Root-level logs were quarantined rather than deleted at `tmp/agent/20260716-score-contract-cleanup/root_logs_quarantine/`.

Latest check: batches 2–5 are completed successfully, batch 1 is an explicit `completed_no_predictions` case, and batches 6–7 are running. The partial run has 416 transformation PDBs and 67 Rosetta PDBs; scoring has not started and no aggregate DockQ/iRMSD report is available.

Scoring validity correction: the earlier PRISM-main analyzer contract was rejected for multichain scoring after direct audit. The known-invalid dependent job was cancelled. `benchmark/scripts/score_bijective_benchmark_models.py` passed the 1AHW smoke with explicit `LHA:ABC` assignment, cross-interface DockQ selection, and grouped iRMSD; corrected dependent job `1363282` is queued after generation `1361684`.

## Context and Orientation

The working repository is `/scratch/rshadi25/GitHub/PRISM-prescript`. The reference scoring directory is `/scratch/rshadi25/GitHub/PRISM-main/benchmark`. Current generation is Slurm array `1360392`; dependent scorer `1360454` must not be treated as authoritative until the contract audit passes. Important local scripts are `benchmark/scripts/score_single_prism_pair.py`, `benchmark/scripts/score_comparison_models.py`, `benchmark/scripts/irmsd.py`, and the `benchmark/prism_processed` native-complex tree.

## Plan of Work

First compare file hashes, imports, interfaces, command construction, and native/chain mapping across the two repositories. Then reproduce the mismatch on one model and run an independent positive-control self-score through both implementations. Repair only the demonstrated divergence, add a test, and submit a small Slurm smoke through the production launcher/scorer. Only a smoke that yields valid non-null DockQ and iRMSD permits a replacement benchmark submission. Cleanup is based on a reviewed manifest and uses quarantine rather than deletion for pre-existing root-level logs.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, compare `benchmark/scripts` scoring files to `/scratch/rshadi25/GitHub/PRISM-main/benchmark` and record exact differences.
2. Inspect the current scorer's model/native chain declarations and validate a positive-control self-score using each implementation.
3. Run a one-model current-output score under both contracts and capture model, native, chain mapping, return code, DockQ, iRMSD, and raw errors.
4. Patch only the identified divergence; run focused tests and a Slurm smoke on the smallest representative pair/template panel.
5. Build a `tmp/agent/20260716-score-contract-cleanup/cleanup_manifest.tsv`, review classifications, and move only confirmed root-level log files to a dated quarantine directory.
6. After successful smoke metrics, submit a new benchmark root with the fixed scorer and explicit provenance; do not overwrite `tmp/agent/20260715-benchmark55-full`.

## Validation and Acceptance

- The two score implementations either have identical relevant code/behavior or every justified difference is documented.
- A positive control returns valid DockQ and iRMSD through the adopted contract.
- A current-pipeline small smoke returns valid non-null DockQ and iRMSD with explicit model/native chain declarations.
- Focused tests, Python compilation, and shell syntax pass.
- Replacement submission records exact environment, commands, chain mapping, and score outputs.
- Cleanup manifest shows no tracked, source, raw-data, validated-result, memory, or active-run path was deleted.

## Idempotence and Recovery

All new outputs use `tmp/agent/20260716-score-contract-cleanup/` or a new timestamped run root. The active `20260715` benchmark root is read-only evidence. Quarantine moves preserve relative paths and can be restored by moving files back from the quarantine root. Any failed smoke is retained with its log and exit record; no benchmark is resubmitted unless the scoreability acceptance check passes.

## Artifacts and Notes

- Active generation: `tmp/agent/20260715-benchmark55-full/runs/current/`.
- Interim scores: `tmp/agent/20260715-benchmark55-full/scoring/progress_scores2.csv`.
- Cleanup manifest/quarantine: `tmp/agent/20260716-score-contract-cleanup/`.
- Reference implementation: `/scratch/rshadi25/GitHub/PRISM-main/benchmark`.

## Interfaces and Dependencies

- Pipeline Python: `/home/rshadi25/.conda/envs/gtalign_env/bin/python`.
- DockQ Python: `/scratch/tmp/prism-dockq-env/bin/python`.
- External Rosetta: module `rosetta/2022.42` inside Slurm jobs.
- Replacement score entry point: `benchmark/scripts/score_benchmark55_main_contract.sbatch`, which calls `/scratch/rshadi25/GitHub/PRISM-main/benchmark/scripts/rosetta_output/analyze_prism_rigid_results.py` directly after `benchmark/scripts/stage_current_models_for_main_benchmark.py` creates symlinks.
- Slurm execution uses `cosbi` account/partition currently; verify live scheduler state before any replacement submission.
