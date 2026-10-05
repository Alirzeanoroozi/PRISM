# Evaluate ranking integration and pipeline effects

This ExecPlan records the reproducible evaluation of the optional candidate-ranking path in PRISM-prescript. It is a living document and should be updated as validation proceeds.

## Purpose / Big Picture

Determine whether ranking is implemented in the current pipeline, whether baseline and ranked arms complete successfully, and what measurable effect ranking has on refinement load, runtime, and native-like accuracy.

## Progress

- [x] Read project guidance and project memory.
- [x] Inspect the CLI, audit writer, selector, scorer, and benchmark evaluators.
- [x] Complete focused regression and pipeline health checks.
- [x] Record final evidence and limitations.

## Surprises & Discoveries

- The retained paired Slurm smoke has identical baseline/ranked inputs and audit records, with 2 transformed candidates in both arms; ranking selects 1 for refinement when top-k is 1.
- The paired smoke demonstrates refinement-load reduction and successful completion, but it does not contain a controlled unranked-vs-ranked native score comparison.
- The retained 10-row labeled batch shows 3 rows with at least one native-like candidate, while the deterministic baseline selects a native-like top-1 for 1 row.
- The current-tree rerun (job 1400311) confirmed the same ranking load reduction, but the ranked arm took 62.7s in refinement versus 56.3s baseline; Rosetta runtime variation dominates this one-pair comparison.
- Re-ranking the retained labeled batch with the current scorer changed the pilot result to 3/10 native-like top-1 groups; the prior CSV was generated with an older, unversioned score.

## Decision Log

- Decision: Treat the paired smoke as implementation/load evidence and the labeled batch as ranking-quality evidence; do not combine them into a causal accuracy claim.
  Rationale: They use different artifacts and the paired smoke has no native evaluation.
  Date/Author: 2026-07-28 / Codex
- Decision: Do not change production defaults during this evaluation.
  Rationale: Ranking remains opt-in and the available quality evidence is pilot-scale.
  Date/Author: 2026-07-28 / Codex

## Outcomes & Retrospective

Ranking is truly implemented and operational in the current pipeline. The
current evidence supports a 2-to-1 refinement-load reduction in the paired
smoke, but not a stable wall-clock speedup or a general accuracy claim. The
current-code pilot retained all three oracle-positive groups at top-1, but it
is GTalign-derived, small, and missing the coverage fields present in the live
audit path; a matched current-TMalign labeled benchmark remains necessary.

## Context and Orientation

The current entry point is `prism.py`. It performs transformation filtering, optionally invokes `src/candidate_selector.py` when `--rank` is enabled, then invokes the selected refiner. Candidate features are appended by `src/candidate_audit.py`; the deterministic baseline is in `src/candidate_ranker.py`. Paired smoke evidence is under `tmp/agent/20260726-isolated-ranked-smoke/runs/1392725/`; labeled ranking evidence is under `tmp/agent/20260726-bm55-canonical-scoring/ranking-batch-0001/`.

## Plan of Work

First verify implementation and audit wiring from source. Then verify the paired baseline/ranked run using stage-status records, return codes, candidate counts, and refinement outputs. Evaluate quality from independent labeled candidates using per-row top-1 versus oracle native-like outcomes. Finally run focused tests and a current-tree paired smoke through Slurm, then update project memory with only durable findings.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, inspect `prism.py`, `src/candidate_selector.py`, `src/candidate_ranker.py`, `src/candidate_audit.py`, and the retained run manifests.
2. From the same directory, inspect exact `results.tsv`, stage-status JSONL, audit JSONL, ranking CSV, and ranking evaluation JSON.
3. From the same directory, run the focused ranking/regression test suite with `/home/rshadi25/.conda/envs/gtalign_env/bin/python` and execute the current-tree paired smoke through Slurm.
4. Update the project memory references with verified conclusions and unresolved quality/causality limits.

## Validation and Acceptance

Acceptance requires: source-level ranking call before refinement; paired baseline and ranked arms both return zero with matching input/audit evidence; ranked arm selects fewer candidates; ranking stage has completed status; focused tests pass; and accuracy is reported as top-1 native-like rate/regret versus an oracle, with denominator and limitations.

Validation result: the focused suite passed `38 passed, 1 skipped` in 21.69s. The current-tree Slurm smoke passed with job 1400311 (`0:0`), identical source/input/audit manifests, and 2 baseline versus 1 ranked refined structure.

## Idempotence and Recovery

This evaluation is read-only except for this plan and durable memory notes. Existing benchmark and smoke artifacts are preserved. Any new temporary evaluator output must be written under `/tmp` or `tmp/agent/` and removed after inspection.

## Artifacts and Notes

- Paired smoke results: `tmp/agent/20260726-isolated-ranked-smoke/runs/1392725/results.tsv`
- Paired stage records: `tmp/agent/20260726-isolated-ranked-smoke/runs/1392725/{baseline,ranked}/status/stages.jsonl`
- Labeled ranking evaluation: `tmp/agent/20260726-bm55-canonical-scoring/ranking-batch-0001/ranking_evaluation.json`
- Current-code ranking evaluation: `tmp/agent/20260728-ranking-pipeline-evaluation/ranking_evaluation_current.json`
- Current-tree paired smoke: `tmp/agent/20260728-ranking-pipeline-evaluation/runs/1400311/`

## Interfaces and Dependencies

The paired smoke used Slurm job 1392725 on `cosbi`, Python 3.11.13 from `gtalign_env`, and external Rosetta 2022.42. The focused tests use the same pipeline interpreter. DockQ labels in the retained batch were produced with the project’s canonical DockQ 2.1.3 scoring environment.
