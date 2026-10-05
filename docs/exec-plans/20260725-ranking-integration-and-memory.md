# Repair Opt-in Ranking Integration and Consolidate Project Memory

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` up to date as work proceeds.

## Purpose / Big Picture

Make `prism.py --rank` a reproducible opt-in pre-refinement selection stage without changing the stable unranked pipeline. Replace the accumulated project memory with concise, evidence-based operational guidance.

## Progress

- [x] Read project-local guidance and all three memory files.
- [x] Trace ranking CLI, transformer audit creation, selector behavior, and focused tests.
- [x] Repair runtime audit routing and fail-open selection behavior.
- [x] Add focused tests and validate the ranking path.
- [x] Rewrite durable memory and record final evidence.

## Surprises & Discoveries

- Observation: `src.transformation` captures `PRISM_CANDIDATE_AUDIT_PATH` at module import, before `prism.py` parses CLI options.
  Evidence: `src/transformation.py` module constant `AUDIT_PATH`; `prism.py` imports `transformer` before `main()`.
- Observation: a missing audit causes a safe skip, but nonempty unmatched/stale audit records can leave a pair group with no selected candidates.
  Evidence: `src/candidate_selector.py:select_top_candidates()` drops a group when `rankable_rows` is empty.

## Decision Log

- Decision: ranking remains disabled unless `--rank`/`PRISM_RANK` is set.
  Rationale: the latest grouped pilot did not establish an improvement over deterministic selection.
  Date/Author: 2026-07-25 / Codex
- Decision: a ranked run without an explicit audit path receives a fresh run-scoped JSONL audit under `processed/candidate_audit/`.
  Rationale: this prevents stale append-only records from influencing a new run.
  Date/Author: 2026-07-25 / Codex

## Outcomes & Retrospective

Completed 2026-07-25. `prism.py --rank` now receives a fresh audit unless the
caller deliberately supplies one, and passes it directly into the transformer.
The transformer resolves the destination at write time rather than module
import. The selector preserves unmatched transformed candidates, preventing a
stale/partial audit from suppressing refinement. The stable no-rank path keeps
its existing transformer call and defaults.

Focused validation passed: 21 tests covering candidate selection, audit
serialization, rank scoring/data/table behavior, and CLI compatibility.

2026-07-25 follow-up: a historical reconstruction was added at
`docs/memory-archive/20260725-pre-consolidation-reconstruction.md`. It is
explicitly non-authoritative and did not change code, result artifacts, or
active memory.

2026-07-25 archive correction: exact pre-consolidation copies were found in
the separate `/home/rshadi25/GitHub/PRISM-prescript` checkout and copied to
`docs/memory-archive/20260725-exact-pre-consolidation/`. SHA-256 checksums are
recorded in that directory's README.

## Context and Orientation

`prism.py` runs download, surface extraction, alignment, transformation, optional ranking, and one selected refiner. `src/transformation.py` emits eligible transformed PDB pairs and can write audit records. `src/candidate_selector.py` maps those pairs back to audit records and retains top-K candidates per receptor/ligand pair.

The stable entry points remain `benchmark/scripts/run_prism_pipeline_smoke.sh` and the documented direct `prism.py` commands. This work changes only the opt-in ranking branch and project-local Markdown memory.

## Plan of Work

Make audit selection explicit at transformation-call time, while retaining the old environment variable as a runtime fallback. Add a CLI audit-path option and a unique default for ranked runs. Preserve unmatched transformed candidates rather than silently deleting them. Test the two regressions, then document verified results.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, patch `prism.py`, `src/transformation.py`, and `src/candidate_selector.py`.
2. Add regression tests in `tests/test_candidate_selector.py`, `tests/test_transformation_audit.py`, and `tests/test_prism_cli.py`.
3. Ran `/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q` against focused selection, audit, ranker/data, CLI, and table tests: 21 passed.
4. Rewrote the three `project-memory/references/*.md` files with only current verified findings and retained blockers.

## Validation and Acceptance

Acceptance requires:

- `--rank` allocates a run-scoped audit path when one is not supplied.
- the transformation stage writes to a path passed after import time.
- an unmatched audit row cannot reduce a nonempty transformed set to zero.
- ranking and existing audit/ranker tests pass.
- normal no-rank invocation does not request ranking or change its default refiner/aligner behavior.

## Idempotence and Recovery

No existing audit file is overwritten. Automatically generated ranked-run paths include timestamp and PID. Explicit audit paths retain append-only behavior. Reverting this change is limited to the touched source, test, documentation, and memory files.

## Artifacts and Notes

- Stable operational guide: `docs/STABLE_PIPELINE.md`.
- Ranking design/evidence: `docs/ML_TRAINING.md` and `docs/biological-ranking-pilot-20260711.md`.
- Memory authority: `.agents/skills/project-memory/references/`.

## Interfaces and Dependencies

- `prism.py` CLI: `--rank`, `--top-k`, `--rank-min-score`; new audit-path option remains opt-in.
- `src.transformation.transformer(templates, alignment_dir=..., audit_path=...)` must stay compatible with callers that omit `audit_path`.
- `src.candidate_selector.select_top_candidates()` returns the same `(left_pdb, right_pdb)` tuple interface consumed by all refiners.
