# DockQ scoring-contract correction

This ExecPlan is a living document. Keep `Progress`, `Decision Log`,
`Validation`, and `Outcomes & Retrospective` current as the work advances.

## Purpose / bounded goal

Repair the ad-hoc PRISM DockQ path so that its reported quality metrics match
DockQ 2.1.3 semantics and the benchmark's receptor-ligand contract. Preserve
the existing raw JSON and incorrect CSV as diagnostic evidence; write any
corrected scores to a new output root.

This plan covers the PRISM-side transformed-model helper and the selected
`PRISM-prescript` evaluator code. It does not submit Slurm jobs, overwrite
benchmark artifacts, change stable pipeline defaults, or treat transformed
models as refined models.

## Framework lifecycle

`Project -> Plan -> bounded Goal -> Agent execution -> Evidence -> Validation
-> Advancement`

Current state is the post-Green validation checkpoint. Production changes
cover both repository parser copies, the transformed-model caller, standalone
benchmark scoring, explicit cross-interface scoring, and notebook Slurm CPU
propagation; benchmark artifacts have not been regenerated.

## Progress

- [x] Recover Valar framework instructions, selected-project guidance, active
  project state, DockQ evidence, and dirty-worktree boundaries.
- [x] Define the metric contract: use bounded `GlobalDockQ` for complete
  mappings; use requested receptor-ligand cross-interface values as the
  primary multichain benchmark metric; keep iRMSD in Angstroms.
- [x] Add isolated RED tests for the ad-hoc parser and explicit CPU control.
- [x] Run the RED tests and observe the expected failures.
- [x] Run the existing canonical scorer/runtime/failure tests: 16 passed.
- [x] Complete three independent read-only audits for parser semantics,
  canonical cross-interface behavior, and command/resource provenance.
- [x] Implement the parser/CPU Green cycle in both parser copies.
- [x] Verify transformed-model caller wiring and row-level DockQ provenance
  with a synthetic end-to-end test.
- [x] Review RED evidence and authorize the minimal Green implementation.
- [x] Implement parser, command, and initial output-contract corrections.
- [x] Run targeted validation, a representative raw JSON replay, and a
  complete retained-JSON range/identity audit.
- [x] Record that pre-existing repository-wide pytest collection conflicts make
  an unscoped full-suite run invalid as a gate for this change.
- [ ] Refresh Graphify after the code change; the attempted update was
  interrupted after traversing the large dirty/evidence tree without
  completing.
- [x] Synchronize the production PRISM runtime parser and add explicit
  requested receptor-ligand cross-interface scoring to the transformed-model
  helper.
- [x] Correct the PRISM standalone benchmark scorer: JSON-first parsing,
  bounded GlobalDockQ, valid-unscored status, no-align guard, collision-safe
  raw JSON, hashes, argv, and explicit CPU provenance.
- [x] Propagate explicit DockQ CPU settings through every notebook scoring and
  validation Slurm entry point.
- [x] Run final focused PRISM and PRISM-prescript tests and compile checks;
  the final focused totals are 50 PRISM tests and 57 prescript tests.
- [x] Complete an independent Terra validation over both absolute repository
  paths after the final no-align and low-level CPU-default corrections.
- [ ] Regenerate corrected scores and EDA artifacts in a new directory.
- [x] Update project memory after validation; keep this plan open for the
  separately authorized corrected-score and EDA regeneration.

## Confirmed correction status

| Recommendation | Evidence status | Action |
|---|---|---|
| Report `GlobalDockQ`, not `best_dockq` or an interface fallback | Confirmed defect: RED tests report `1.4` instead of `0.7` and promote a missing-global interface score | Fix parsers and add bounded-range/valid-unscored regressions |
| Do not promote first-interface components in multichain output | Confirmed defect: RED test receives first-interface values | Return component fields only when their scope is unambiguous |
| Score requested receptor-ligand interfaces separately | Existing canonical path validated: relevant tests pass | Reuse canonical scorer; do not broaden the ad-hoc path silently |
| Preserve explicit partial/failure states | Existing canonical failure tests pass | Filter complete EDA by score status; retain failures explicitly |
| Make CPU use explicit | Confirmed contract gap: `n_cpu` is not accepted or passed | Add an explicit parameter and pass `--n_cpu 1` for the one-CPU job |
| Distinguish transformed from refined models | Confirmed from run metadata and input discovery | Keep separate score labels and output roots; no model substitution |

## Parallel audit findings

Three disjoint read-only audit lanes were dispatched using
`superpowers:dispatching-parallel-agents` and then closed after their
evidence was returned.

- The aggregation audit independently verified the retained raw-record
  identity: `best_dockq` is the per-interface sum and `GlobalDockQ` is the
  bounded average. It also found 4,696 multi-interface records and confirmed
  that the duplicate parser in `benchmark/scripts/dockq.py` must be kept in
  sync with `src/eval/dockq.py`.
- The canonical-path audit confirmed that requested cross-interface metrics,
  complete bijective mappings, GlobalDockQ separation, and explicit failure
  handling are implemented and covered by focused tests. It found that the
  single-pair/comparison path does not enforce the canonical no-align gate.
- The command/provenance audit confirmed that default sequence alignment is
  justified for the inspected inputs, but CPU usage is implicit (`DockQ`
  defaults to eight CPUs while the job allocates one), and the ad-hoc CSV
  lacks exact argv, raw-JSON linkage, tool version, and runtime/resource
  provenance.

These findings refined the Green acceptance criteria and identified a
comparison-path follow-up that must not be silently mixed into the parser-only
fix. The approved Green implementation now also includes the no-align guard
and single-pair provenance fields; the comparison-path guard is covered by
`tests/test_model_output_integrity.py`.

The transformed-model caller now preserves complete-complex `GlobalDockQ` and
separately scores every requested receptor-ligand chain pair. Cross-component
rows include status, mapping, raw JSON path/hash, argv, and CPU count;
cardinality mismatches remain explicit cross-score failures.

The PRISM standalone benchmark helper now follows the same JSON-first
contract. A valid JSON document with no usable complete-complex
`GlobalDockQ` is recorded as `valid_unscored`; it is not counted in the
`scored` total. Notebook scoring scripts pass `--n-cpu` explicitly from
`DOCKQ_N_CPU`/`SLURM_CPUS_PER_TASK`, defaulting to one when no Slurm value is
available.

The canonical staged scorer's recovery path is intentionally different: when
DockQ's complete mapping crashes, it may score requested receptor-ligand
interfaces pairwise. Those rows retain `GlobalDockQ` as unavailable, use the
separate `requested_cross_interfaces_only` scope, and are labeled
`scored_cross_only`.

The low-level parser helpers now default to `--n_cpu 1` when callers omit a
CPU setting. Resource-aware wrappers may still pass the actual Slurm
allocation explicitly; this default prevents an unmodified caller from
silently requesting DockQ's larger internal parallelism.

## RED evidence

New tests:

`tests/test_dockq_parser_contract_regression.py`

Command:

```bash
python -m pytest -q tests/test_dockq_parser_contract_regression.py
```

Observed result: `3 failed`. The failures are expected and identify
production behavior, not test-setup errors:

1. `best_dockq` is used instead of `GlobalDockQ`.
2. First-interface fields are promoted to multichain result fields.
3. The wrapper does not accept or pass an explicit CPU count.

Existing validation command:

```bash
python -m pytest -q \
  tests/test_score_bijective_benchmark_models.py \
  tests/test_model_output_integrity.py \
  tests/test_dockq_runtime.py
```

Observed result: `16 passed`. This supports the conclusion that the canonical
bijective/cross-interface path and its existing failure-state contracts are
already the preferred repair target.

Post-Green focused validation:

```bash
python -m pytest -q \
  tests/test_model_output_integrity.py \
  tests/test_standardized_evaluator.py \
  tests/test_score_bijective_benchmark_models.py \
  tests/test_dockq_parser_contract_regression.py \
  tests/test_dockq_runtime.py
```

Observed result: `45 passed`.

PRISM-side helper validation:

```bash
python -m pytest -q \
  tests/test_dockq_parser_contract.py \
  tests/test_score_transformed_models_contract.py \
  tests/test_compare.py \
  tests/test_transformation.py \
  tests/test_optional_backends.py
```

Observed result: `29 passed` before the final parser/helper edge-case test;
the final parser/helper subset is `5 passed` and the complete command was
rerun after the changes.

Final PRISM regression and script-contract validation:

```bash
python -m pytest -q \
  tests/test_dockq_parser_contract.py \
  tests/test_benchmark.py \
  tests/test_compare.py \
  tests/test_score_transformed_models_contract.py \
  tests/test_scoring_slurm_contract.py \
  tests/test_transformation.py \
  tests/test_optional_backends.py \
  tests/test_prodigy_ranker.py
```

Observed result: `49 passed`. The standalone parser/helper/Slurm subset was
also run independently and passed `16` tests before the final no-align
regressions. Python compilation and
`bash -n` checks passed for the changed entry points.

Final PRISM-prescript focused validation:

```bash
python -m pytest -q \
  tests/test_dockq_parser_contract_regression.py \
  tests/test_model_output_integrity.py \
  tests/test_score_bijective_benchmark_models.py \
  tests/test_score_comparison_models.py \
  tests/test_standardized_evaluator.py
```

Observed result: `53 passed`; changed Python and Slurm entry points compiled
and passed shell syntax checks.

Final post-correction validation after the low-level CPU-default change:

```bash
# PRISM
pytest -q tests/test_dockq_parser_contract.py tests/test_benchmark.py \
  tests/test_compare.py tests/test_score_transformed_models_contract.py \
  tests/test_scoring_slurm_contract.py tests/test_transformation.py \
  tests/test_optional_backends.py tests/test_prodigy_ranker.py

# PRISM-prescript
pytest -q tests/test_dockq_parser_contract_regression.py \
  tests/test_model_output_integrity.py tests/test_score_bijective_benchmark_models.py \
  tests/test_score_comparison_models.py tests/test_standardized_evaluator.py
```

Observed results: `50 passed` in PRISM and `57 passed` in PRISM-prescript.
Changed Python files compiled, every notebook/benchmark Slurm wrapper passed
`bash -n`, and scoped `git diff --check` passed. Pinned DockQ help exposes
`--json`, `--mapping`, `--no_align`, and `--n_cpu`.

Independent Terra validation returned `PASS WITH CAVEATS`: the five final
scope/status claims were supported by the specified source files and focused
tests. Its caveats are intentional: real production scoring and corrected
BM55/EDA regeneration were not run in this contract-correction phase.

Retained raw-JSON audit over 14,470 records:

- parser errors: `0`
- `GlobalDockQ` mismatches: `0`
- corrected scores outside `[0,1]`: `0`
- records with multiple interfaces: `4,696`

The representative retained record changed from the invalid sum `13.0026` to
the valid `GlobalDockQ` value `0.684348921892922`.

The repository-wide `python -m pytest -q` command is not a valid gate in the
current checkout: collection reports duplicate module-name conflicts,
package-relative import failures, nested-test conftest collisions, and the
run was interrupted before completion. Those failures are pre-existing
repository-layout issues, not failures of the scoped DockQ contract tests.

## Plan of work after review

1. Implement the smallest parser change in both parser copies: select
   `GlobalDockQ` as the bounded complete-mapping scalar, preserve `best_dockq`
   only as an explicitly named diagnostic sum, and do not fall back to an
   interface `DockQ` value when `GlobalDockQ` is absent.
2. Prevent ambiguous single-interface fields from being used as grouped
   multichain metrics.
3. Add explicit `n_cpu` propagation and keep it equal to the Slurm allocation;
   retain exact argv and raw-JSON provenance where the ad-hoc path is kept.
4. Run the focused tests, then the relevant full suite.
5. Replay the retained raw JSON case with both old and corrected parsing;
   verify no corrected DockQ or CAPRI score is outside `[0, 1]`.
6. Use the canonical bijective scorer for receptor-ligand cross metrics and
   retain `GlobalDockQ`, internal-best, cross-best, and cross-mean as separate
   columns. The transformed-run helper now implements the same separation for
   its assembled models; the staged benchmark path remains canonical for
   confirmatory BM55 scoring.
7. Regenerate downstream EDA only from rows with the appropriate successful
   score status, retaining missing and failed rows in audit tables.
8. Separately audit or route `score_comparison_models.py` through the
   canonical no-align validation before using it for confirmatory results.

## Acceptance criteria

- Every complete-mapping scalar DockQ reported by the corrected path is in
  `[0, 1]` and equals raw `GlobalDockQ` when that field is present.
- `best_dockq` is never exposed as an unqualified `dockq` field.
- Multichain component fields are not silently taken from an arbitrary first
  interface.
- The exact DockQ command and CPU setting are retained in evidence.
- Receptor-ligand cross-interface metrics are separate from global/internal
  metrics.
- Failed or partial records remain explicit and are excluded from complete
  quality summaries.
- Refined and transformed model results remain separate.
- Existing raw JSON, benchmark inputs, validated results, and unrelated dirty
  worktree changes are untouched.

## Decision Log

- **2026-09-17 — Use TDD before production edits.** The new regression tests
  must fail against the current ad-hoc parser before implementation begins.
- **2026-09-17 — Prefer the existing canonical benchmark scorer.** Its focused
  tests already cover complete bijections, requested cross interfaces, and
  explicit failure states; duplicate scoring logic should not be invented.
- **2026-09-17 — Do not clamp invalid DockQ values.** Values above one are
  evidence of wrong aggregation and must be corrected at the source.
- **2026-09-17 — Reconcile duplicate parsers.** Fixing only `src/eval/dockq.py`
  would leave `benchmark/scripts/dockq.py` and the production PRISM runtime
  parser scientifically inconsistent.
- **2026-09-17 — Keep no-align follow-up separate.** The canonical scorer
  validates residue correspondence before `--no_align`; the comparison path
  requires its own routing fix before confirmatory use.
- **2026-09-17 — Preserve exact single-pair scoring provenance.** When a
  persistent JSON directory is requested, retain its path, SHA-256, exact
  argv, mapping-validation status, and CPU setting in the score row.
- **2026-09-17 — Treat valid-but-unscored JSON as a first-class status.** A
  successful DockQ process is not sufficient for a quality score; complete
  mappings without usable `GlobalDockQ` remain auditable but are excluded from
  scored totals.
- **2026-09-17 — Propagate CPU settings at submission boundaries.** The
  notebook Slurm wrappers derive `DOCKQ_N_CPU` from the actual allocation and
  pass it explicitly, preventing a hidden DockQ default from diverging from
  the requested resources.
- **2026-09-17 — Require strict no-align correspondence.** All low-level
  helpers that expose `--no_align` require an explicit bijective chain map and
  matching residue numbering, insertion codes, and residue identities before
  launching DockQ.
- **2026-09-17 — Label pairwise cross-only recovery explicitly.** A
  pairwise recovery after a complete-mapping crash may retain requested
  receptor-ligand component scores, but it must use the cross-only scope and
  `scored_cross_only` status; it must never imply a complete-complex
  `GlobalDockQ` score.
- **2026-09-17 — Default low-level DockQ calls to one CPU.** Direct helper
  calls pass `--n_cpu 1` unless a caller supplies an explicit resource-matched
  override; this prevents hidden parallelism outside Slurm-aware wrappers.

## Validation

Authoritative evidence paths and sources:

- `tmp/prism_all_pipelines/runs_gtalign/20260916T034139Z/bm55_full/full/gtalign_gpu/processed/dockq_irmsd/`
- `src/eval/dockq.py`
- `benchmark/scripts/score_bijective_benchmark_models.py`
- `docs/exec-plans/20260717-dockq-irmsd-scoring-contract-audit.md`
- `docs/exec-plans/20260726-bm55-canonical-scoring.md`
- [DockQ v2 paper](https://academic.oup.com/bioinformatics/article/40/10/btae586/7796530)

## Outcomes & Retrospective

The parser, CPU, no-align, row-provenance, and transformed-helper
cross-interface corrections are implemented and validated by the scoped tests
and retained-JSON audit. Corrected benchmark scores and EDA have not been
regenerated; the canonical staged scorer remains the required next path for
confirmatory BM55 results. Low-level omitted-CPU calls now use one CPU, while
Slurm wrappers preserve the allocation explicitly. No Slurm job, raw model,
native structure, or existing scoring artifact was modified. The plan remains
open only for the separately authorized corrected-score and EDA regeneration.
