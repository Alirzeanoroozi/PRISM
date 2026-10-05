# PRISM validation report

## IMPLEMENTED

Canonical `prism.py` supports TMalign, GTalign CPU/GPU, MultiProt, NACCESS,
FreeSASA, external Rosetta, PyRosetta, FiberDock, and optional ranking paths.
USalign and PRODIGY implementations exist in isolated candidates.

## TESTED

The latest retained evidence includes canonical focused tests (22/22), isolated
USalign tests (20/20), isolated PRODIGY tests (9/9), full canonical tests
(344 passed, 6 skipped), a real USalign transform probe, real PRODIGY success
and no-contact probes, current stable execution smokes, and the generated
comparison notebook's static/portable checks.

The notebook passed nbformat/AST validation, loaded all three CSVs plus the
four provenance JSON artifacts through `PRISM_ARTIFACT_ROOT`, and passed a
direct successful-versus-failed-or-missing filter check (15/44 of 59 matrix
rows). Kernel-backed headless execution completed with zero cell errors using
the explicit non-interactive `Agg` plotting backend; the validation record is
`NOTEBOOK_VALIDATION.json`. The default login-node Matplotlib backend is
environment-sensitive, so reproducible headless execution should set
`MPLBACKEND=Agg` and a writable `MPLCONFIGDIR`.

The isolated shared-contract candidate also passed 34 focused contract/
adapter/import/ledger tests and 37 tests in its complete clean-baseline suite.
The benchmark-side ledger is a separate, additive projection from resolved
candidate manifests; this is software evidence only and does not validate a
provider, a full pipeline, or scientific output.

## VALIDATED (bounded scope only)

USalign transform parsing/application, PRODIGY state/failure preservation,
1gte geometry-only attrition accounting, the historical 946 alignment pair
coverage, the exact current checked-prefix/calculated-panel alignment-only
timings for GTalign GPU, TMalign, and USalign, and historical 19,855-run
output/quality subsets have quantitative evidence within their stated scopes.
None establishes general scientific superiority or promotion readiness.

## REVIEWED

The canonical guidance/memory, source call paths, isolated worktrees, benchmark
manifests, notebooks, and both VALAR namespaces were reviewed. The review found
missing uniform score/status/no-drop contracts, incomplete refiner observability,
and known CPU/GPU/panel/evaluator mismatches. Independent Luna High review
confirmed that canonical `run_evidence`, `investigation_lineage`,
`investigation_contracts`, and `pipeline_completion_contract` remain the
authoritative layers; the isolated candidate is TESTED only and is not a
replacement.

## BLOCKED / UNKNOWN

The current exact panels have matched GTalign GPU/TMalign/USalign alignment-only
runs, but no clean matched full-pipeline comparison through transformation/filtering,
ranking, refinement, and evaluation; current exact TMalign timing is also
present only for alignment-stage runs. PRODIGY lacks end-to-end timing and quality/regret evidence, and no
ranking top-k matrix is complete. Job 1657005 was canceled without execution;
run-2688 had two failed worker attempts without a scientific report. These
remain explicit in the matrix.

Successful Slurm or process return codes are not being promoted to scientific
validation without stage counts, provenance, quantitative outputs, and review.
