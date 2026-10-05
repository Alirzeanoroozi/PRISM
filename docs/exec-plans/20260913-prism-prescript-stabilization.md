# Assess and stabilize comparable PRISM pipeline variants

This ExecPlan is a living document. Keep `Progress`, `Surprises &
Discoveries`, `Decision Log`, and `Outcomes & Retrospective` current while
working. The selected project is the single underlying
`PRISM-prescript`/`prism-refactoring` project; the two VALAR evidence
namespaces remain separate and their provenance is not rewritten.

## Purpose / Big Picture

Produce an auditable current-state report, promotion-ready isolated pipeline
changes, controlled alignment/ranking/refinement comparisons at the exact
small and large template scales, and a portable interactive notebook. The
canonical `/scratch/rshadi25/GitHub/PRISM-prescript` checkout, raw data,
canonical outputs, environments, and protected files remain untouched. Final
claims must distinguish IMPLEMENTED, TESTED, VALIDATED, and REVIEWED.

Hourly monitoring rule: at each hourly checkpoint while this goal remains
active, inspect the newest VALAR status and authoritative run/job artifacts,
reconcile scheduler state and durable manifests, record the current state and
exact next action in the evidence/plan trail, and continue only within the
selected project and bounded goal. This interface exposes no native recurring
prompt/scheduler, so the rule is durable and must be applied at each available
continuation/checkpoint rather than represented as an unbounded background
process.

## Progress

- [x] Recover both VALAR namespaces, canonical guidance/memory, current plan,
      stable-pipeline contract, worktrees, notebooks, and recent evidence.
- [x] Reconcile job `1657005` and the `run-2688...` job/state mismatch.
- [x] Freeze exact template panels and build a provenance/no-drop comparison
      artifact inventory.
- [x] Review or implement isolated alignment/ranking/stage-contract changes
      with focused tests and small deterministic smoke coverage.
- [x] Run permitted Slurm benchmarks at both panel sizes, or preserve explicit
      BLOCKED/UNKNOWN states with scheduler evidence.
- [ ] Complete the same-set no-ranking, baseline-ranking, and PRODIGY
      comparison with quality and end-to-end timing; the current package
      preserves explicit `NOT_RUN` rows because this evidence does not exist.
- [x] Update and headlessly validate the authoritative PRISM notebook.
- [x] Produce the requested reports/matrices in an isolated project artifact
      area and update framework status without merging or promoting code.
- [ ] Promote a shared AlignmentResult/stage-manifest contract or call any
      pipeline variant stable; current work only identifies the required
      isolated implementation boundary.

## Surprises & Discoveries

- Observation: the PRISM-prescript status is `PARTIAL / NOT_PROMOTION_READY`;
  job `1657005` was canceled before execution and has no durable result.
  Evidence: `status/prism-prescript/26-09-28-LATEST.md` and its Step 5 evidence.
- Observation: prism-refactoring's direct run manifest points to retry job
  `1657894`, while the delegation summary reconciles the same run to
  `1657891`/`RELEASED/FAILED`.
  Evidence: `evidence/prism-refactoring/run-2688ab679e0f47409f2a19e182483a75/manifest.json`
  and `qwen-delegation-summary.json`.
- Observation: the canonical checkout is dirty and the requested
  `.github/copilot-instructions.md` is absent at the stated path.
  Evidence: canonical `git status` and `.github` listing.
- Observation: the previously selected `/scratch/tmp` snapshot has an
  incomplete Git object database (`fatal: bad object HEAD`), so it cannot be
  treated as a valid worktree without independent source provenance.
  Evidence: read-only `git -C /scratch/tmp/prism-prescript-refactoring-light.zypqDR
  status`.
- Observation: live Slurm accounting supersedes the stale 2026-09-12
  orientation text for PRISM-prescript job `1657005`: it is terminal
  `CANCELLED by 1365446`, with zero runtime and no durable output.
  Evidence: `sacct -j 1657005` on `ai01` at the 2026-09-13 checkpoint.
- Observation: both jobs associated with prism-refactoring
  `run-2688ab679e0f47409f2a19e182483a75` are terminal failures: `1657891`
  failed with exit `127:0`, and retry `1657894` failed with exit `1:0`.
  The direct manifest remains stale and points at `1657894` as pending;
  accounting and the delegation summary establish that no worker report was
  produced. Evidence: `sacct -j 1657891,1657894`, the run manifest/logs, and
  `evidence/prism-refactoring/qwen-delegation-summary.json`.
- Observation: scheduler controllers are currently up, but the live AI
  allocation is fragmented and unrelated user jobs occupy resources. The
  right-sized exact-panel TMalign jobs were admitted in available capacity;
  no existing user job was changed. Evidence: `gpu-status.sh`, `scontrol
  ping`, `sinfo`, `squeue -u rshadi25`, and `sacct`.
- Observation: the current on-disk template inputs have distinct exact
  denominators: `new_template/template/checked_templates.txt` contains
  19,948 entries, while `calculated_templates.txt` contains 19,062; the
  historical 2026-07-18 alignment panel contains 946 and the historical
  2026-07-19 20K run used a separate 19,855-entry list. These must not be
  silently conflated. Evidence: line counts and SHA-256 manifests.
- Observation: the isolated comparison package contains 59 pipeline rows, 60
  performance/optimization rows, and 13 ranking rows, plus exact panel,
  no-drop, source, and reconciliation manifests. Missing clean comparisons
  remain explicit `NOT_RUN`/`NOT_COMPARABLE` rows.
  Evidence: `evidence/prism-prescript/20260913-comparison/`.
- Observation: the nine-cell artifact-driven notebook passed Jupyter/JSON/AST
  checks and portable headless execution (`rc=0`) with
  `PRISM_ARTIFACT_ROOT` set to the package directory; it launched no pipeline
  arm. The shell `python3` lacks `nbformat`, so the verifier used the Jupyter
  `nbconvert` environment plus standard-library JSON/AST checks. The final
  notebook also loads the four provenance JSON artifacts, and its outcome
  filter partitions the 59 matrix rows into 15 successful and 44
  failed-or-missing rows. Evidence: `evidence/prism-prescript/20260913-comparison/NOTEBOOK_VALIDATION.json`
  and disposable `/scratch/tmp/prism-notebook-headless-6uzhch/executed-final.ipynb`.
- Observation: the isolated GTalign GPU runner completed the 946-prefix as
  job `1659356` and the materialized 19,062-entry calculated panel as job
  `1659358`; the all-entry 19,948 checked panel failed closed as job `1659357`
  when staging found a missing/empty interface. Submission `1659355` failed
  immediately because `/etc/bashrc` was sourced after `set -u`; that runner
  defect was fixed and the failed attempt remains preserved.
- Observation: the isolated USalign runner completed the same 946-prefix and
  materialized calculated panel as jobs `1659360` and `1659361`, respectively,
  with eight independent CPU workers and `USalign QUERY REFERENCE -outfmt -1
  -m -`. These are alignment-only measurements; downstream candidate-set,
  transform/filter, refinement, and quality parity remain unvalidated.
  Evidence: the two runner scripts, `sacct`, and the disposable run roots.
- Observation: final accounting is terminal and consistent with the durable
  run manifests: `1659356`/`1659358` and `1659360`/`1659361` completed; `1659355`
  and `1659357` failed for their recorded reasons. The package now contains
  59 matrix rows, 60 performance rows, 13 ranking rows, and five explicitly
  distinct panel definitions, including the current checked-prefix 946 hash.
- Observation: current exact TMalign references completed as jobs `1659448`
  (946-prefix: 1,892/1,892 successful sides, 12.31 s total) and `1659449`
  (19,058 materialized templates: 38,116/38,116 successful sides, 157.18 s
  total). The records include three-row transform validation and explicit
  per-side no-drop status, but remain alignment-only.
- Observation: final package/schema checks passed, all notebook code cells
  parsed as Python AST, and the generated notebook executed headlessly with
  no error outputs in disposable `/scratch/tmp` (`rc=0`) using the explicit
  non-interactive `Agg` backend and writable Matplotlib cache. The default
  host backend remains environment-sensitive; CSV loading, filtering,
  provenance output, and the dependency-light views remain functional without
  optional plotting support.
- Observation: the isolated shared-contract candidate now includes a
  provider-neutral event schema, alignment-result adapter, additive
  `prism.py` stage-event integration, and guarded optional-backend imports.
  Its focused contract/adapter/import suite passed 28 tests and its complete
  isolated suite passed 31 tests. The default parser imports on the pinned
  clean baseline; missing optional backends fail explicitly when requested,
  and unknown/rejected alignment statuses remain explicit non-success records
  even when score or mapping fields are present. This remains TESTED isolated
  candidate evidence, not canonical pipeline validation or promotion. Evidence:
  `/scratch/tmp/prism-prescript-contract-20260913` and
  `20260913-comparison/CANDIDATE_CHANGES.md`.
- Review finding: the contract candidate remains partial for end-to-end
  no-drop provenance because it has no event/parent correlation, does not wrap
  every pipeline stage, and its adapter is not wired into provider writers.
  Optional run/attempt fields are present, but no bridge to canonical
  benchmark foreign keys exists. The independent Luna High review also
  identified and prompted the corrected status-classification bug; no claim of
  VALIDATED or promotion follows from the passing suite.
- Observation: the candidate contract now carries optional `run_id` and
  `attempt_id` fields and retains provider-specific numeric metrics such as
  MultiProt RMSD without treating them as TM-score. Event/parent correlation,
  complete stage wrapping, and provider-writer integration remain open.
- Review finding: canonical `src/provenance/run_evidence.py`,
  `benchmark/scripts/investigation_lineage.py`,
  `benchmark/scripts/investigation_contracts.py`, and
  `benchmark/scripts/pipeline_completion_contract.py` remain authoritative
  for run identity, benchmark lineage, evaluator/ranking contracts, and
  completion classification. The isolated candidate is an additive diagnostic
  projection; it must not replace those layers or be transplanted wholesale
  into the dirty canonical `prism.py`.
- Observation: the rebuilt comparison notebook was freshly executed
  headlessly from the package with `MPLBACKEND=Agg`, zero cell errors, and no
  interactive-session dependency. Evidence: disposable
  `/scratch/tmp/prism-notebook-headless-final-sY3Nlw/executed.ipynb` and
  `20260913-comparison/NOTEBOOK_VALIDATION.json`.
- Observation: a bounded Luna High worker added a benchmark-side
  `build_alignment_event_ledger.py` projection and six new regression tests
  to the isolated contract worktree. The ledger emits one deterministic JSONL
  alignment event per resolved manifest row, retains missing/not-run states,
  and rejects missing identity or duplicate candidate IDs before writing.
  The focused contract/adapter/import/ledger suite passed 34 tests and the
  complete isolated suite passed 37 tests. This is TESTED synthetic-manifest
  evidence only; the ledger is not connected to canonical provider writers or
  lineage and does not upgrade any pipeline arm to VALIDATED.

## Decision Log

- Decision: use an isolated copy/worktree under `/scratch/tmp` for all project
  edits and derived benchmark artifacts; do not edit the canonical checkout.
  Rationale: the canonical worktree contains extensive user changes and the
  project contract forbids overwriting canonical results or protected files.
  Date/Author: 2026-09-13 / Codex.
- Decision: preserve `evidence/prism-prescript/` and
  `evidence/prism-refactoring/` as distinct provenance namespaces while
  treating their project identity as one for synthesis.
  Rationale: the user requires unified project handling without provenance
  collapse.
  Date/Author: 2026-09-13 / Codex.
- Decision: never pool MultiProt's RMSD-derived proxy with TM-scores and never
  call candidate reduction alone a speedup.
  Rationale: existing evidence documents incompatible score semantics and
  incomplete end-to-end timing/quality evidence.
  Date/Author: 2026-09-13 / Codex.
- Decision: do not automatically resubmit `1657005` or retry the failed
  worker review. The former was a canceled, never-executed stale validation
  attempt; the latter already has two terminal failures and no scientific
  report. Reuse existing evidence and make any new execution a separately
  manifested bounded run after exact inputs and safe capacity are frozen.
  Date/Author: 2026-09-13 / Codex.
- Decision: treat retained 946/19,855 measurements as observational. Record
  the new exact 946-prefix and 19,062 materialized-panel GTalign/USalign runs
  as alignment-only evidence, while keeping current exact TMalign and all
  downstream transformation/ranking/refinement/evaluator arms `NOT_RUN` until
  a single query/template/evaluator manifest is frozen. This prevents a
  superficially successful but non-causal end-to-end benchmark.
  Date/Author: 2026-09-13 / Codex.
- Decision: submit only one-GPU GTalign alignment-only jobs for the exact
  current checked-panel prefix/all-entry runs, with explicit staging/search
  timing and hashes. Also submit the bounded eight-worker USalign counterpart
  only after live resource inspection. Do not infer downstream quality,
  candidate-set agreement, or promotion from alignment-only jobs; reconcile
  their manifests before adding results.
  Date/Author: 2026-09-13 / Codex.
- Decision: add the exact TMalign reference as a separately manifested
  alignment-only arm rather than infer it from the canceled Step 5 attempt.
  Keep its complete pairwise output and three-row transform checks separate
  from downstream transformation/filtering and quality claims.
  Date/Author: 2026-09-13 / Codex.
- Decision: after the user explicitly directed that Qwen/Copilot not be used,
  use bounded Codex Luna High read-only audits for package/evidence, notebook,
  downstream-contract, and architecture review. No worker launched Slurm,
  changed project files, or made a cross-workstream scientific decision.
  Date/Author: 2026-09-13 / Codex.
- Decision: continue the isolated contract repair with Codex Luna High only;
  Qwen/Copilot are excluded from this goal. Keep the optional-import repair,
  alignment adapter, and stage-event integration isolated until an independent
  review and downstream matched runs establish whether they can be promoted.
  Date/Author: 2026-09-13 / Codex.
- Decision: preserve the adapter's fail-closed semantics for unknown and
  explicitly non-success statuses, even when legacy records contain scores or
  residue mappings. Add regression coverage for optional-backend guards and
  those status classes before considering any integration. Date/Author:
  2026-09-13 / Codex, following Luna High review.
- Decision: keep the canonical provenance/lineage/evaluator/completion layers
  as the integration authority. If integration is later approved, emit a
  separate alignment-event projection only after the benchmark-side resolver
  has exact dataset/pair/template/orientation identities and alignment hashes;
  do not alter provider writers or `transformation.load_alignment()` in this
  candidate. Date/Author: 2026-09-13 / Codex, following Luna High review.
- Decision: accept the isolated manifest-to-alignment-event ledger as the
  next safe integration boundary, but only as a tested additive projection.
  It may be used for future reconciliation after exact benchmark foreign keys
  and raw tool-output hashes are resolved; it must not replace canonical
  `run_evidence`, lineage, evaluator, or completion contracts. Date/Author:
  2026-09-13 / Codex, following bounded Luna High implementation and local
  review.

## Outcomes & Retrospective

Recovery is reconciled through the live 2026-09-13 checkpoint. An isolated,
provenance-linked current-state package and artifact-driven notebook are now
generated and validated. No project code, canonical result, raw dataset,
environment, or protected file has been changed by this goal. Exact
alignment-only GTalign, TMalign, and USalign arms are complete within their
manifests;
the all-entry checked-panel staging failure and all downstream
transformation/ranking/refinement/evaluator gaps remain explicit. Candidate
throughput is not claimed where the provider output was not normalized, and
retained runs do not satisfy exact current panel/evaluator parity. The package
records the safe next benchmark protocol.

## Context and Orientation

Framework root: `/home/rshadi25/valar-agent-framework`.

Canonical project: `/scratch/rshadi25/GitHub/PRISM-prescript`.
The maintained path is `prism.py` plus `src/`, with benchmark scripts under
`benchmark/scripts/`, current notebooks under `notebooks/`, and isolated
artifacts under `tmp/agent/`. Stable defaults and operational commands are in
`docs/STABLE_PIPELINE.md`. The active project planning state is under
`.planning/`; the canonical checkout is already dirty and must be treated as
user-owned state.

The framework status indexes are:
`status/prism-prescript/26-09-28-LATEST.md` and
`status/prism-refactoring/26-09-13-LATEST.md`. Their authoritative evidence roots are
`evidence/prism-prescript/` and `evidence/prism-refactoring/`.

## Plan of Work

First reconcile durable status, scheduler/accounting state, and source
provenance. Then inspect existing isolated candidates and retained benchmark
artifacts to avoid repeating conclusive work. Establish a new isolated
artifact root and exact template-panel manifests. Add only the smallest
contract/test changes needed for stage-level comparability and no-drop
accounting. Prepare bounded Slurm runs using live resource inspection and
internal parallelism where validated. Analyze outputs with explicit
NOT_RUN/BLOCKED/FAILED/NOT_COMPARABLE/COMPLETED states. Update the notebook to
consume stored artifacts, validate it headlessly, review critical diffs, and
write the requested reports. Do not merge, push, promote, or alter canonical
project state.

## Concrete Steps

1. Working directory `/home/rshadi25/valar-agent-framework`: read status,
   evidence pointers, and mismatch artifacts; record hourly checkpoint state.
2. Working directory `/scratch/rshadi25/GitHub/PRISM-prescript`: read source,
   planning, stable-pipeline, benchmark, notebook, and isolated-run metadata
   read-only; classify canonical dirtiness.
3. Working directory `/home/rshadi25/valar-agent-framework`: create a
   timestamped evidence index and isolated-worktree provenance record.
4. Working directory `/scratch/tmp/<isolated-project-root>`: implement and
   test contract/notebook changes only after a failing public-interface test
   establishes each needed behavior.
5. Login node: run `gpu-status.sh`, `hostname`, Slurm identity, `squeue`,
   `sinfo`, scheduler/account/QOS checks; compute nodes only for heavy tests,
   alignment, scoring, and notebook execution.
6. Compute node via `sbatch`: run one bounded job or controlled array per
   benchmark arm with logs under the isolated run's `logs/` directory and
   machine-readable manifests for every stage and candidate disposition.
7. Working directory `/scratch/tmp/<isolated-project-root>`: generate the
   requested CSV/JSON/Markdown/notebook deliverables and validate schemas,
   imports, AST, headless execution, tests, and diff cleanliness.

## Validation and Acceptance

- Job/run reconciliation is supported by both scheduler/accounting output and
  durable manifests; absent scheduler data remains UNKNOWN.
- Exact small and large template counts, paths, hashes where practical,
  exclusions, and duplicate/leakage notes are recorded.
- Every comparison row has source revision/worktree, executable/version,
  parameters, template/query identity, stage counts, timing, status, and
  failure/no-drop reasons.
- USalign uses the transform-producing invocation equivalent to
  `-outfmt -1 -m -`; its output is tested independently before comparison.
- MultiProt proxy scores remain separate from TM-score contracts.
- PRODIGY conclusions include ranking overhead, end-to-end wall time, quality,
  regret, and no-contact behavior; no speedup claim is made from load
  reduction alone.
- Focused tests, deterministic smoke, real Slurm evidence where required,
  machine-readable manifests, quantitative output checks, and independent
  critical-region review are present; missing checks are explicit.
- Final response names the fastest defensible, best-quality, and recommended
  trade-off configurations plus evidence required before promotion.

## Idempotence and Recovery

All derived data goes to a new isolated run root and is append-only or
content-addressed where practical. Existing evidence is reused only when
inputs, code, evaluator, and provenance match. A failed or canceled Slurm
task remains in the matrix and is never inferred complete from sibling output.
If scheduler state is unavailable, stop new expensive submissions, record
UNKNOWN/BLOCKED, and resume from the same plan/checkpoint after recovery.
Uncertain temporary artifacts are quarantined or preserved; no destructive
cleanup is performed.

## Artifacts and Notes

- Living plan: `docs/exec-plans/prism-prescript/20260913-prism-prescript-stabilization.md`.
- Deliverables: isolated project artifact root, with copies/indexes linked from
  the corresponding framework evidence namespace.
- Hourly rule: re-read both latest status indexes and the newest authoritative
  manifests/logs, reconcile state, write the next action, then continue.
- Canonical status fingerprints, raw data, and result roots are read-only.

## Interfaces and Dependencies

The project uses Python 3.11 in
`/home/rshadi25/.conda/envs/gtalign_env/bin/python`, external alignment tools,
optional NACCESS/FreeSASA, Rosetta 2022.42, optional PyRosetta/FiberDock,
DockQ 2.1.3, and Slurm. The framework uses JSON/TSV/CSV evidence manifests
and status indexes. Required notebook dependencies must be imported or
parameterized explicitly; no author-session state is allowed.
