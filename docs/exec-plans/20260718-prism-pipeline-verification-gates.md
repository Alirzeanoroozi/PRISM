# PRISM Pipeline Verification Gates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` to implement this plan task-by-task. Do not use parallel agents unless the user explicitly authorizes delegation. Stop at every go/no-go gate before submitting the next Slurm scale.

**Goal:** Replace the invalid v3 completion and scoring path with fail-closed, source-gated, full-complex PRISM verification, then prove each alignment/refinement arm on matched inputs before any benchmark-quality claim is released.

**Architecture:** The work is divided into independent gates: immutable baseline evidence, structured stage/completion status, refined full-pose integrity, bijective scoring and provenance, template leakage control, published-filter fidelity, backend regression/load testing, and staged HPC validation. Existing provenance, evaluator, and source-gate modules are reused as the authority; ad-hoc scoring scripts are retained only as negative controls.

**Tech stack:** Python 3.11 from `gtalign_env`, pytest, Bash/Slurm, DockQ 2.x from `/scratch/tmp/prism-dockq-env`, Biopython, TM-align, GTAlign 0.19.00, external Rosetta 2022.42, and the opt-in PyRosetta 2026.3 adapter.

## Purpose / Big Picture

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` current during implementation.

The user should ultimately receive a verification bundle in which every reported model can be traced from a benchmark row and eligible template through alignment, transformation, refinement, full-complex validation, explicit model-to-native chain mapping, raw DockQ JSON, grouped iRMSD, and aggregation. A cancelled or partial run must never be called complete, a transformation half must never be accepted as a docking model, and a target-derived/self-homologous template must never enter a confirmatory denominator unnoticed.

Global constraints:

- Work only in `/scratch/rshadi25/GitHub/PRISM-prescript`; do not read PDF files.
- Treat `benchmark/originals/`, `new_template/template/`, retained run roots, and validated results as read-only. The first derived verification assets go to `tmp/agent/20260718-prism-verification/`; retries add a numbered `retry-N/` child.
- Preserve the invalid v3 run unchanged as a negative control. Supersede it through a claim ledger; do not delete or rewrite it.
- Use `gtalign_env` for pipeline-side Python and `/scratch/tmp/prism-dockq-env/bin/python` for DockQ.
- Run only lightweight unit/static checks on a login node. Execute alignment, Rosetta/PyRosetta, large hashing/sequence comparisons, and benchmark loops through Slurm.
- Keep TM-align plus external Rosetta as the baseline. GTAlign and PyRosetta remain diagnostic until matched downstream evidence passes.
- Keep pipeline completion, prediction production, scoreability, and prediction quality as separate fields.
- Do not use native-derived scores for candidate selection. Report invalid/missing models as explicit rows with null metrics.
- Do not open the unresolved 17-row source-gate cohort. The largest allowed confirmatory cohort is the currently frozen 240 strict rows unless authoritative resolution changes the policy.

## Progress

- [x] Gate 0: freeze a claim/evidence baseline and register the v3 outputs as invalid scientific evidence.
- [x] Gate 1: implement and test structured stage and cancellation-aware completion status.
- [x] Gate 2: enforce refined, assembled, full-complex model integrity before staging.
- [x] Gate 3: make bijective DockQ/grouped-iRMSD scoring provenance-complete and fail closed.
- [x] Gate 4: freeze the template panel and enforce self-hit/homology exclusion per benchmark row.
- [~] Gate 5: restore or explicitly gate published hotspot/contact filtering with validated assets.
- [x] Gate 6: regression-test GTAlign filtering and bound TM-align concurrency/resource use.
- [~] Gate 7: pass cancellation, no-prediction, and positive full-pose canaries.
- [ ] Gate 8: complete the matched 12-row four-arm pilot and adjudicate all failures.
- [ ] Gate 9: run the 240-row confirmatory comparison only after all prior gates pass.
- [~] Gate 10: publish the claim ledger, validation report, artifact hashes, and durable memory update.

## Surprises & Discoveries

- Observation: v3 Slurm jobs `1368048` and `1368049` were cancelled, but `status/exit.json` reports `scientific_status=completed` because the current `EXIT` trap treats any transformation PDB as success.
  Evidence: `tmp/agent/20260718-benchmark20k-v3/smoke/*/slurm-*.err` and `benchmark/scripts/submit_comparison_batches.sbatch`.
- Observation: `score_v3.py` scores `_L.pdb`/`_R.pdb` transformation halves. The representative score `0.955` measures the receptor's native internal `AB` interface and omits ligand chain `C`.
  Evidence: `tmp/agent/20260718-benchmark20k-v3/score_v3.py` and the retained `1cfzAE_1vfaAB_8lyz_o2_L.pdb` case.
- Observation: reusable strict pieces already exist in `benchmark/scripts/investigation_provenance.py`, `investigation_contracts.py`, `standardized_evaluator.py`, `stage_current_models_for_main_benchmark.py`, and `score_bijective_benchmark_models.py`.
- Observation: the 19,855-entry expanded list has complete interface PDB and interface-list assets, but at least the sampled expanded template `104lAB` lacks current JSON contact and hotspot assets even though legacy text assets exist. Published-filter asset coverage must be measured before claiming PRISM fidelity.
- Observation: baseline hash job `1368154` was submitted to KUTEM on 2026-07-18 and is pending by scheduler priority. It reads 178,217 retained files (approximately 2.1 GB) and writes only `tmp/agent/20260718-prism-verification/baseline/`.
  Evidence: `benchmark/jobs/pipeline_verification_baseline.sbatch` and Slurm job `1368154`.
- Observation: the split AI/COSBI evidence baseline completed successfully. AI job `1368194` hashed the v3 shard in 53 seconds (about 52 MB MaxRSS); COSBI job `1368195` hashed the smoke shard in 3 minutes 23 seconds (about 130 MB MaxRSS); dependent AI merge `1368198` completed in 2 seconds. The merged manifest has 179,225 artifact rows plus a header.
  Evidence: `tmp/agent/20260718-prism-verification/baseline/` and `tmp/agent/20260718-prism-verification/logs/`.
- Observation: the cancelled v3 run roots have no refined-output directories, so guarded staging correctly emits `no_current_models_found` rather than accepting transformation halves. Minimal integrity jobs `1368213` (AI, 3 s, 11.8 MB MaxRSS) and `1368214` (COSBI, 2 s, 11.8 MB MaxRSS) completed successfully.
  Evidence: `tmp/agent/20260718-prism-verification/gate2/retry-2/*/stage_manifest.csv` and `tmp/agent/20260718-prism-verification/logs/stage-integrity-136821{3,4}.*`.
- Observation: the frozen source policy remains `blocked_source_authority` with `confirmatory_run_authorized=false`; the new source gate rejects audit-only and strict rows until authorization changes.
- Observation: TM-align previously submitted all futures eagerly. The bounded queue now limits pending work to `2 * workers`; GTAlign acceptance requires both parsed TM scores and minimum matches without fabricating filtered-hit JSON.
- Observation: `hotspot_analysis` previously returned `True` unconditionally. `PRISM_FILTER_MODE` now distinguishes published protocol (fail-closed assets) from geometry-only experimental behavior.
- Observation: isolated completion canaries ran successfully on both AI (`1368299`, 3 s, 6.6 MB MaxRSS) and COSBI (`1368300`, 2 s, 6.7 MB MaxRSS); both produced identical cancellation/no-prediction/positive-status classifications.
- Observation: real retained GTAlign output staging with the cohort-manifest join produced 30 valid refined full-pose candidates and three explicit schema/filename failures. COSBI scoring retry `1368308` completed in 2:30 (549 MB MaxRSS): 22 models scored, 8 retained as explicit DockQ failures, and 58 global/component interface records were written as true TSV.
- Observation: AI modern-template preflight `1368312` completed in 3:08 (238 MB MaxRSS): all 19,855 listed IDs were syntactically valid, but only 19,005 were fully resolvable; 850 require missing assets. Summary: `tmp/agent/20260718-prism-verification/template-gate/preflight-ai/source_gate_summary.json`.
- Observation: legacy preflight `1368319` completed in 1:22 (152 MB MaxRSS) and covered all 19,855 listed contact profiles; legacy hotspot/contact filename intersections are also complete for the frozen list. A derived JSON tree was built in `protocol-assets-ai`; retry `1368334` adds a per-asset SHA-256 manifest.
- Observation: hash-bound protocol asset retry `1368334` completed in 2:15 (43.8 MB MaxRSS) with 39,710 manifest rows, two assets per each of 19,855 templates, and 64-character SHA-256 values for every derived file.
- Observation: `published_protocol` now loads those derived assets by template ID and records their hashes; missing/invalid assets fail closed. The geometry-only path remains explicitly labeled experimental.
- Observation: corrected GTAlign contract smokes completed on AI (`1368352`, 4 s) and COSBI (`1368353`, 3 s) with 1 CPU/2 GB/5 minutes. Each produced 2/2 TM-align and GTalign pairs, with no missing backend records. The first submission (`1368349`/`1368350`) intentionally failed because an empty workspace override skipped the seed copy; it is retained as a launcher-negative-control artifact.
- Observation: the minimum 25-pair alignment-only backend comparison completed on AI (`1368355`, 2 s) and COSBI (`1368356`, 2 s), each with 25/25 TM-align pairs, GTalign return code 0, and five output files. These are execution/parity evidence only; they contain no DockQ, iRMSD, refinement, or source-gate claim.
- Observation: the actual 25-template panel smoke completed on AI (`1368362`, 12 s, 9.9 MB MaxRSS) and COSBI (`1368363`, 5 s, 10.0 MB MaxRSS). The single 1fgnH query generated 50/50 TM-align and GTalign chain-pair records on each node; no backend-only pair was lost. This remains alignment-stage evidence, not a pipeline-quality claim.
- Observation: the frozen 946-template alignment panel completed on AI (`1368365`, 2:05, 76,224 KB MaxRSS) and COSBI (`1368366`, 1:41, 68,780 KB MaxRSS). The single 1fgnH query generated 1,892/1,892 TM-align and GTalign chain-pair records on each node, with zero backend-only records. This confirms execution/count reconciliation at the diagnostic panel scale only.
- Observation: a GTalign-only search over the complete interface directory completed on AI (`1368375`, 3:08, 655,356 KB MaxRSS) and COSBI (`1368376`, 2:59, 568,004 KB MaxRSS), but GTalign searched 39,878 of 39,896 files because 18 empty PDBs are outside the frozen list. This is retained as an explicit inventory diagnostic, not as frozen-panel coverage.
- Observation: the exact frozen-list GTalign search completed on AI (`1368398`, 5:57, 583,636 KB MaxRSS) and COSBI (`1368399`, 5:46, 576,832 KB MaxRSS). Both staged 19,855 templates/39,710 nonempty chain interfaces and raw GTalign reported exactly 39,710 structures searched, 1,829,409 total residues, and one output file. This closes the one-query/19,855-template backend count reconciliation, while remaining alignment-only evidence.
- Observation: the corrected all-chain preflight was rerun on AI (`1368389`, 4:10, 224,060 KB MaxRSS) and COSBI (`1368390`, 4:45, 237,676 KB MaxRSS). It writes 79,420 asset rows for 19,855 IDs and still reports 19,005 fully resolvable modern profiles/850 missing modern profiles; the missing set is unchanged because the 18 empty interface files are not in `final_list.txt`. Preflight now validates every chain-specific PDB instead of only the first.
- Observation: protocol-asset parity jobs AI `1368413` and COSBI `1368414` produced identical 19,855-row audits: 19,005 modern/derived disagreements, 850 modern-asset-missing rows, and zero exact rows. The complete derived legacy asset tree is therefore usable for a labeled canary but cannot be claimed semantically equivalent to the modern JSON assets.
- Observation: published-filter replay jobs AI `1368423` and COSBI `1368424` agreed on 22,704 retained-alignment rows: 4,259 published-protocol passes, 18,445 failures, and 96 rows across four templates with missing protocol assets. This is filter-replay evidence, not a confirmatory run.
- Observation: the first real one-row published canary exposed two launcher/integration defects: the launcher created empty template coordinate/profile directories, and the published threshold gate checked absent GTalign hotspot fields instead of loaded protocol assets. The launcher now links read-only template asset directories and records filter/cutoff parameters; `src/transformation.py` forwards protocol hotspots into the final threshold gate.
- Observation: diagnostic published-protocol canaries completed on AI `1368455` and COSBI `1368456` with minimum 8 CPU/40 GB/1 hour requests. Using explicitly relaxed, non-confirmatory thresholds and template `1bjaAB`, both reached transformation and external-Rosetta refinement: one paired transformation, terminal stage events, and `scientific_status=completed` (AI retained one structure; COSBI retained one structure). The canary audit from the preceding retry shows orientation o2 generated and o1 rejected; no score was claimed because the native complex PDB was not staged.
- Observation: the query-first staging parser now accepts the canary grammar and stages one COSBI refined model when joined to the frozen cohort manifest. Canonical scoring then failed closed with `non-bijective group: model='C' native='AB'`; the model is retained as an explicit `score_failed` row with no interface rows, so the canary does not support a quality claim.
- Observation: a chain-compatible single-chain diagnostic canary (`1BU6_O`/`1F3Z_A`, `1g60AB`) completed on AI `1368460` (3:07) and COSBI `1368461` (4:46), each with 8 CPU/40 GB/1 hour. Both staged one refined pose and passed the canonical scorer: one model row, two interface rows, DockQ 2.1.3, complete `OA:GF` mapping, raw JSON/hash, and grouped iRMSD (~15.16 Å). AI/COSBI refined-model hashes differ, as expected for independent Rosetta runs. The run used relaxed diagnostic thresholds and therefore is wiring/contract evidence only, not a quality estimate.
- Observation: the optimized full source-gate matrix completed on AI `1368528` (7:01) and COSBI `1368530` (7:14), using 8 CPU/40 GB/1 hour. Each processed 12 rows against all 19,855 frozen templates, emitted 1,270,720 similarity rows, found zero missing interface assets, and produced byte-identical similarity and eligible-template TSVs (SHA-256 `f6ab1745...cf3063` and `097ea0b7...f5258e6b`). The source policy remains `blocked_source_authority`, so every row has zero confirmatory-eligible templates despite 649,749 similarity-level eligible records; this is an intentional fail-closed result.
- Observation: the canonical scorer now hashes the native PDB before scoring so `score_failed` rows retain source-model and native-input provenance. COSBI scoring retry `1368545` was submitted with 2 CPU/8 GB/30 minutes to regenerate the retained positive-canary ledger without overwriting retry-1.
- Observation: scoring retry `1368545` completed in 2:29; it retained 30 model rows, 22 scored models, 8 explicit `dockq_runtime_buffer_dimensions` failures, and 58 interface rows. The regenerated failure adjudication has non-empty native hashes for all eight failures and `retry_authorized=false` for each. Gate 3 is complete as a provenance/fail-closed contract; the eight difficult structural cases remain excluded from quality denominators.
- Observation: normal-threshold published-protocol canaries AI `1368555` and COSBI `1368556` completed in 5 seconds with `completed_no_predictions` and zero paired transformations for the chain-compatible `1BU6_O/1F3Z_A`–`1g60AB` case. A replay-confirmed protocol-pass row (`1DQQ_CD/3LZT_`–`1axcAB`) was also tested on AI `1368560` and COSBI `1368561`; both completed in 5 seconds with `completed_no_predictions`. This confirms strict fail-closed behavior but does not satisfy the positive-full-pose canary requirement. Gate 7 remains partial until a strict published-protocol positive is found or source/filter authority is revised.
- Observation: a strict GTAlign candidate search based on prior replay records found no hit with both normalized GTalign scores ≥0.4 plus the published filter and match-percentage thresholds. New strict GTAlign canaries AI `1368567` and COSBI `1368568` therefore also completed in 5–6 seconds with `completed_no_predictions`; the current `parse_gtalign_hits()` contract intentionally requires both scores to meet threshold. This rules out treating prior reference-normalized-only GTalign records as strict positives.
- Observation: added `benchmark/scripts/validate_confirmatory_prism_run.py` and regression tests. Running it against the frozen source policy and the AI source-matrix eligible list returns `status=blocked`, `validator_rc=2`, with no template/hash inconsistencies; the sole blocker is the authoritative `blocked_source_authority` decision. This prevents accidental Gate 8/9 submission while preserving a deterministic readiness check for a future authorized policy.
- Observation: modern protocol JSON assets use chain-keyed hotspot records with three-letter residue names and numeric contact pairs, whereas the derived legacy tree uses flattened chain-qualified one-letter records. `src/template_filtering.py` now normalizes both forms, preserves `asset_format` and raw asset hashes, and exposes `hotspots_by_chain`; `src/transformation.py` selects the correct hotspot lists and reverses contact orientation for o2 instead of applying one flattened hotspot list to both partners. This fixes an integration/semantic bug but does not resolve the independent 19,005-row modern/derived parity disagreement.
- Verification step: added regression coverage for modern hotspot/contact normalization and chain/orientation selection, then ran the complete focused suite (`51 passed in 3.28s`), Python compilation, Slurm-script syntax checks, and `git diff --check`. The first three checks passed; `git diff --check` reported pre-existing trailing whitespace in unrelated `src/rosetta_refinement.py` and `src/surface_extract.py`, which was left untouched.

## Decision Log

- Decision: retain v3 as an immutable negative control and exclude it through a machine-readable superseded-run ledger.
  Rationale: its contradictory scheduler/status evidence and half-model scores are useful regression fixtures.
  Date/Author: 2026-07-18 / Codex plan.
- Decision: classify the four 946-template smoke runs as supported execution evidence only, and classify both v3 completed/quality claims as unsupported.
  Rationale: the smoke status records are internally consistent, whereas both v3 Slurm logs report cancellation despite completed status JSON.
  Date/Author: 2026-07-18 / Codex implementation.
- Decision: `completed` requires an observed normal pipeline return, terminal records for every expected stage, and output-integrity validation. A Slurm signal always wins over shell return code or file counts.
  Rationale: stage directories and partial files survived cancellation in v3.
  Date/Author: 2026-07-18 / Codex plan.
- Decision: only outputs from a refinement directory may become score candidates. Transformation halves are diagnostic intermediates, never models.
  Rationale: DockQ can assign a misleadingly high score to a native internal interface in a half-complex.
  Date/Author: 2026-07-18 / Codex plan.
- Decision: aligned DockQ with an explicit complete mapping is the default. `--no_align` is allowed only when `standardized_evaluator.validate_pdb_mapping()` proves exact chain, residue identity, and numbering correspondence.
  Rationale: strict no-align validity is model-specific.
  Date/Author: 2026-07-18 / Codex plan.
- Decision: the confirmatory template report hard-excludes exact self-hits and uses a preregistered primary homology exclusion of greater than 50% global sequence identity with at least 70% shorter-chain coverage. Sensitivity tables use 30%, 40%, 50%, 70%, and 100% cutoffs.
  Rationale: the local Tuncbag Markdown explicitly reports performance after removing templates above 50% sequence similarity and distinguishes 100% native/self templates. Recording multiple thresholds prevents dependence on one arbitrary cutoff.
  Date/Author: 2026-07-18 / Codex plan.
- Decision: publish two clearly separated transformation modes. `published_protocol` fails closed on missing hotspot/contact assets and applies the historical filters; `geometry_only_experimental` preserves current behavior but cannot support a claim of faithful PRISM reproduction.
  Rationale: silently treating `hotspot_analysis()` as true changes candidate selection and false-positive risk.
  Date/Author: 2026-07-18 / Codex plan.

## Outcomes & Retrospective

Gate 0 is complete. `tmp/agent/20260718-prism-verification/baseline/claims.tsv` records four supported bounded execution claims and two unsupported v3 claims. `artifact_manifest.tsv` contains 179,225 hashed evidence rows.

Gates 1 and 2 are complete. Current-pipeline runs now write opt-in stage JSONL and are classified fail-closed by signal, normal return marker, terminal stages, and refined-output presence. Staging accepts only paired refined output names, validates raw partner-chain integrity before symlinking, retains a manifest row for every candidate, and records model hash/chain/CA provenance. The cancelled v3 negative control stages zero models on both AI and COSBI.

Gates 3 and 4 are complete at the contract/provenance level. The scorer emits bijective mappings, raw JSON hashes, grouped iRMSD, and native/model hashes even on explicit failures; the full 12-row × 19,855-template source matrix reconciles byte-for-byte on AI/COSBI and records fail-closed authorization. Gate 5 remains unresolved: modern JSON and derived legacy filter assets disagree for 19,005 templates, with 850 modern profiles missing. Gate 6 is complete for its defined alignment-stage scope: GTAlign filtering and bounded concurrency pass focused tests, and the 25-template, 946-template, and exact 19,855-template Slurm diagnostics reconcile on both AI and COSBI.

Gate 7 has partial evidence: cancellation/no-prediction canaries and real published-protocol refinement canaries pass on both AI and COSBI. A chain-compatible single-chain diagnostic also passes canonical staging/scoring on both nodes, but all positive canary results use explicitly relaxed diagnostic cutoffs; they cannot authorize the matched pilot or confirmatory run.

The retained full-pose scoring canary still supplies separate downstream evidence: 22 models produced canonical DockQ/grouped-iRMSD rows, while eight difficult rows failed explicitly inside DockQ (`Buffer has wrong number of dimensions`) and remain excluded from quality metrics pending adjudication. The new published-protocol canary validates stage/refinement wiring only; it is not merged with those scoring claims.

Gate 10 is partially implemented as a fail-closed preflight: `validate_confirmatory_prism_run.py` verifies source authorization, row count, and template-list hashes before any confirmatory submission. Its current report is intentionally `blocked` because the authoritative source policy has not changed; no Gate 8 or Gate 9 artifacts have been generated.

Continuation verification (2026-07-18): modern JSON filter assets are now normalized by chain and residue identity, and orientation-specific contact/hotspot selection is covered by tests. The focused verification suite increased from 49 to 51 passing tests. Python compilation and Slurm shell syntax checks passed. `git diff --check` still reports unrelated pre-existing whitespace in `src/rosetta_refinement.py` and `src/surface_extract.py`; those files were not changed in this continuation. Gate 5 remains partial because normalization improves correctness but cannot establish semantic equivalence between the modern and derived asset populations.

## Context and Orientation

The primary launcher is `benchmark/scripts/submit_comparison_batches.sbatch`. It runs `prism.py`, which calls alignment, transformation, and either external Rosetta or PyRosetta. Current transformations are paired files named `{template}_{left}_{right}_{orientation}_L.pdb` and `_R.pdb`; refined outputs live under `processed/rosetta_refinement/` or `processed/pyrosetta_refinement/structures/`.

`benchmark/scripts/stage_current_models_for_main_benchmark.py` discovers refined outputs and creates symlinks plus a staging manifest. Its current filename parser accepts external-Rosetta grammar but rejects retained PyRosetta names. `benchmark/scripts/score_bijective_benchmark_models.py` already performs complete sequence-based chain assignment, requests cross-partner DockQ components, and calculates forward/reverse grouped iRMSD, but it needs stronger hash/version/validation records and fixture coverage.

The source/provenance authorities are:

- `benchmark/scripts/investigation_provenance.py` for hashes, effective configuration, executables, versions, and template assets.
- `benchmark/scripts/investigation_contracts.py` for frozen chain mappings, foreign keys, raw DockQ JSON verification, and native-independent ranking.
- `benchmark/scripts/standardized_evaluator.py` for raw PDB chain/residue integrity, DockQ normalization, and the fail-closed `--no_align` gate.
- `tmp/agent/20260713-investigation-implementation/source-gate-aggregate-final/source_gate_policy.json` for the 240 strict / 17 audit-only source policy.

Published-filter reference behavior is available in Markdown only: `references/nprot.2011.367.md` requires matching thresholds, hotspots, clash rejection, and at least five complementary contacts; the legacy executable reference is `working_version/Multiprot-new/prism-fiberdock-cli/run_files/transformationFiltering.py`.

## Plan of Work

Implement each gate as a small testable unit. First freeze the evidence and statuses so later runs cannot be mislabeled. Next restrict staging to complete refined poses and harden the canonical scorer. Freeze source/template exposure before touching scientific filters. Restore protocol filters as pure functions with parity fixtures, keeping the current geometry-only behavior as an explicitly experimental mode. Only then test backend filtering and concurrency. Run three canaries, followed by the existing preregistered 12-row matched pilot. The full 240-row comparison is the final gate, not a debugging vehicle.

### Task 1: Freeze baseline claims and invalid-run exclusions

**Files:**

- Create: `benchmark/configs/prism_pipeline_verification.json`
- Create: `benchmark/scripts/build_pipeline_verification_baseline.py`
- Create: `tests/test_pipeline_verification_baseline.py`
- Produce: `tmp/agent/20260718-prism-verification/baseline/claims.tsv`
- Produce: `tmp/agent/20260718-prism-verification/baseline/artifact_manifest.tsv`
- Produce: `tmp/agent/20260718-prism-verification/baseline/superseded_runs.tsv`

**Interfaces:**

- Produces `claims.tsv` columns: `claim_id`, `claim`, `classification`, `scope`, `evidence_paths`, `reason`, `superseded_by`.
- Produces a SHA-256 manifest for every retained status, log, script, model, and score file used by later gates.
- Classifications are exactly `supported`, `unsupported`, `exploratory`, or `unresolved`.

- [ ] Write a failing fixture test proving that cancelled v3 evidence is classified `unsupported` even when `exit.json` says completed.

```python
def test_cancelled_scheduler_evidence_overrides_completed_json(tmp_path):
    run = make_run(tmp_path, exit_status="completed", slurm_stderr="JOB 42 CANCELLED")
    row = classify_retained_run(run)
    assert row["classification"] == "unsupported"
    assert row["reason"] == "scheduler_cancelled_status_contradiction"
```

- [ ] Implement deterministic classification and SHA-256 collection. Reject missing paths rather than silently omitting them.
- [ ] Run from the repository root:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q tests/test_pipeline_verification_baseline.py
/home/rshadi25/.conda/envs/gtalign_env/bin/python benchmark/scripts/build_pipeline_verification_baseline.py \
  --config benchmark/configs/prism_pipeline_verification.json \
  --output tmp/agent/20260718-prism-verification/baseline
```

Expected: the v3 mean/high-quality claims are `unsupported`; the 946-run model counts are `supported`; the earlier 18-model scores are `exploratory`; no source artifact changes.

### Task 2: Implement cancellation-aware stage and completion contracts

**Files:**

- Create: `benchmark/scripts/pipeline_completion_contract.py`
- Modify: `benchmark/scripts/submit_comparison_batches.sbatch`
- Modify: `prism.py`
- Create: `tests/test_pipeline_completion_contract.py`
- Create: `tests/test_submit_comparison_signal_status.py`

**Interfaces:**

- `prism.py` appends JSONL records to `PRISM_STAGE_STATUS_PATH` with fields `stage`, `event`, `timestamp`, `return_code`, and `detail`.
- Required current-pipeline stages are `input`, `alignment`, `transformation`, and `refinement`.
- `classify_run(run_root: Path, process_return_code: int, termination_signal: str | None) -> dict[str, object]` returns `process_status`, `scientific_status`, stage completion, paired-transform count, refined-model count, valid-model count, and reasons.
- Allowed scientific states are `completed`, `completed_no_predictions`, `cancelled`, `failed`, and `incomplete`.

- [ ] Write table-driven failing tests for normal completion, completed/no-prediction, missing refinement terminal event, nonzero exit, `SIGTERM`, unpaired halves, zero-byte files, and invalid refined PDBs.
- [ ] Add a signal-safe launcher state variable and traps:

```bash
termination_signal=""
on_term() { termination_signal="SIGTERM"; exit 143; }
on_int()  { termination_signal="SIGINT";  exit 130; }
trap on_term TERM
trap on_int INT
trap finish_status EXIT
```

- [ ] Make `finish_status` call the Python classifier and atomically replace `status/exit.json`; do not infer success with `find ... transformation | wc -l`.
- [ ] Write `status/pipeline_returned.json` only after `prism.py` returns normally. Stage terminal events and this marker must agree.
- [ ] Run:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_pipeline_completion_contract.py \
  tests/test_submit_comparison_signal_status.py
bash -n benchmark/scripts/submit_comparison_batches.sbatch
```

Expected: simulated v3 artifacts classify `cancelled` or `incomplete`, never `completed`.

### Task 3: Enforce assembled refined-pose integrity

**Files:**

- Modify: `benchmark/scripts/stage_current_models_for_main_benchmark.py`
- Reuse: `benchmark/scripts/standardized_evaluator.py`
- Create: `tests/test_stage_current_models.py`

**Interfaces:**

- Extend the parser to return a typed record containing `template_1`, `template_2`, `target_left`, `target_right`, `orientation`, `refinement_backend`, and source model path for both external-Rosetta and PyRosetta grammars.
- Staging records add `source_model_sha256`, observed chain order, receptor/ligand groups, CA counts, integrity status, and integrity reason.
- Discovery accepts only `processed/rosetta_refinement/` and `processed/pyrosetta_refinement/structures/`. `_L.pdb` and `_R.pdb` inputs are rejected as `transformation_intermediate`.

- [ ] Write failing tests using minimal synthetic PDBs for an `_L` half, overlapping partner groups, one-chain refined output, duplicated residue numbering, valid two-partner output, external filename, and PyRosetta filename.
- [ ] Call `validate_raw_pdb_chain_contract()` before creating a symlink. Preserve one manifest row for every discovered candidate, including failures.
- [ ] Run:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_stage_current_models.py \
  tests/test_model_output_integrity.py
```

Expected: the v3 run yields zero stageable models; all 18 retained PyRosetta names are parseable or receive a precise non-schema failure; a valid assembled fixture stages once with disjoint partner chains.

### Task 4: Harden bijective scoring and provenance

**Files:**

- Modify: `benchmark/scripts/score_bijective_benchmark_models.py`
- Reuse: `benchmark/scripts/investigation_contracts.py`
- Reuse: `benchmark/scripts/investigation_provenance.py`
- Reuse: `benchmark/scripts/standardized_evaluator.py`
- Create: `tests/test_score_bijective_benchmark_models.py`

**Interfaces:**

- Per-model output must include source/model/native hashes, benchmark row ID, explicit complete mapping, assignment diagnostics, mapping-validation status, DockQ version and argv, raw JSON path/hash, `GlobalDockQ`, every requested cross-interface row, grouped forward/reverse/min iRMSD, and failure reason.
- Emit `scores_models.tsv` with one row per candidate and `scores_interfaces.tsv` with one row per global/component result.
- Default scoring mode is aligned DockQ. Optional `--no-align` calls `validate_pdb_mapping()` and `no_align_is_safe()` first.

- [ ] Write failing tests for incomplete mappings, overlapping groups, missing native chains, unsafe no-align, absent cross interfaces, corrupt raw JSON hash, and a valid multichain mapping.
- [ ] Replace direct JSON field extraction with `evaluate_dockq_with_frozen_mapping()` and `standardize_dockq_json()`; retain exceptions as `score_failed` rows.
- [ ] Add a negative-control assertion that the v3 `_L` model cannot enter the score manifest.
- [ ] Run:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_score_bijective_benchmark_models.py \
  tests/test_investigation_contracts.py \
  tests/test_standardized_evaluator.py
```

Expected: every staged model has exactly one terminal model row; scored rows have a verifiable raw JSON hash and complete mapping; invalid rows have null structural metrics.

### Task 5: Freeze template provenance and exclude leakage

**Files:**

- Create: `benchmark/scripts/build_template_source_gate.py`
- Create: `tests/test_template_source_gate.py`
- Produce: `tmp/agent/20260718-prism-verification/template_gate/template_assets.tsv`
- Produce: `tmp/agent/20260718-prism-verification/template_gate/target_template_similarity.tsv`
- Produce: `tmp/agent/20260718-prism-verification/template_gate/eligible_templates_by_pair.tsv`
- Produce: `tmp/agent/20260718-prism-verification/template_gate/source_gate_summary.json`

**Interfaces:**

- Inputs are the frozen benchmark manifest, source-gate policy, template list, template asset root, and staged target structures.
- Sequence rows record identities, aligned residues, both coverages, assignment orientation, exact PDB/chain identity, and exclusion reason.
- A hard self-hit is same source PDB/chain or 100% sequence identity with at least 95% coverage of both chains.
- Primary homology exclusion is identity greater than 50% with at least 70% shorter-chain coverage for either mapped partner; sensitivity membership is also written at 30%, 40%, 50%, 70%, and 100%.

- [ ] Write tests for exact self-hit, swapped partner assignment, close homolog, low-coverage fragment, unrelated template, missing sequence, duplicate template, and deterministic output order.
- [ ] Reuse `read_template_ids()`, `preflight_template_assets()`, and provenance hashing. Do not alter `new_template/template/final_list.txt`.
- [ ] Fail a benchmark row closed if either target sequence or required template asset is unresolved.
- [ ] Run the full 19,855-template similarity matrix only through Slurm, writing a unique run root and recording `MaxRSS` and elapsed time.

Expected: every eligible pair has a hash-bound template list; exact/self-homologous templates are absent; all four pipeline arms receive identical per-pair template exposure.

### Task 6: Restore and verify published hotspot/contact filtering

**Files:**

- Create: `src/template_filtering.py`
- Create: `benchmark/scripts/build_protocol_filter_assets.py`
- Modify: `src/transformation.py`
- Create: `tests/test_template_filtering.py`
- Extend: `tests/test_transformation_thresholds.py`
- Extend: `tests/test_transformation_audit.py`

**Interfaces:**

- `load_filter_assets(template_id, root) -> TemplateFilterAssets` returns chain hotspot records, complementary contact pairs, asset hashes, and validation status.
- `evaluate_hotspots(match_dict, hotspots, criterion=2, minimum=1) -> FilterDecision` reproduces legacy same-residue hotspot criterion 2.
- `count_matched_complementary_contacts(left_match, right_match, contacts) -> int` reproduces legacy `interfaceMatchCheck` semantics.
- `evaluate_protocol_candidate(...)` requires normal alignment thresholds, at least one matching hotspot per partner, at least five matched complementary contacts, and the existing clash gate.
- `PRISM_FILTER_MODE` is exactly `published_protocol` or `geometry_only_experimental`; confirmatory launchers require `published_protocol`.

- [ ] First inventory contact/hotspot coverage for all 19,855 templates. Emit `asset_coverage.tsv`; do not infer missing assets from interface PDBs without a recorded derivation method.
- [ ] Parse legacy `template_old/template/contact/*.txt` and `hotspot/hotspot*` into a derived run-local JSON tree. Compare every overlapping current JSON asset against the legacy parse and stop on semantic disagreement.
- [ ] Write legacy-parity fixtures for hotspot pass/fail, residue-type mismatch, exactly 4 versus 5 contacts, missing assets, swapped orientation, and templates with no validated hotspot records.
- [ ] Replace the unconditional `hotspot_analysis(): return True` with the pure filter module. Candidate audit rows must record each filter count and terminal reason.
- [ ] Run:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_template_filtering.py \
  tests/test_transformation_thresholds.py \
  tests/test_transformation_audit.py \
  tests/test_transformation_clashes.py
```

Expected: `published_protocol` fails closed for missing/invalid assets and matches the legacy fixture decisions; `geometry_only_experimental` is visibly labeled in parameters and cannot be collected into confirmatory results.

### Task 7: Verify GTAlign filtering and bound TM-align concurrency

**Files:**

- Modify: `src/alignment_gtalign.py`
- Modify: `src/alignment.py`
- Extend: `tests/test_runtime_and_gtalign_contracts.py`
- Create: `tests/test_alignment_backend_filtering.py`
- Create: `tests/test_alignment_concurrency.py`

**Interfaces:**

- Factor GTAlign parsing into a pure `parse_gtalign_hits(raw_output, query_ids, template_ids) -> list[AlignmentRecord]` function.
- Write JSON only for actual parsed hits; filtered or missing pairs appear in an alignment summary TSV with explicit status, not as fabricated `no_hit` JSON records.
- Replace eager creation of every TM-align future with `iter_bounded_results(tasks, worker, workers, max_pending)` where `max_pending` defaults to `2 * workers`.
- Record aligner version, exact command, raw-output hash, requested/parsed/filtered/error counts, worker count, elapsed time, and peak RSS.

- [ ] Write GTAlign fixtures containing a passing hit, below-prescore hit, below-postfilter hit, malformed hit, duplicate hit, and absent hit. Assert no all-pairs fallback JSON is written.
- [ ] Capture the GTAlign help/version text and test that the production `--pre-score` argument matches the installed 0.19.00 contract.
- [ ] Write concurrency tests with an instrumented fake worker. Assert maximum active workers is bounded, maximum pending work is bounded, each task runs once, failures remain explicit, and one-worker/eight-worker output records are identical after sorting.
- [x] Run the actual 25-template and frozen 946-template alignment-only Slurm contract smokes on both AI and COSBI; retain scheduler/resource manifests and backend output summaries. The one-pair/19,855-template workload remains gated.
- [ ] Run:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_runtime_and_gtalign_contracts.py \
  tests/test_alignment_backend_filtering.py \
  tests/test_alignment_concurrency.py
```

- [x] Through Slurm, run successively: a 25-template contract smoke, the frozen 946-template panel, and one query against exactly 19,855 templates. AI/COSBI raw searched-structure counts reconcile to 50/1,892/39,710 chain interfaces respectively. Record `sacct` elapsed/MaxRSS, raw hits, retained JSONs, transformations, refined models, and failures; only the alignment-stage fields are applicable to these diagnostics.

Expected: counts reconcile at every boundary; TM-align memory is bounded by worker configuration rather than total pair count; GTAlign sensitivity/speed claims remain scoped to the matched panel and recorded hardware.

### Task 8: Run three end-to-end canaries

**Files:**

- Create: `benchmark/scripts/prepare_pipeline_verification_canaries.py`
- Create: `benchmark/scripts/validate_pipeline_verification_canaries.py`
- Create: `tests/test_pipeline_verification_canaries.py`
- Produce: `tmp/agent/20260718-prism-verification/canaries/`

**Canaries:**

1. `cancelled`: launch an agent-owned one-task job, wait until it enters `RUNNING`, then send `SIGTERM`. Acceptance is `process_status=cancelled`, `scientific_status=cancelled`, no score rows.
2. `no_prediction`: use a valid bounded input/template combination known to produce no candidate. Acceptance is normal stage completion and `completed_no_predictions`.
3. `positive_full_pose`: use a previously verified candidate-producing pair/template under `published_protocol`. Acceptance is one or more valid refined full poses and canonical bijective score rows.

- [ ] Before submission, run on the login node:

```bash
hostname
printf 'SLURM_JOB_ID=%s\n' "${SLURM_JOB_ID:-<not-inside-slurm>}"
squeue -u "$USER"
sinfo -o '%P %a %l %D %G'
```

- [ ] Submit with unique roots. Use general `ai`, account/QOS `ai`, and one GPU only for GTAlign GPU; use the minimum live suitable CPU placement for TM-align/Rosetta after checking availability.
- [ ] Validate all canary manifests and hashes before continuing.

Expected: the three status classes are mutually distinguishable; the positive model is assembled/refined and scoreable; the v3 failure mode cannot recur.

### Task 9: Execute the matched 12-row four-arm pilot

**Files:**

- Reuse: `tmp/agent/20260715-matched-benchmark-pilot/analysis_plan.json`
- Reuse/modify only if tests demand: `benchmark/scripts/submit_variant_matrix.sh`
- Modify: `benchmark/scripts/collect_matched_benchmark.py`
- Create: `tests/test_collect_pipeline_verification.py`
- Produce: `tmp/agent/20260718-prism-verification/pilot/`

**Arms:**

- TM-align + external Rosetta baseline.
- GTAlign + external Rosetta, isolating aligner effect.
- TM-align + seeded PyRosetta, isolating refiner effect.
- GTAlign + seeded PyRosetta, completing the interaction cell.

- [ ] Freeze one 12-row manifest with four rigid, four medium, four difficult cases and balanced single-/multichain context.
- [ ] Apply the same 50%-primary eligible template list to all arms per row; verify list hashes before execution.
- [ ] Keep candidate generation, filtering, refinement, and scoring counts by pair and arm. Never replace a failure/null with zero quality.
- [ ] Select candidates using only preregistered native-independent fields. Score native metrics after selection.
- [ ] Report paired differences with denominators and bootstrap confidence intervals, but label the 12-row pilot as calibration—not a general quality conclusion.

Expected go/no-go criteria:

- All task statuses reconcile with Slurm and stage events.
- Every scored model passes full-pose integrity and bijective mapping.
- Every pair/arm has identical eligible-template exposure.
- Raw DockQ JSON hashes and grouped iRMSD are complete.
- Any backend-specific failure is understood and represented before expansion.

### Task 10: Execute the 240-row confirmatory run and publish the evidence

**Files:**

- Create: `benchmark/scripts/validate_confirmatory_prism_run.py`
- Create: `tests/test_validate_confirmatory_prism_run.py`
- Create: `docs/validation/prism-pipeline-verification-20260718.md`
- Create: `docs/validation/prism-pipeline-claims-20260718.tsv`
- Update after verified completion: `.agents/skills/project-memory/references/summary.md`
- Update after verified completion: `.agents/skills/project-memory/references/decisions.md`
- Update unresolved items: `.agents/skills/project-memory/references/open_questions.md`

- [ ] Open the 240 strict rows only after the pilot report passes all gates. Keep the 17 audit-only rows out unless the frozen source policy is authoritatively revised.
- [ ] Use new output roots for every arm and immutable per-task retry IDs. Retry scheduler/infrastructure failures only; never overwrite scientific failures.
- [ ] Aggregate coverage, prediction production, scoreability, and quality separately for all rows, single-chain rows, and multichain rows.
- [ ] Require the final validator to prove foreign-key completeness, unique model identity, template-list hash equality, no self/homology leakage, stage/status reconciliation, full-pose integrity, mapping completeness, raw JSON hash integrity, and denominator accounting.
- [ ] Run focused checks first, then the stable suite:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_pipeline_completion_contract.py \
  tests/test_stage_current_models.py \
  tests/test_score_bijective_benchmark_models.py \
  tests/test_template_source_gate.py \
  tests/test_template_filtering.py \
  tests/test_alignment_backend_filtering.py \
  tests/test_alignment_concurrency.py \
  tests/test_validate_confirmatory_prism_run.py
bash benchmark/scripts/run_stable_checks.sh
```

- [ ] Publish claims only from validator-approved rows. The final claims table must state effect estimate, denominator, uncertainty, scope, confounders, and evidence paths.

Expected: the report can say which arms executed, predicted, produced valid models, and achieved which quality distribution without conflating those stages. If an arm fails a gate, the report states `blocked` or `unresolved`; it does not infer performance.

## Concrete Steps

All commands below use `/scratch/rshadi25/GitHub/PRISM-prescript` as the working directory.

1. Create a unique run root without modifying retained evidence:

   ```bash
   run_id=20260718-prism-verification
   mkdir -p "tmp/agent/$run_id"
   ```

2. Implement Tasks 1–7 in order, following the failing-test, minimal-implementation, passing-test cycle in each task. After each task, update this plan's `Progress`, `Surprises & Discoveries`, and `Decision Log`.

3. Before any Slurm work, run the focused unit suite and syntax checks. Expected result is all tests passing and `bash -n` returning zero.

4. Run Task 8 canaries. Stop if scheduler status, stage status, and scientific status disagree or if the positive model cannot be scored canonically.

5. Run Task 9's 12-row pilot. Stop if template exposure differs across arms, any score lacks provenance, or failure accounting is incomplete.

6. Review the pilot artifact manifest and report. Only a recorded go decision permits Task 10.

7. Run the 240-row strict cohort, validate it independently, publish the report/claim ledger, and update project memory with only durable verified outcomes.

Frequent commits are recommended at the end of each accepted gate, but do not commit unrelated pre-existing changes. Suggested messages are `test: define pipeline completion contract`, `fix: reject partial pipeline completion`, `fix: require full refined poses for scoring`, `feat: freeze benchmark template source gate`, `fix: restore protocol candidate filters`, `perf: bound alignment task concurrency`, and `docs: publish PRISM verification evidence`.

## Validation and Acceptance

The complete plan passes only if all of the following are observable:

- A cancelled canary is classified cancelled even if partial transformation files exist.
- A normal zero-prediction run is distinguished from cancellation and failure.
- No `_L.pdb` or `_R.pdb` transformation half is stageable or scoreable.
- Every scored model contains both disjoint partners, has valid residue numbering, and maps bijectively to the native complex.
- Every score retains model/native hashes, exact command and version, complete mapping, raw DockQ JSON/hash, cross-interface values, `GlobalDockQ`, and grouped iRMSD.
- `--no_align` is impossible without strict identity/numbering validation.
- Every benchmark row has an immutable, per-pair eligible template list; exact self-hits and primary >50% homologs are excluded and sensitivity memberships are reported.
- Confirmatory candidates use validated hotspot/contact assets, at least one matching hotspot per partner, at least five matched complementary contacts, and the clash gate.
- GTAlign output counts reconcile from raw hits through filters; no all-pairs placeholder explosion occurs.
- TM-align concurrency is bounded and deterministic across worker counts, with peak memory recorded at scale.
- The 12-row pilot uses identical rows, templates, filters, refiners, scoring, and failure accounting across the intended paired comparisons.
- The 240-row run excludes the unresolved 17 source-gate rows and reports all denominators explicitly.
- The v3 means/high-quality counts remain marked unsupported and are absent from new aggregates.

No single DockQ value, successful Slurm exit, or model count is sufficient for acceptance.

## Idempotence and Recovery

The first execution uses `tmp/agent/20260718-prism-verification/` and refuses nonempty destinations. Baseline and template manifests are content-hashed; reruns compare hashes before reusing inputs. Task outputs are append-free and written atomically through temporary files followed by rename. Scheduler retries use `tmp/agent/20260718-prism-verification/retry-N/` and never overwrite a scientific attempt.

If a gate fails, retain its manifest, logs, status, and failure rows, update this plan, and rerun only that gate in a new retry root. Do not cancel unrelated jobs. The intentional-cancellation canary may be cancelled only after confirming the exact agent-owned Slurm ID. Derived protocol-filter assets can be regenerated from their immutable sources; `new_template/template/`, `template_old/template/`, benchmark originals, and v3 artifacts remain untouched.

Temporary unit-test files remain under pytest-managed directories. Any task-created disposable `/tmp` files are removed at the end of that task. No cleanup of pre-existing repository artifacts is part of this plan.

## Artifacts and Notes

Primary negative controls:

- `tmp/agent/20260718-benchmark20k-v3/`
- `tmp/agent/20260717-benchmark55/scoring/score_final.py`

Reusable evidence/contracts:

- `docs/exec-plans/20260717-dockq-irmsd-scoring-contract-audit.md`
- `docs/exec-plans/20260716-score-contract-and-cleanup.md`
- `benchmark/scripts/investigation_provenance.py`
- `benchmark/scripts/investigation_contracts.py`
- `benchmark/scripts/standardized_evaluator.py`
- `benchmark/scripts/score_bijective_benchmark_models.py`
- `tmp/agent/20260713-investigation-implementation/source-gate-aggregate-final/source_gate_policy.json`

Reference Markdown, without reading PDFs:

- `references/nprot.2011.367.md`
- `references/Proteins - 2011 - Tuncbag - Fast and accurate modeling of protein protein interactions by combining.md`

Keep short command transcripts and summary TSV/JSON artifacts in the run root. Do not copy giant alignment logs into this plan.

## Interfaces and Dependencies

- Pipeline interpreter: `/home/rshadi25/.conda/envs/gtalign_env/bin/python`.
- DockQ interpreter: `/scratch/tmp/prism-dockq-env/bin/python`.
- GTAlign executables: `/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu` and `gtalign_gpu`; capture `--version`/help and binary hash per run.
- TM-align executable: resolve from the effective runtime manifest and record its hash.
- External Rosetta: load `rosetta/2022.42` inside Slurm and record prepack/dock/database paths.
- PyRosetta: opt-in adapter with frozen `-mute all -constant_seed -jran 12345`; no fallback to external Rosetta.
- Scheduler: verify live partition/account/QOS/GPU compatibility before each submission. Current known GPU placement is partition/account/QOS `ai/ai/ai` with `--gres=gpu:1`; do not assume that availability remains unchanged.
- Native data: frozen Benchmark 5.5 role files and the source-gate policy. No derived native complex may be written into `benchmark/originals/`; any needed assembly belongs in the run root with source hashes.
- Template assets: `new_template/template/` and legacy `template_old/template/` are read-only inputs with different schemas. Derived adapters must record source and output hashes and may not claim equivalence without semantic comparison.
