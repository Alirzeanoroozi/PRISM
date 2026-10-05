# Historical PRISM vs Current Pipelines Investigation

## Objective

Make the historical `benchmark/prism_processed_results/benchmark_reports/cross_benchmark_report.md` reproducible as an evaluation target, while preventing unequal template assets, stage attrition, mapping drift, and multichain score aggregation from being interpreted as causal pipeline effects.

This plan deliberately starts with provenance, template compatibility, evaluator contracts, and lineage. It does not authorize a full 257-pair rerun until those contracts pass diagnostic validation.

## Repository state and constraints

- `PRISM-prescript` is the canonical benchmark/scoring repository.
- The worktree is dirty with unrelated current-pipeline and benchmark artifacts; existing changes are preserved.
- New run-specific outputs belong under `tmp/agent/<timestamp>-historical-current-investigation/`.
- `benchmark/README.md` and `benchmark/prism_processed/README.md` remain workflow references.
- Single-chain and multichain results must remain distinguishable in every summary.
- The two unavailable native inputs (`1erk`, `4zai`) remain explicit failures, not silently removed rows.

## Evidence currently accepted as input, not causal proof

- The historical/working-version template manifest has 21,072 entries; the current run used 946.
- The July legacy arm listed many templates but had complete staged assets for only one template.
- Current and working TM-align binaries are identical (20220412), while the current invocation/parser and runtime assets differ.
- Current comparison outputs are unpaired: 143 current structures across 56 pairs versus two legacy structures for one disjoint pair.
- Current multichain scoring previously retained a final component DockQ value rather than an explicitly stored GlobalDockQ.

These observations motivate the controls but do not establish that any one aligner or refiner is intrinsically better.

## Artifact contracts

Each investigation run must contain:

1. `provenance.json`: repository state, code/config/tool/database/template hashes, resolved executables, package versions, redacted environment exports, seeds, Slurm resources, and exact argv commands.
2. `template_assets.tsv`: one row per manifest entry, including duplicate/unique status, required assets, format, per-asset hashes, and usability classification.
3. `lineage.tsv`: one row per attempted pair/template/orientation lineage with stage statuses and one terminal status/failure reason.
4. `poses.tsv` and immutable PDB files: pose ID, coordinate hash, pair/template/orientation, chain/residue provenance, and native-independent ranking inputs.
5. `scores_global.tsv` and `scores_interfaces.tsv`: raw DockQ JSON hash/path, mapping declaration, GlobalDockQ, and one row per component interface.
6. `pair_summary.tsv`: all-pair unconditional utility (zero when no eligible scoreable model), conditional quality, coverage, model counts, interface measurements, and efficiency.
7. A report that separates observations, controlled estimates, unresolved hypotheses, and stage-by-stage before/after changes.

## Implementation sequence

### Phase 1 — provenance and template preflight

- Build deterministic stdlib helpers for SHA-256 hashing, safe environment capture, package/executable resolution, command capture, and effective configuration capture.
- Preflight every listed template before execution and report listed, unique, valid, fully resolvable, and missing counts.
- Support both legacy and modern asset roots without assuming that a filename alone proves compatibility.
- Fail closed for a requested execution arm when any required template asset is missing or format-incompatible.

### Phase 2 — standardized evaluation

- Parse DockQ 2.1.3 JSON while retaining raw JSON and its hash.
- Store GlobalDockQ separately from component interfaces; never overwrite repeated interface records.
- Preserve DockQ iRMSD/LRMSD/fnat/F1/clashes and custom grouped iRMSD as distinct fields.
- Freeze benchmark-role chain mapping before scoring. Use symmetric-chain equivalence only when declared in advance.
- Allow `--no_align` only after residue identity and numbering correspondence validation.

### Phase 3 — lineage and common pose contract

- Make stage transitions explicit: candidate generation → alignment → filter → pose → refinement → energy gate → scoring.
- Assign every attempted lineage exactly one terminal status, including missing input, filter rejection, refinement failure, score failure, and no-scoreable-model.
- Copy accepted pose PDBs into an immutable run directory and record coordinate hashes.
- Enforce a native-independent ranking key allowlist for equal-budget selection.

### Phase 4 — controlled diagnostics

- Run one-template fixtures and the equal 946-template arm only after preflight passes.
- Compare TM-align and MultiProt raw candidate tables on identical query/template/orientation inputs.
- Replay native filters and one common canonical filter without using native scores for selection.
- Send byte-identical accepted poses through no-refinement, Rosetta, and FiberDock crossover arms.
- Exercise adversarial multichain fixtures: second-chain-only match, reversed chain order, overlapping IDs, split partners, and insertion/gap numbering.

### Phase 5 — confirmatory benchmark

- Freeze code, templates, mappings, scoring, analysis, and seeds before running the 255 input-available pairs.
- Primary outcome: all-pair best GlobalDockQ@20, with zero for eligible pairs lacking a scoreable model.
- Secondary outcomes: success at DockQ ≥0.23/0.49, conditional quality, shared-success comparisons, model/pair/template coverage, failure reasons, wall time, CPU-hours, memory, and call counts.
- Use paired pair-level bootstrap/permutation inference, clustered by native complex or sequence family where applicable.

## Acceptance checks

- Every execution template is fully resolvable before launch.
- Native self-score, rigid-transform invariance, chain permutation, insertion/gap, symmetric-chain, and multichain evaluator tests pass.
- Mapping is frozen before DockQ values are inspected.
- All stage counts are pair-specific and batch completion is never treated as pair success.
- Refiner crossover inputs have identical bytes and coordinate hashes.
- Equal-budget ranking contains no native-derived metric.
- Conditional scores are never used as unconditional superiority claims.
- A factor is called causal only if changing that factor alone changes downstream all-pair utility on the frozen confirmatory set.

## Current risks and explicit non-goals

- This implementation does not itself claim a causal ranking of TM-align, MultiProt, Rosetta, or FiberDock.
- Full historical asset validation may be storage- and I/O-heavy; preflight should be run on the login node, while alignment/refinement workloads belong in Slurm.
- External binary behavior remains an empirical question until the diagnostic panel is run with matched inputs.

## Implementation status — 2026-07-13

Completed infrastructure and validation evidence:

- The repository benchmark cohort is frozen as the existing 257 CSV rows: 162 rigid, 60 medium, and 35 difficult. `dataset_row_id` is the primary key; normalized PDB pairs are audit fields only.
- Curated archive roles are explicit and fail closed: `r_u`/`l_u` are pipeline inputs and `r_b`/`l_b` are native truth. Files under `benchmark/data/pdbs` remain audit-only selectors and are never substituted.
- Biopython validation records parser/version, models, all-chain and polymer-chain IDs, residue-ID and sequence hashes, altloc/disorder counts, duplicate identifiers, warnings, and chain-set status.
- Source-gate tasks were executed as 26 isolated KUTEM arrays (25 arrays of 10 plus one array of 7). Each task has its own inputs, outputs, logs, `exit.json`, Slurm identity, hashes, and scientific status.
- Current root pipeline smoke passed after fixing canonicalized alignment filenames in `src/transformation.py`; the smoke is plumbing evidence only and produced zero accepted pairs.
- The MultiProt compatibility snapshot is hash-pinned and non-destructive. The initial standalone Python 2.7.15/PyMySQL environment was superseded by the project-local NumPy stage and complete derived toolchain gate documented below; the legacy native dependency limitations remain explicit.

Source-gate rerun status:

- The first source-array run recorded 225/257 successful tasks and 32 deliberate source-gate failures. Scheduler/task accounting reconciled all 257 IDs; the 31 archive failures were resolver/chain-contract artifacts in the pre-correction run, and `medium:000055` was a real curated-ligand chain/sequence disagreement (`3CPH_l_u` chain C versus CSV selector `1G16_A`).
- The parser gate was corrected to compare polymer chains while retaining all-chain metadata, and the runner/aggregator now hash all consumed artifacts, bind summaries to task directories and manifests, reject duplicate/missing row identities, and reject stale retry outputs.
- The final corrected 257-task source-array rerun completed under `tmp/agent/20260713-investigation-implementation/source-array-submission-v3/`. Aggregation under `tmp/agent/20260713-investigation-implementation/source-gate-aggregate-final/` reconciled all 257 task identities and showed 257/257 curated role sets present and staged, 240/257 chain/parse-clean rows, and 17 explicit chain-contract failures. No model confirmatory arrays were authorized.

The first source-array outputs remain immutable observational evidence. They are not merged with the corrected rerun and are not used for confirmatory estimates.

## Implementation status — legacy tool environment and controller smoke (2026-07-14)

Completed:

- Installed NumPy 1.16.6 into the project-local Python 2 site directory with
  PyMySQL 0.9.3 already available; the shared `tmalignRosetta` environment was
  not modified. The activation script now resolves its own path correctly when
  sourced from Bash.
- Staged all checked-out MultiProt, NACCESS, POPS, and FiberDock payloads in a
  derived environment with per-file hashes, modes, archive comparisons, Python
  dependency probes, `ldd` records, and explicit NACCESS profile selection.
- The compatibility NACCESS profile uses the repository current binary
  (`f05779ca...`, `libgfortran.so.5`) and is operational. The historical
  working-version profile is retained but blocked by missing `libgfortran.so.3`.
- KUTEM array `1355097` independently passed NACCESS, POPS, MultiProt, and
  FiberDock probes. FiberDock produced `resFile.ref` in an energy-only smoke;
  this does not validate hydrogen/NMA refinement.
- The first controller smoke exposed and fixed a compatibility-adapter defect:
  `TemplateChecker` accepted `work_path` but did not store it. The rebuilt
  snapshot was tested directly and the controller rerun was isolated by profile.
- KUTEM array `1355109` reached preprocessing, surface extraction, alignment,
  transformation filtering, and refinement setup under the compatibility
  profile. It produced zero filtered candidates and is classified explicitly as
  `pipeline_plumbing_complete_no_candidates`, not as a successful docking run.
  The historical profile used a recorded derived wrapper-relocation patch and
  then failed at `accall` because `libgfortran.so.3` is unavailable; missing
  intermediate files are explicitly classified as
  `pipeline_blocked_missing_intermediate`. Both results have task-local logs,
  commands, parameters, hashes, and `exit.json` records.

Current derived artifacts:

- `tmp/agent/20260713-investigation-implementation/legacy-tool-environment-v5/`
- `tmp/agent/20260713-investigation-implementation/legacy-tool-environment-historical-v4/`
- `tmp/agent/20260713-investigation-implementation/legacy-tool-probes-v5/`
- `tmp/agent/20260713-investigation-implementation/legacy-pipeline-smoke-v9/`
- `tmp/agent/20260713-investigation-implementation/multiprot-compat-v5/`

Remaining blockers:

- No candidate reached FiberDock from the selected synthetic smoke pair/template;
  a matched positive template/pair is required to validate the controller's
  full hydrogen/NMA/FiberDock path.
- The historical NACCESS and bundled 32-bit helper binaries remain unavailable
  on the host. The compatibility profile must remain analytically separate from
  the historical arm.
- The exact BM3.0 88-case source list and the 17 source-gate orientation
  decisions remain unresolved; no confirmatory benchmark arrays are authorized.

## Validation update — 2026-07-14

The remaining runtime and controller boundaries were exercised with isolated KUTEM arrays under the fixed `array-kutem` profile.

- Historical NACCESS is operational in `legacy-tool-environment-historical-v5/` by exposing the cluster-provided GCC 6
  `libgfortran.so.3` directory through task-local `activate.sh`. Historical tool array `1355128` passed NACCESS, POPS, MultiProt,
  and the energy-only FiberDock probe. The probe does not establish full refinement readiness.
- The positive plumbing fixture is benchmark row `T_Rigid.csv:47`, selectors `1RGH_B` and `1A19_B`, template `1b27AD` with chains
  A/D. At the reference 50% filter, MultiProt's best interface-A match is 22/45 (48.9%), so both adapters produce zero candidates.
  The existing TMalign artifact accepts the same pair/template, establishing an alignment-stage difference but not a causal
  production estimate.
- A recorded one-factor diagnostic at 40% produced one candidate in both current and historical controller tasks. KUTEM array
  `1355189` reached FiberDock intermediate files and was classified `pipeline_blocked_full_refinement_capability`; the override is
  diagnostic-only and excluded from benchmark conclusions.
- The staged environment records the full-refinement blockers as 32-bit `nma`, `reduce.2`, and `reduce.3` binaries. No host
  mutation, binary modification, or library substitution was performed.
- `benchmark/scripts/freeze_source_gate_policy.py` generated
  `tmp/agent/20260713-investigation-implementation/source-gate-aggregate-final/source_gate_policy.json` (SHA256
  `2d652ad7ba4ebd4b7dfbad168458af85cbb309edf177d46c2a12a8cd5a5767db`). It records 240 strict rows, 17 audit-only rows,
  `blocked_source_authority`, and prohibits automatic orientation swaps, full-PDB substitution, or repair.

These results close the historical NACCESS blocker and the missing-positive-fixture blocker for plumbing. They do not authorize
the 257-row confirmatory run: the exact BM3.0 source list remains unavailable, 17 BM5/5.5 rows remain unresolved, and complete
FiberDock refinement is not executable with the staged 32-bit helper payloads.

## Implementation update — 2026-07-14

Implemented additive validation contracts:

- Added `environment.yaml` and `runtime_manifest.json`; the existing `environment.yml` remains preserved as a compatibility recipe.
- Added `benchmark/scripts/validate_runtime_manifest.py` and focused tests for binary hashes, environment contracts, and GTalign provenance.
- Changed root TM-align execution to use a unique per-call scratch directory, explicit subprocess return-code checking, and missing-output rejection. Shared `matrix.out`/`out.tm` files are no longer used by the current adapter.
- Changed GTalign runs to use unique run directories and reject nonempty stale output directories. Each parsed record retains query- and reference-normalized TM-scores plus the raw output hash.
- Corrected legacy capability classification so 32-bit architecture alone is not treated as execution failure. `nma` is recorded as observed-working; full FiberDock remains fail-closed until an end-to-end refined-output probe passes.
- Added a non-destructive `benchmark/scripts/build_cleanup_manifest.py` and generated `tmp/agent/20260714-pipeline-validation/cleanup_manifest.tsv`.
- Added `docs/pipeline-validation-report-20260714.md` and the NotebookLM research handoff at `research/refinement-tools-nlm-handoff-20260714.md`.

Focused validation passed 14 tests in `gtalign_env`, Python compilation, runtime-manifest validation, and Slurm-script syntax. Full
benchmark arrays and the exact `reduce.2` runtime experiment remain gated on smoke validation and live Slurm/network setup.

## Implementation update — option-1 recovery, option-2 probe, and benchmark replay preparation (2026-07-14)

Progress:

- [x] Recover and hash the missing historical model from the archived joblist.
- [x] Re-score the recovered original model with the fail-closed raw-PDB contract.
- [x] Stage a distinct-chain FiberDock probe with separate A/B inputs and isolated KUTEM tasks.
- [ ] Complete the patched distinct-chain FiberDock probe and classify its terminal boundary.
- [ ] Run the hardened observational benchmark replay and compare it with the previous report.
- [ ] Validate optional PyRosetta adapter availability without changing the default Rosetta CLI backend.

Evidence and decisions:

- The historical original model is recoverable from `benchmark/joblist_001-045_20260223.tar.gz` and has SHA-256
  `44eef786f101e232f293afafad133bd39519a2ee9f35fe941a3c150908ab26b0`; its raw ATOM records contain only chain `B` and
  96 residues. The recovered model is staged under `tmp/agent/20260714-historical-artifact-recovery/recovered/` and the
  archived chain-fixed counterpart remains absent.
- The fail-closed scorer rejects the recovered original when the declared model roles are `B:B`, with explicit errors for
  overlapping model roles, fewer than two raw chains, and a residue-number reset. The previous report's `DockQ=0.854188`/
  `iRMSD=0.612` therefore depends on the missing chain-fixed artifact and is not reproducible from the original model bytes alone.
- Distinct-chain FiberDock tasks `1355834` and `1355835`/`1355842` demonstrated, respectively, a launcher-parent failure and
  a compatibility TemplateChecker path failure. Task `1355856` is the first run after the nested workspace-manifest lookup
  correction; its output must be inspected before any claim about FiberDock chain preservation.
- The final batch will be labeled `observational_replay` unless all source-gate, evaluator, and paired-score requirements pass.
  It must preserve unavailable inputs and non-scoreable models as explicit rows and may not use the previous report's aggregate
  means as a causal method comparison.

## Implementation update — final benchmark replay and status (2026-07-14)

Progress:

- [x] Complete the patched distinct-chain FiberDock probe and classify its terminal boundary.
- [x] Run the hardened observational benchmark replay and compare it with the previous report.
- [x] Validate optional PyRosetta adapter availability without changing the default Rosetta CLI backend.

Results:

- Distinct-chain probe KUTEM job `1355856` reached Flexible Refinement but produced zero transformations/candidates and no
  FiberDock output. It is inconclusive for FiberDock chain preservation and is retained as a pre-FiberDock failure boundary.
- Strict replay KUTEM job `1355887` processed 145 retained model rows across ten isolated tasks. All rows failed the explicit
  complete-correspondence/no-align contract and remain null rather than being promoted to zero structural quality.
- Alignment-enabled audit KUTEM job `1355897` processed the same rows. It reproduced the two legacy FiberDock scores to
  rounding and yielded 119 current DockQ values, but cannot be compared directly with the previous current report because the
  evaluator regime differs.
- PyRosetta probe reports `ModuleNotFoundError: No module named 'pyrosetta'`; the optional adapter remains unavailable and the
  external Rosetta module remains the stable default.

The final replay is therefore complete as an observational validity check, not as a causal 257-row benchmark comparison.
Confirmatory generation remains blocked by source-gate disagreements, missing chain-fixed historical artifacts, incomplete
strict mappings, and unverified historical FiberDock full-refinement capability.

The hardened rerun supersedes those initial replay IDs: strict job `1355946` and aligned job `1355947` retained raw DockQ JSON,
verified output hashes and model coverage at collection, and corrected iRMSD best-value direction (`min`, not `max`). Final
summaries are under `tmp/agent/20260714-observational-score-replay-{strict,aligned}-v4/collected/`.
