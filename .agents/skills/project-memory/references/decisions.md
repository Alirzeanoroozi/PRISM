# Verified decisions

## Verified decisions

### Stable execution boundary

Decision:
- `docs/STABLE_PIPELINE.md` is the operational source for the current pipeline.
- Keep the current `prism.py` implementation separate from the legacy Python 2
  tree and from diagnostic/rescue configurations.

Reason:
- Mixing legacy assets, relaxed thresholds, or unpaired benchmark outputs
  invalidates causal comparison.

Consequences:
- Use the standard `gtalign_env` interpreter; load Rosetta explicitly only for
  `--refiner external_rosetta`; run substantial work through Slurm.

### Resource discovery boundary

Decision:
- Treat the runtime-resource table in `summary.md` and
  `docs/STABLE_PIPELINE.md` as the first lookup for interpreters, executables,
  staged assets, output roots, and environment controls.
- Do not infer an executable, asset directory, or Python environment from
  `PATH`, a stale graph node, or a sibling checkout.

Reason:
- Isolated workspaces use symlinks, compute nodes may not inherit interactive
  paths, and the repository includes large legacy/vendor trees with similarly
  named tools.

Consequences:
- Verify a required path before launching. Record the resolved executable,
  interpreter, staged template root, and output root in run provenance.
- Graphify remains navigation-only until rebuilt with path-qualified IDs and a
  deliberate corpus scope. That scope must retain curated `tmp/agent` run
  records (notes, manifests, status files, and logs); it must exclude only
  copied environments, vendored dependencies, caches, and Graphify's own
  generated corpus trees within those directories.

### Candidate ranking boundary

Decision:
- Ranking is opt-in (`--rank` or `PRISM_RANK`) and operates only after
  transformation, before refinement.
- The deterministic baseline may reduce refinement load, but learned ranking
  remains disabled in production.
- Ranked runs default to a fresh audit; explicitly named audit paths remain
  append-only.
- Coverage contributes to the deterministic score only when both template
  partner denominators are known. Otherwise use TM-score evidence alone; do
  not substitute a fixed `match_count / 50` pseudo-coverage.

Reason:
- The available grouped pilot is too small and showed no top-1 native-like
  improvement over the deterministic baseline. Fresh audits prevent prior runs
  from changing a new selection.

Consequences:
- Rank only complete transformed candidates. Failed or unlabeled candidates
  remain in the audit but are not rankable/training negatives.
- A candidate whose canonical score is explicitly `not_scoreable` remains in
  candidate/audit tables for provenance but is excluded from the rankable
  candidate-versus-score identity-set equality. Do not treat its missing label
  as a pipeline loss or invent a negative label.
- Top-K selection is per receptor/ligand pair and rejects nonpositive K.
- Offline candidate tables must rank and evaluate each `dataset_row_id`
  independently. Use `native_complex_id` only for legacy tables without row
  identity, and never report a single global top-1 across unrelated complexes.
- Job 1392725 verifies that deterministic top-1 selection can reduce a
  two-candidate current-pipeline case to one refinement with identical
  pre-ranking provenance. This authorizes a resource-saving experiment only;
  it does not authorize production-by-default ranking or a quality claim.
- Future pipeline audits receive actual template-chain sizes from
  `transformation.template_size`, so real match coverage is available without
  affecting the stable unranked path. Historical reconstructed rows lacking
  denominators remain rankable from TM-score alone.

### Benchmark and scoring boundary

Decision:
- Keep benchmark identity as `dataset_row_id`, preserve raw selectors, and use
  curated archive role files as native truth. Do not substitute a full PDB for
  a chain-qualified source without a separate validated materialization step.
- Report global model scores separately from per-interface artifacts. Do not
  treat invalid multichain `HLA:AC`/`HLA:BC` mappings as canonical metrics.

Reason:
- Normalized PDB pairs collapse distinct benchmark rows; chain-role ambiguity
  changes biological interpretation.

Consequences:
- Unresolved rows stay explicit audit failures and cannot enter the strict
  confirmatory denominator.
- Canonical current-output scoring uses
  `benchmark/scripts/stage_current_models_for_main_benchmark.py`,
  `assemble_bm55_native_complexes.py`, and
  `score_bijective_benchmark_models.py` in that order.
- Select one most-refined external-Rosetta artifact per generated pose. Derive
  model partner-group sizes from query receptor/ligand selectors, never from
  the two template-interface chains.
- Use complete within-partner model:native bijections. Preserve GlobalDockQ,
  but use requested receptor-ligand cross-interface DockQ as the ranking label;
  preserve grouped forward/reverse iRMSD separately.
- If complete-mapping DockQ 2.1.3 fails with the verified empty-array signature,
  recover only requested receptor-ligand interfaces. Set score scope to
  `requested_cross_interfaces_only`, leave GlobalDockQ unavailable, and retain
  every requested pair as `scored` or `no_native_interface` with component
  argv, JSON, hashes, and the original complete-mapping traceback. Do not use
  this fallback for unrelated errors.
- Keep `benchmark/scripts/irmsd.py` unchanged. Use
  `benchmark/scripts/irmsd_grouped_safe.py` for repaired multichain scoring;
  it retains legacy residue/interface/RMSD calculations while replacing the
  mutable symmetry dictionary and incorrect permutation-count construction.
- Preserve the DockQ virtual-environment entry path exactly in batch launchers;
  do not resolve its Python symlink into the base interpreter.
- Attach ranking labels only when `dataset_row_id` and model SHA256 both match.
  Use `dockq_cross_mean` and `irmsd_grouped_min`; exclude non-scored and
  audit-only rows. Path-only joins and GlobalDockQ labels are legacy-only.

### Protocol and completion boundary

Decision:
- Keep `published_protocol` and `geometry_only_experimental` as explicitly
  separate transformation modes.
- Treat a result as scientifically complete only after normal pipeline return,
  terminal stage records, refined-pose integrity, and canonical scorer output
  agree; scheduler cancellation overrides partial outputs.

Reason:
- Missing protocol assets and cancellation can otherwise look like successful
  predictions through surviving directories or transformation halves.

Consequences:
- Do not score `_L.pdb`/`_R.pdb` transformation intermediates as docked models.
- Preserve raw DockQ JSON, model/native hashes, complete model-to-native chain
  mapping, and explicit failures in every benchmark collector.

### Template exposure and comparison boundary

Decision:
- Freeze source inputs, template panel, candidate budget, ranking key, and
  evaluator before comparing aligners or refiners.
- Hard-exclude exact self hits; retain the documented >50% global sequence
  identity / ≥70% shorter-chain-coverage homology exclusion as the primary
  comparison setting, with sensitivity reporting when applicable.

Reason:
- Existing July current-versus-legacy aggregates use unequal staging and are
  observational rather than causal.

Consequences:
- Current/legacy, TM-align/GTalign, and external-Rosetta/PyRosetta quality
  claims require a newly paired run; do not infer them from unrelated counts.

### Legacy compatibility execution boundary

Decision:
- Treat `Bad system call`/exit 159 from legacy 32-bit MultiProt or Reduce as a
  potential seccomp execution-context failure until it is reproduced in an
  approved unrestricted or Slurm context.

Reason:
- The archived `test3-escalated-20260717` validation produced parsed MultiProt
  output outside the restricted sandbox, while direct restricted executions
  failed with exit 159.

Consequences:
- Do not replace legacy binaries or diagnose an intrinsic alignment failure
  from that signature alone. Preserve the execution context and validate raw
  plus parsed alignment output before changing compatibility tooling.

### FiberDock provenance boundary

Decision:
- `prism.py --refiner fiberdock` is an available current-pipeline backend, but
  it is not evidence of historical MultiProt/FiberDock equivalence.

Reason:
- The historical helper stack still requires a permitted compatible 32-bit
  runtime, and current-versus-legacy inputs, staging, and evaluation are not
  yet paired.

Consequences:
- Do not substitute `reduce.3` for the historical `reduce.2` helper or use
  current FiberDock outputs as legacy-equivalent benchmark evidence. Any such
  comparison requires matched provenance and evaluator contracts.

### DockQ evaluation boundary

Decision:
- DockQ evaluation is opt-in (`--compare` or `PRISM_COMPARE`), runs after
  refinement as the final pipeline stage.
- Multi-worker parallel evaluation via `--compare-jobs`.

Reason:
- Requires native PDBs and DockQ dependencies not needed for prediction-only
  runs. Keeps the prediction-only fast path intact.

Consequences:
- The comparison adapter parses output filenames to reconstruct metadata.
  Filename convention (`{template}_{receptor}_{ligand}_{orientation}_{L|R}.pdb`)
  is the contract between `transformation.py` and `compare.py`.
- `compare_pairs_from_outputs()` handles the current 2-tuple passed_pairs
  format; `compare_and_summarize()` handles legacy 4-tuple format.
## Solved questions

- GTalign’s earlier `--pre-score` behavior was double filtering. The prefilter
  now has a separate `PRISM_GTALIGN_PRE_SCORE` control; the remaining sparse
  hit rate on the broad template panel is expected biology/coverage behavior.
- The current pipeline now has a real MultiProt backend; TM-align is not used
  as a hidden fallback.
- FiberDock is integrated as `prism.py --refiner fiberdock`; older notes that
  describe it as unavailable in the current pipeline are superseded.
- PyRosetta is ready in `gtalign_env` and remains an explicit refiner choice,
  not a fallback for external Rosetta.
- Ranking audit path loss and stale-audit candidate loss were corrected on
  2026-07-25 and covered by focused regression tests.
- The repaired deterministic ranking branch completed a paired stable-threshold
  smoke in job 1392725: identical two-candidate audits entered both arms and
  top-1 ranking selected/refined one candidate while baseline refined two.
- The “all 18 pipeline variants verified” statement is not retained as a
  decision: the exact archive confirms nine combinations, while later
  18-variant notes describe incomplete execution.
- The canonical multichain score contract is resolved for strict-clean rows:
  row-specific curated bound-role native assembly, complete bijective mapping,
  raw DockQ 2.1.3 JSON, requested cross-interface metrics, and grouped
  forward/reverse iRMSD. Job 1392974 validates it on 235 final poses.
- External-Rosetta suffixes are sequential refinement artifacts, not three
  independent candidate poses. The double-suffixed artifact is preferred;
  lower suffixes are fallback only when a later artifact is absent.
- A virtual-environment executable path is part of provenance. Job 1392965
  proved that resolving `/scratch/tmp/prism-dockq-env/bin/python` to its base
  interpreter removes DockQ in clean Slurm shells; job 1392974 verified the
  former corrected launcher. The currently verified replacement is
  `benchmark/prism_processed/env/prism_score_env/bin/python` (DockQ 2.1.3),
  because the old `/scratch/tmp` entry now lacks DockQ.
- The retained GTalign JSON is sufficient to reconstruct current-run ranking
  features without rerunning alignment. Batch 1 recovered and hash-joined all
  235 candidates; ranking identity is `dataset_row_id` plus model SHA256 and
  ranking order is local to that row.

### Safe grouped-iRMSD repair contract

Decision:
- Keep `benchmark/scripts/irmsd.py` unchanged. Retry only the affected
  identities with `benchmark/scripts/irmsd_grouped_safe.py`, then merge by
  `(dataset_row_id, source_model_sha256)` into a fresh output root.

Reason:
- The 96 retry failures were legacy symmetric-chain `KeyError`s, not DockQ
  score failures. Reusing the legacy script reproduces the defect and leaves
  valid DockQ rows without auxiliary iRMSD.

Consequences:
- The original dependency chain `1394212` -> `1394260` -> `1394261` is
  historical: `1394212` completed, `1394260` failed on 22 timed-out auxiliary
  iRMSD rows, and `1394261` never ran because its dependency became
  unsatisfiable; that obsolete pending job was canceled after final2 passed.
- The optimized safe wrapper and repository-local DockQ environment were used
  for the successful repair job `1404093`; its 22-row output was merged into
  `full-v2-safe-irmsd-final2/` and passed the final audit.
- Final acceptance requires and now has 6,539 model rows, 5,827 score-bearing
  rows, 712 `not_scoreable` rows, 192 cross-only scopes, zero failed auxiliary
  iRMSD rows, and zero hash or interface-contract failures.
- The optimized wrapper caches aligned chain-pair results and residue-index
  interface masks across symmetry orders. This is a performance repair around
  the safe path; `benchmark/scripts/irmsd.py` remains unchanged.

### Current canonical ranking result

Decision:
- Treat the final2 ranking audit as the current verified ranking result, while
  keeping ranking opt-in in the stable pipeline.

Evidence:
- `audit.json` passed with 5,828 candidates, 5,827 eligible/rankable labels,
  one retained explicit non-scoreable candidate, 155 groups, and no provenance
  or identity errors.
- `baseline-evaluation.json` reports 58/155 native-like top-1 groups versus
  68/155 oracle-positive groups; median top-1 DockQ is 0.02234852 and median
  best DockQ is 0.09339036.

Consequence:
- This validates the scoring/ranking contract and provides a baseline result;
  it does not establish a quality improvement or change the production default.

### Ranking eligibility contract

Decision:
- Retain explicit `not_scoreable` candidates in ranking tables for provenance,
  but exclude their identities from rankable candidate-vs-score equality. Fail
  on any other missing or unexpected label.

Reason:
- `rigid:000023` is a genuine non-bijective partner-cardinality case and has
  no canonical score by design; treating it as a missing score falsely failed
  the full ranking audit.

### Pipeline comparison evidence boundary

Decision:
- Use `tmp/agent/20260727-pipeline-comparison/pipeline_stage_ledger.csv` as
  an artifact-count and PDB-integrity ledger, not as a causal quality
  comparison.

Reason:
- Retained artifacts include current TMalign+external-Rosetta and
  TMalign+FiberDock arms, mixed current MultiProt experiments, and a legacy
  MultiProt+FiberDock workspace, but not a clean matched run across all
  inputs, templates, thresholds, source gates, and evaluators.

Consequences:
- Report observed stage attrition and execution limitations, but do not
  attribute quality differences to alignment, surface, or refinement until the
  paired contract is regenerated.

### Matched stage diagnosis boundary

Decision:
- Use `matched_candidate_ledger.json` and `multiprot_gate_diagnosis.json` as
  the detailed diagnostic evidence for the retained current runs. Treat the
  broad PDB counts as artifact counts only.

Evidence:
- The current TMalign external-Rosetta and FiberDock arms have byte-identical
  inputs, template manifests, alignment JSON, and transformed PDBs. Their ASA
  surface PDBs are identical; RSA differences are limited to NACCESS total-line
  formatting/rounding. Candidate-level divergence begins in refinement output:
  Rosetta has 5 canonical, 1 partial, and 1 missing candidate; FiberDock has 7
  energy PDBs, with one missing `zero-trial` error in its log.
- The retained MultiProt run has 39 successful individual alignments, only 2
  paired successful orientations, and 0 paired orientations passing the
  current transformation gate. Its Kabsch-RMSD-derived TM-score proxy is not
  calibrated to the TMalign threshold of 0.5.

Consequence:
- Do not attribute the current MultiProt loss to refinement or declare a
  refiner-quality difference. First calibrate/validate the MultiProt score
  contract, then replay matched transformed candidates through both refiners.
- Preserve the external-Rosetta partial/missing cases as an observability
  issue until per-candidate subprocess return codes and score-gate reasons are
  retained.

### Current-tree ranking evaluation boundary

Decision:
- Keep ranking opt-in. Accept the ranking integration as operational and
  report candidate/refinement-load reduction, but do not claim general speedup
  or production accuracy improvement from the current evidence.

Reason:
- Current-tree job 1400311 completed both baseline and ranked TMalign+
  external-Rosetta arms with identical source/input/audit manifests; ranking
  selected 1 of 2 candidates and added a completed ranking stage. However,
  ranked refinement took 62.7s versus 56.3s baseline in that one stochastic
  Rosetta run. The current-code reconstruction of the 10-group labeled pilot
  selected all 3 oracle-positive groups, but the table is GTalign-derived and
  lacks the current audit coverage fields.

Consequences:
- Treat `tmp/agent/20260728-ranking-pipeline-evaluation/` as current evidence
  and preserve the older ranking CSV as historical evidence only. Require
  repeated matched runs with frozen seeds and current-pipeline labels before
  changing defaults or making causal accuracy claims.

### MultiProt score-contract repair boundary

Decision:
- Dispatch alignment eligibility by aligner contract. MultiProt uses its
  native match-count/coverage gates; TMalign and GTalign continue to use the
  shared TM_SCORE_THRESHOLD gate. Keep MultiProt's legacy RMSD-derived
  tm_score field only as a compatibility/diagnostic field and label its
  contract explicitly.

Reason:
- MultiProt does not emit the TMalign TM-score quantity. Applying the shared
  0.5 threshold to max(0, 1 - Kabsch_RMSD/10) rejected successful MultiProt
  records before transformation and produced 0/200 paired eligible
  orientations in the retained panel.

Evidence:
- The fix is implemented in src/transformation.py and the MultiProt JSON
  writer labels tm_score_contract and score_gate_contract.
- The retained 100-template diagnostic under
  tmp/agent/20260728-multiprot-gate-fixed-test/ reports 2/200 paired
  orientations eligible and 37/39 successful sides passing the native gate.
- The focused transformation/alignment regression suite passes 21 tests;
  modified files also pass py_compile.

Consequence:
- Stable TMalign defaults and ranking behavior are unchanged. The fix removes
  only the incompatible cross-aligner score comparison. The two eligible
  MultiProt orientations still require transformation/clash and native DockQ
  evaluation before any quality claim.

### Exact matched replay and FiberDock output-contract boundary

Decision:
- Treat `tmp/agent/20260728-matched-align-refiner-replay/results-v2/` and job `1404934` as the current causal stage evidence. Do not call the replay a quality benchmark or change stable thresholds from it.

Evidence:
- The MultiProt and TMalign arms used the same `1fgnHL`/`1tfhA` inputs, `1ahwAF/o1` candidate, source PDBs, transform implementation, clash settings, and Rosetta/FiberDock entry points. TMalign generated the common pair; MultiProt returned `clash_rejected` after writing valid intermediates.
- MultiProt-vs-TMalign transformed CA direct RMSDs were `59.8668` and `45.8653` Å for 428 and 202 common CA atoms. MultiProt had 16 cross-partner CA clashes with minimum `1.478` Å; TMalign had zero with minimum `6.015` Å.
- The diagnostic MultiProt `1ahwBC/o1` replay produced valid FiberDock and raw Rosetta PDBs. Its FiberDock `.ref` file is `fiberdock_energies.ref` and contains `glob = 0.00`, but `src/fiberdock_refinement.py` parses `fd_params.ref`, so the returned energy is incorrectly `-`. Its Rosetta raw structure had interaction score `0.0`, so the existing `-5.0` gate did not retain a canonical flat output.

Consequence:
- The first matched divergence is the alignment-to-transform geometry/clash gate; refiner quality is not separable for the common candidate because MultiProt never reaches refinement.
- The FiberDock energy parser is a confirmed output-contract bug in the current implementation. Repair it only in an isolated branch with a focused test before changing the stable pipeline.
- Preserve the failed setup job `1404921`; use `1404934` for current evidence.

### FiberDock output-contract diagnostic boundary
- The current `src/fiberdock_refinement.py` parser derives `fd_params.ref` from the parameter-file stem, but FiberDock uses the declared `energiesOutFileName` stem and writes `fiberdock_energies.ref`. The isolated diagnostic found valid declared solution files and valid refined PDBs, so this is confirmed as a parser-path defect. Do not patch the stable source until a focused regression test and isolated corrected replay are available.

The corrected FiberDock replay completed in job 1405134 with both candidates returning code 0 and valid 6,765-atom outputs. The focused fixture regression `tests/test_fiberdock_output_contract.py` passes. This validates the isolated output-path correction as a reproducible contract fix, but does not promote it into stable source or establish a structural-quality result; source integration remains an explicit opt-in change.

### Broader matched transformed-candidate panel boundary

Decision:
- Use jobs 1405218 and 1405237 under
  `tmp/agent/20260728-multiprot-fiberdock-broader-replay/` as the current
  broader diagnostic panel, not as a method-quality benchmark.

Evidence:
- The seven selected candidates have matching input and source-PDB inventory
  hashes, but different MultiProt/TMalign alignment and surface inventories.
  Direct CA geometry differs for 3/7 left partners and 1/7 right partners;
  total cross-partner CA clashes are 454 versus 385 at 3 Angstrom.
- MultiProt FiberDock produced seven valid declared-energy PDBs. The current
  parser returned `-` for all seven because it still looks for
  `fd_params.ref`; the declared output is `fiberdock_energies.ref`.
  TMalign external Rosetta produced five canonical, one partial, and one
  missing candidate.

Consequence:
- The panel localizes observed differences to aligner-dependent transform
  geometry plus downstream output acceptance, but does not establish causal
  structural quality. A native DockQ comparison requires a validated
  MultiProt gate, per-candidate Rosetta return/score observability, matched
  refinement controls, and a frozen evaluator.

### Phase 1 run identity and manifest planning boundary

Decision:
- Give each operational attempt a readable `run_id` and each declared run
  contract a canonical manifest hash. Store the contract in canonical JSON
  with a TSV artifact ledger keyed by durable row identity, scientific role,
  and run-relative path.
- Permit dirty worktrees but record exact Git state fingerprints. Preserve raw
  and normalized selectors, hash all materialized scientific artifacts,
  preserve symlink metadata while hashing target bytes, append artifact records
  during stages, and finalize them before downstream consumption.
- Fail closed on hash or row-identity mismatch, preserve diagnostics, and
  represent retries as linked immutable attempts. Expose validation as both a
  reusable library and a CLI gate.
- Capture a secret-safe allowlist, package/tool identity, seeds, structured
  argv, launcher hash, and a unified local/Slurm execution context.

Reason:
- Phase 1 must make provenance and artifact identity trustworthy before later
  stage, candidate, evaluation, HPC, or evidence-bundle contracts consume it.

Canonical context:
- `.planning/phases/01-run-identity-and-manifest/01-CONTEXT.md`

### Phase 1 grilling decisions: contract, attempt, and closure boundary

Decision:
- Separate `contract_hash` from `run_id`: the contract hash covers declared
  source/tool/config/input identity, while execution facts and evolving
  artifacts belong to the attempt and closure records.
- Require explicit `dataset_row_id` for scored/benchmark runs; permit only
  labeled synthetic IDs for exploratory runs.
- Drive the artifact ledger from an expected inventory, record scientific
  extras as unclassified, use row/role/path identity, and enforce validation
  for new manifest-aware consumers with an explicit legacy-unverified bypass.
- Preserve TSV inspection with per-row and whole-ledger digests. Operational
  retries retain the contract hash; corrected declarations receive a new hash.

Canonical records:
- `docs/adr/0001-contract-and-attempt-identity.md`
- `docs/adr/0002-row-aware-artifact-ledger-and-consumer-gate.md`
- `docs/glossary.md`

### Phase 1 evidence-ledger deepening boundary

Decision:
- Make `src/provenance/run_evidence.py` the canonical standard-library module
  for contract/attempt records, row-aware artifact observations, closeout, and
  the fail-closed consumer gate.
- Keep `benchmark/scripts/investigation_*` as compatibility adapters during the
  migration; preserve direct path-based CLI execution and existing output
  formats.
- Prove the seam through an isolated CLI/temporary-run fixture before adding an
  opt-in `prism.py` hook. Do not merge candidate lineage, ranking, evaluation
  joins, or stage lifecycle into the Phase 1 module.

Evidence:
- Focused Phase 1 and compatibility tests pass: 42 tests.
- Direct `python benchmark/scripts/investigation_provenance.py --help` exits 0.
- `tests/test_run_evidence.py` covers contract/attempt identity, secret
  redaction, symlink target hashing, mutation rejection, missing records,
  duplicate keys, row mismatch, closeout, and the CLI wrapper.

Consequence:
- The next migration decision is whether template preflight and the remaining
  benchmark artifact writers should move into `src/provenance` or remain
  format-specific adapters. Do not expand the current Phase 1 core until that
  ownership is explicitly chosen.


### PRODIGY ranking boundary
**Date:** 2026-07-30
**Type:** integration
**Status:** accepted

Decision:
- PRODIGY is an explicit opt-in ranking method selected with
  `--rank --rank-method prodigy`; it does not change the deterministic baseline
  or stable unranked defaults.
- PRODIGY runs through an externally installed executable/environment and its
  combined inputs, stdout/stderr, command, return code, and score metadata are
  retained under `processed/ranking/prodigy/`.
- If any candidate in a receptor/ligand group cannot be scored, the complete
  group remains selected so partial external-tool output cannot become an
  implicit scientific decision.

Reason:
- PRODIGY 2.4.0 requires FreeSASA and NumPy >=2, while the verified DockQ
  environment has a separate pinned dependency contract.
- Its predicted affinity is exploratory ranking evidence, not a native DockQ
  label or a demonstrated quality improvement.

Consequences:
- Installation is documented but not performed by the pipeline implementation.
- Independent frozen native-complex evaluation is required before enabling this
  method for production claims.

### Learning workflow boundary
**Date:** 2026-07-30
**Type:** workflow
**Status:** accepted

Decision:
- Use `mentoring-juniors` for Socratic guidance and `teach` for persistent
  lessons, reference documents, and learning records when learning ongoing
  PRISM concepts.
- Do not use the deleted learnship learning skill for this teaching workflow.

Reason:
- The user explicitly selected these two teaching tools as the replacement
  learning workflow.

Consequences:
- Teaching sessions should preserve learning state through the `teach`
  workspace and use questions and progressive clues from `mentoring-juniors`.
- This decision concerns concept learning only; it does not supersede the
  project's task-planning or repository-routing instructions.

### MultiProt true TM-score contract and calibrated gates
**Date:** 2026-07-30
**Type:** architecture
**Status:** accepted

Decision:
- MultiProt alignment JSON now contains both `tm_score` (proxy: 1-RMSD/10) and `true_tm_score` (standard length-normalized TM from match_dict via Kabsch alignment of matched CA pairs).
- Transformation gate for MultiProt uses calibrated thresholds: `true_tm_score >= 0.3 AND match_count >= 10 AND match_pct >= 30%`.
- TMalign/GTalign retain the shared `TM_SCORE_THRESHOLD` (default 0.5) on their native TM-scores.
- Diagnostic thresholds do not become production defaults.

Reason:
- The proxy TM-score was incomparable with TMalign scores and caused 97% candidate loss at the transformation gate (0/200 pairs passed in 100-template test).
- Diagnostic SBATCH (job 1416942, 8 CPUs, 1 ai QoS slot) with seccomp bypass yielded 47/56 successful pairs and exposed true TM distribution (median 0.016, max 0.604).
- True TM-score formula: TM = (1/L) Σ 1/(1+(d_i/d0)²) where d0 = 1.24(L-15)^(1/3)-1.8, computed from Kabsch-aligned match_dict residue pairs.

Consequences:
- MultiProt candidates can now pass transformation and reach refinement for the first time.
- Any causal pipeline comparison MUST use matched inputs, assets, refiners, and native DockQ—alignment gate passage alone is not a quality signal.
- Stable TMalign thresholds remain unchanged.

### Graphify analysis findings

Decision:
- Graphify analysis (2026-07-28) identified `split_target_id()` in `pdb_download.py` as a god node (12 edges, cross-community hub connecting download, compare, transformation).
- Low-cohesion communities: TMalign.cpp (0.10), compare.py (0.07), hotspot.py (0.07) suggest refactoring opportunities.
- No import cycles detected in codebase.

Reason:
- Architectural bottlenecks identified for targeted refactoring.
- God node indicates tight coupling between download, compare, transformation.
- Low cohesion suggests modules doing too many unrelated things.

Consequences:
- Consider refactoring `split_target_id()` usage to reduce coupling.
- Investigate TMalign.cpp, compare.py, hotspot.py for cohesion improvements.
- Graphify output preserved under `graphify-out/` for future reference.

### Feature branch workflow

Decision:
- Major features (ranking, prescript integration) developed on `feature/ranking-progression-system` branch (PRISM-prescript) and `pre_scripts` branch (PRISM-main-archive).
- Main branches kept clean. Merge only after end-to-end validation.

Reason:
- Isolates experimental work from stable main branch.
- Enables parallel development and clean integration.

Consequences:
- PRISM-prescript `feature/ranking-progression-system` branch contains ranking implementation (commit 33b0551).
- PRISM-main-archive `pre_scripts` branch contains prescript integration (commit 2571180).
- Both main branches remain clean.

### PRISM-prescript to PRISM-main-archive integration

Decision:
- 9 prescript modules integrated into PRISM-main-archive `pre_scripts` branch: `candidate_audit.py`, `candidate_ranker.py`, `candidate_selector.py`, `ranking_data.py`, `ranking_metrics.py`, `alignment_multiprot.py`, `pyrosetta_refinement.py`, `fiberdock_refinement.py`, `template_filtering.py`.
- PRISM-main-archive `prism.py` updated with full CLI integration.
- Missing imports fixed in PRISM-main-archive: `normalize_target_id` in `pdb_download.py`, `extract_chain_and_res_ids`/`_has_ca_atoms`/`_write_empty_alignment` in `alignment.py`.

Reason:
- PRISM-main-archive needed prescript's advanced features for benchmark comparison work.
- Integration preserves PRISM-prescript as scoring/comparison layer.

Consequences:
- PRISM-main-archive `pre_scripts` branch now has full prescript feature set.
- PRISM-prescript `feature/ranking-progression-system` branch also updated.
- Both main branches remain clean.

### Cross-repository CLI compatibility policy
**Date:** 2026-07-31
**Type:** compatibility
**Status:** accepted

Decision:
- PRISM-prescript accepts the `/PRISM` option names while preserving existing
  underscore/hyphen spellings as argparse aliases.
- Prescript refinement remains enabled when neither `--refine` nor
  `--no-refine` is supplied; `--no-refine` is the explicit skip control.
- Backend paths and runtime settings are passed to the consuming functions,
  rather than treated only as import-time environment configuration.

Reason:
- Commands should transfer between the repositories without invalid-option
  failures, but changing prescript's established default refinement behavior
  would invalidate stable run expectations.

Consequences:
- Option spelling is compatible, while repository-specific defaults remain
  documented rather than silently unified.
- New tests must cover both legacy prescript forms and `/PRISM` forms.
