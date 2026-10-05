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

### Curated architecture graph scope - 2026-09-17

Decision:
- Use the maintained `src/` tree as the bounded architecture-navigation scope
  for the current Graphify artifact. Keep generated benchmark/history trees,
  copied environments, caches, and Graphify-generated corpus files outside the
  architecture graph.
- Treat the resulting `graphify-out/graph.json` as navigation evidence only;
  exact execution and scientific decisions remain governed by source,
  manifests, status records, and validated run artifacts.

Reason:
- The full checkout is dominated by historical/generated evidence and did not
  finish a useful corpus audit within the bounded inspection window. A broad
  graph would obscure current module relationships and repeat the stale-graph
  problem already documented in project memory.

Evidence:
- The 2026-09-17 source graph contains 584 nodes, 1,113 edges, and 26
  communities across 47 code files. Graphify diagnostics report zero missing,
  dangling, self-loop, or collapsed edges.

Consequence:
- Refresh the source graph incrementally when architecture changes. Build a
  separate deliberately filtered evidence/chronology graph when dated run
  provenance is needed; do not silently merge it into the source graph.

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
- Implemented on `feature/prism-cli-parity` with reversible checkpoint commits
  followed by `20ebf4c3cc3`; the focused compatibility/backend suite passed
  31 tests. Unrelated dirty worktree files were intentionally not staged.

### Benchmark readiness boundary
**Date:** 2026-08-13
**Type:** workflow
**Status:** accepted

Decision:
- Use current TMalign + NACCESS + external Rosetta as the stable baseline for
  the next benchmark plan.
- Freeze source state, template/interface assets, binaries, environments,
  configuration, input manifest, and evaluator before launching full runs.
- Treat GTalign CPU/GPU, current MultiProt end-to-end, and current-versus-
  legacy FiberDock comparisons as distinct validation arms rather than silently
  combining them into one benchmark result.

Reason:
- The repository currently contains extensive generated/untracked benchmark
  state, and several tools are available without having a fully causal,
  end-to-end comparison contract.

Consequences:
- Do not claim all pipeline variants are benchmark-validated from binary
  presence or partial retained outputs.
- Benchmark manifests must record source/tree identity, assets, executable
  hashes, environment paths, thresholds, Slurm job IDs, and evaluator version.

### Graphify main-files corpus
**Date:** 2026-08-21
**Type:** tooling/provenance
**Status:** accepted

Decision:
- Keep the navigation graph in `graphify-main-files/` as a curated, code-only
  merge of the live current TMalign/Rosetta pipeline and the retained legacy
  MultiProt/FiberDock pipeline.
- Include current `prism.py`, all `src/*.py`, stable-pipeline documentation,
  and the legacy CLI plus `run_files/*.py` and launch/configuration files.
- Exclude generated results, copied environments, caches, benchmark artifacts,
  and external binary payloads.

Reason:
- A full repository scan is contaminated by nested historical runs, vendored
  environments, and generated outputs; it produced an unusable oversized graph.
- The graph is for architecture navigation, not source-of-truth provenance or
  scientific validation.

Consequences:
- The merged graph currently contains 556 nodes and 1,100 edges.
- Refresh it by rebuilding the two code-only arm graphs and merging them; do not
  run incremental update against the older pre-#1504 archived graph.

### Feature-bundle demo notebook scope
**Date:** 2026-08-21
**Type:** workflow
**Status:** accepted

Decision:
- Treat `tmp/prism-prescript-pipeline-extension-clean/notebooks/prism_pipeline_demo.ipynb`
  as a read-only diagnostic/demo artifact for the isolated feature bundle.
- Use current source files, stable-pipeline documentation, and terminal stage
  evidence as the authority for implementation or validation claims.

Reason:
- The notebook inspects retained runs, backend outputs, and source differences;
  its expensive execution cells are explicitly opt-in and the bundle contains
  intentionally Git-ignored external payloads.

Consequence:
- The notebook is useful for demonstration and orientation, but successful
  notebook inspection alone is not pipeline acceptance evidence.

### Transformation gate diagnostics preserve independent outcomes
**Date:** 2026-08-23
**Type:** scientific workflow
**Status:** accepted

Decision:
- Diagnose transformation attrition with an independent gate ledger before
  interpreting a cumulative filter order.
- Use the cumulative diagnostic order: record availability, transform fields,
  minimum matches, interface coverage, aligner-specific score contract,
  hotspot mapping, complementary contacts, transformation materialization,
  then C-alpha clashes.
- Keep a separate replay of the actual production order and leave all stable
  thresholds and filter-mode defaults unchanged.

Reason:
- Cumulative attrition is order-dependent. If the TM-score gate removes the
  cohort early, later protocol and geometry conditions cannot be observed even
  when they are independently measurable.

Consequences:
- Use `notebooks/transformation_gate_ablation.ipynb` for interactive,
  one-variable-at-a-time diagnosis.
- Treat missing assets as `not_evaluable`, and protocol gates as not applied in
  geometry-only mode; neither state is biological success.
- Do not infer threshold calibration from the single retained demonstration
  case, and do not compare TMalign scores directly with MultiProt's native
  match/coverage contract.

### ProInterVal integration boundary - 2026-08-25

Decision:
- Treat ProInterVal as an analysis-only, opt-in auxiliary interface score.
  Do not enable it as a default template curation filter, pre-refinement gate,
  or replacement for physics-based ranking without a matched validation.

Reason:
- ProInterVal validates learned interface plausibility, whereas PRISM generates
  and refines template-derived complexes. Its reported training labels and
  decoy distributions do not establish performance on PRISM-generated rigid or
  flexibly refined candidates. Recent scoring literature also shows that hit
  rates can be unstable under model sampling and that score-to-DockQ
  correlations degrade as conformational change increases.

Consequences:
- Any future adapter must preserve candidate identity, model/native hashes,
  chain roles, preprocessing parameters, model version, score semantics, and
  explicit unavailable/failed states.
- Validation must use frozen PRISM inputs and candidates with native DockQ/iRMSD
  labels, family/structure-disjoint splits, calibration, per-target metrics,
  and false-negative auditing. No default threshold is authorized by this
  literature review.

### VALAR bounded worker run blocked before project work - 2026-09-07

Decision:

- Preserve the PRISM-prescript working tree and canonical outputs unchanged.
- Do not claim USalign integration or transformation attrition findings from
  run `20260907-prism-prescript-usalign-transform`.
- Treat the local Copilot adapter issue as an execution blocker: Copilot 1.0.80
  registered the run-scoped offline session but did not dispatch the required
  `-i` prompt while its auth state was `Logged out`. Switching to `-p`, adding
  credentials, or enabling remote authentication was not authorized.

Evidence:

- Durable blocker diagnosis:
  `tmp/agent/20260907-prism-prescript-usalign-transform/orchestrator/run-050ecc271b334e3eb935363a0f45d7c8/evidence/copilot-interactive-blocker.md`
- Two same-session retry jobs loaded the validated A40 Q6 profile and passed
  local Qwen health/model checks, but recorded zero Copilot turns and zero
  project edits before cancellation.

Consequence:

- The maintained NACCESS/TMalign/default path and all pre-existing dirty state
  remain authoritative. A future bounded run requires an explicitly approved
  offline-compatible Copilot invocation or supported local authentication
  mechanism before scientific work is scheduled.

### VALAR prompt-mode smoke blocked before project work - 2026-09-07

Decision:

- Preserve the fresh bounded run
  `tmp/agent/20260907-prism-prescript-usalign-attrition-prompt/` and do not
  launch USalign or transformation workers from it.
- Treat Qwen health/model discovery as preflight evidence only. Do not claim
  local-provider usage or `/goal` activation when Copilot has produced no
  model turn.
- Do not retry automatically, resume the failed interactive `-i` launch, or
  fall back to remote/authenticated/alternate models.

Evidence:

- Slurm job `1653435` on `ai26` loaded validated A40 Q6, passed health, and
  returned the actual `/v1/models` ID on port `18192`.
- Copilot runtime `1.0.80` exited status 1 after only
  `session.mcp_servers_loaded` and `session.skills_loaded`; its session
  database contains zero `turns` rows and the Qwen log contains no inference
  request.
- Run-scoped provider binding and exact hook-disabled settings are recorded,
  but they are configuration evidence rather than successful inference.

Consequence:

- USalign availability, integration, and the retained 1gte transformation
  attrition figures remain unrun and unverified against current files.
- The next action is a separately authorized smoke-only investigation of the
  Copilot 1.0.80 batch/provider dispatch, with this run retained as evidence.

### VALAR repaired-adapter provider smoke blocked before project work - 2026-09-08

Decision:

- Preserve the fresh run `tmp/agent/20260908-prism-prescript-usalign-attrition-prompt/` and do not launch USalign or transformation workers from it.
- Treat Qwen health, `/v1/models`, and a server-side inference request as preflight evidence only. The repaired batch-safe Copilot launch is not validated because the exact direct-provider sentinel was not proven.
- Do not retry automatically, alter provider-smoke validation, or fall back to remote/authenticated/alternate models in this run.

Evidence:

- Slurm job `1653627` on node `ai25` loaded validated `a40-q6-262k-1gpu` Q6 on port `18682`.
- `logs/qwen-server.log` records one direct local inference request after model load.
- `logs/provider-smoke.json` records `FAILED` with `error_type: ValueError`; response content was intentionally not persisted. Copilot was not started, so no compute-node Copilot runtime or model-turn evidence exists.
- Parent and worker states are `BLOCKED`/`FAILED`; no project worker was launched and prior blocked runs remain preserved.

Consequence:

- USalign availability/integration and current transformation-file attrition remain unrun and unverified in this cycle.
- The next action is a separately authorized minimal provider-smoke diagnosis that preserves response-body non-persistence, followed by a fresh single-worker smoke gate.

### Explicit orientation-comparison run controls - 2026-09-09

Decision:

- Keep `native` as the default/no-option arm and define it as evaluation of
  both implicit template-chain assignments. Retain `o1` and `o2` as explicit
  fixed-orientation comparison arms.
- Add explicit per-run input/template selection and typed threshold overrides
  without changing the established defaults or the separate MultiProt native
  match/coverage score contract.
- Keep notebook execution opt-in and isolate each arm in a fresh run root.

Rationale:

- Comparing orientation requires the same input CSV, template panel, aligner,
  and thresholds while changing only the branch selection. The no-option arm
  must therefore preserve both branches rather than silently select one.
- Environment variables remain supported for compatibility, but CLI values
  must be observable in the run configuration and reach the actual gates.

Consequence:

- Candidate-yield changes from relaxed thresholds are diagnostic observations,
  not validated scientific improvements. A benchmark ledger must still retain
  raw alignment, qualified, contact, transformation, clash, and downstream
  outcomes by orientation.
- MultiProt coverage is measured against the corresponding template-chain
  interface size when that asset is loaded; falling back to mapped-residue
  count is reserved for direct records without a registered template size.

### Stepwise reliability evidence boundary - 2026-09-10

Decision:

- Keep the orientation comparison notebook as the controlled entry point for
  stepwise diagnosis. It must freeze the stable thresholds (including the
  5.0 Å surface scaffold), write per-arm provenance manifests, and retain
  explicit missing/unknown statuses.
- Treat gate ablation, clash grids, refinement checks, ranking comparisons,
  and US-align checks as diagnostic or separate validation arms. None changes
  production defaults or converts survivor counts into quality claims.

Rationale:

- Candidate-yield differences cannot be interpreted while input/interface
  assets, alignment score contracts, or per-candidate downstream evidence are
  unresolved.
- TMalign and MultiProt now retain raw-output hashes and return codes in their
  alignment records, allowing the notebook to distinguish complete provenance
  from a merely parseable JSON file.

Consequence:

- A completed notebook run is still required before making claims about
  orientation-induced candidate reduction, clash-threshold validity, or
  native-like recovery. US-align remains unavailable as a current `prism.py`
  aligner until its executable, output contract, and matched downstream arm
  are validated.

### Audit-contract hardening after independent review - 2026-09-10

Decision:

- Include interface-list JSONs in consumed-asset provenance and classify asset
  origin by both resolved path and content hash. If any required consumed asset
  is missing or cannot be matched to current/legacy reference bytes, preserve
  an unresolved status.
- Treat alignment return/status/hash contract failures as unknown evidence for
  every alignment-dependent gate. Numerical fields in a parseable but failed
  alignment record are not sufficient for a cumulative pass.
- Suppress native-like ranking recovery unless candidate panels are exactly
  matched and the label mapping has independent source, source hash, and a
  deterministic mapping hash.

Rationale:

- Transformation uses interface-list JSONs to define the coverage denominator;
  omitting them could hide a mixed or stale template panel. Alignment output
  can contain plausible numbers after a failed subprocess or truncated raw
  output. Ranking recovery without exact panels or linked label provenance
  confounds selection with quality.

Consequence:

- Existing/legacy runs lacking the new raw hashes or return codes remain useful
  for inventory, but their alignment-dependent gate outcomes are explicitly
  unknown. A biological conclusion requires rerunning or repairing the
  provenance contract before threshold decisions.

### Cross-repository boundary and validation hardening - 2026-09-12

Decision:

- Keep `PRISM-prescript` as the maintained pipeline and scientific-evidence
  authority; accept features from `PRISM` only through reviewed adapters and
  prescript contracts.
- Preserve explicit `warn` results for artifacts recorded as missing or
  unavailable, but require validation CLI consumers to return nonzero for
  `warn` and `fail` so incomplete evidence cannot proceed silently.
- Treat duplicate-ledger structure errors as controlled `fail` validation
  results with machine-readable output rather than an unstructured traceback.

Reason:

- The two repositories have different runtime contracts and the experimental
  branch lacks equivalent durable benchmark/provenance validation.
- Existing Phase 1 schemas distinguish `pass`, `warn`, and `fail`; changing
  the library status semantics would conflict with the accepted contract, while
  nonzero CLI behavior preserves fail-closed consumption.

Consequences:

- `docs/adr/0003-prism-repository-boundary.md` is the durable cross-repository
  ownership record.
- `src/validation_gate.py` now catches duplicate-ledger errors and supports
  `--output` as an alias for `--out`.
- The validation CLI now accepts a declared expected inventory for detecting
  rows absent from the ledger; pipeline producers still need to emit and pass
  that inventory consistently on every benchmark run.

### Smoke lifecycle and provenance boundary - 2026-09-12

Decision:

- Treat skipped as a terminal stage state when refinement is disabled or no
  candidates pass transformation; completion classification must require a
  terminal event, not merely a directory.
- Keep declared inventories separate from observed ledgers. The validation
  CLI supports TSV/JSON expected inventories and returns nonzero for warning
  or failure so missing rows cannot reach scoring silently.
- Represent dirty source state with bounded status entries and hashes rather
  than expanding every file in copied environments. Runtime observations
  remain outside the immutable contract hash.
- Keep CPU and GPU GTalign runs as separate arms until their raw-output,
  parameter, and parser differences are explained.

Evidence:

- Jobs 1658925, 1658928, and 1658929 have isolated stage logs and completion
  JSON under tmp/agent/.
- Job 1658926 has PyRosetta/FiberDock status records and an explicit DockQ
  ABI failure under tmp/agent/20260912-optional-backend-smoke/.

Consequence:

- Zero-pair smokes establish plumbing and backend observability only. They do
  not resolve orientation, threshold, ranking, or biological-quality claims.

### Evaluator output-contract hardening - 2026-09-12

Decision:

- Make explicit model/native chain selectors authoritative in the benchmark
  scorer; infer chain order only when selectors are absent.
- Validate the raw model chain contract before launching iRMSD or DockQ.
- Anchor helper-script resolution to the repository root rather than the
  caller's current directory.
- Give every DockQ invocation an isolated JSON output path and retain the
  interface count; when more than one interface is present, keep detailed
  interface metrics null at the model level.

Reason:

- A malformed legacy PDB could previously raise during chain inference or
  reach an external scorer, while a caller-provided mapping was ignored. A
  fixed/stale JSON path could also cause output reuse or overwrite.

Evidence:

- `tests/test_model_output_integrity.py`: 6 passed after the repair.
- Full prescript suite: 344 passed, 6 skipped.

Consequence:

- This is evaluator integrity hardening only. It does not make any current
  benchmark row scoreable and does not resolve the DockQ environment ABI
  blocker.

### DockQ runtime selection and replay boundary - 2026-09-12

Decision:

- Use `benchmark/prism_processed/env/prism_score_env/bin/python` as the
  repository-local DockQ scoring interpreter for isolated replays until a
  fresh environment matching `environment.yaml` is built and independently
  verified.
- Invoke DockQ as `bin/python -m DockQ` or through
  `benchmark/scripts/score_single_prism_pair.py`; do not call the stale
  repository-local `bin/DockQ` shebang directly.
- Keep `gtalign_env` unchanged. Its failed DockQ ABI attempt remains preserved
  as a blocked run, while the compatible repository-local replay is a separate
  evaluator arm.

Reason:

- The repository-local interpreter successfully imports DockQ `2.1.3` with
  NumPy `1.26.4`, whereas `gtalign_env` previously failed through a compiled
  extension against NumPy `2.4.6`. The two environments must not be conflated.

Evidence:

- Slurm job `1658942`, run root `tmp/agent/20260912-dockq-repo-env/`.
- Slurm job `1658943` passed through `src.eval.dockq` using the real
  `gtalign_env` pipeline interpreter plus the repository-local override.
- Slurm job `1658944` passed through the CLI adapter with
  `--dockq-json-dir`, retaining exactly one raw DockQ JSON artifact in the
  isolated run root.
- `tests/test_dockq_runtime.py`, `tests/test_compare.py`, and
  `tests/test_model_output_integrity.py`: focused runtime/compare/evaluator
  checks passed after the override and CLI addition.
- Raw and adapter return codes were both zero and both reported
  `GlobalDockQ/DockQ=0.2116967149685021` for mapping `OA:GF`.

Consequence:

- DockQ evaluator wiring is available for isolated scoring, but no cohort
  result, ranking conclusion, or biological quality claim is authorized until
  the benchmark denominator, native mappings, and per-row provenance gates are
  frozen.

### USalign parser and manual pilot - 2026-09-21

Decision:

- Recognize both TMalign `Chain_1/Chain_2` and USalign
  `Structure_1/Structure_2` score labels in the shared parser.
- Carry explicit `aligner_name` provenance into parsed records and set it to
  `USalign` only for the USalign runner branch; retain `TMalign` as the
  compatibility default.
- Keep full-panel USalign execution stopped until a paired pilot resolves
  speed flags, input shape, and score normalization.

Evidence:

- The old production records had nonzero mappings but zero TM-scores because
  the parser ignored Structure labels.
- Three single-chain manual pairs ran successfully through TMalign, USalign
  default/fast, and MultiProt, with no large USalign speed advantage.

Consequence:

- Parser validity is improved and tested, but the `tm_score=max(score_1,
  score_2)` compatibility contract must be explicitly accepted or replaced
  before scientific USalign comparison.

### PRISM matched-aligner recovery - 2026-09-28

Decision:

- Keep the historical 19,855-template TMalign/MultiProt ledgers as a labelled
  reusable lane and retain exact 19,948/19,062/19,058 panel evidence as a
  separate alignment-only lane. Never pool the panels silently.
- Reuse completed GTalign/TMalign/USalign exact alignment artifacts and do not
  repeat GTalign GPU search. Continue only from the first missing common
  downstream stage after panel and contract checks.
- Require a fresh corrected-parser USalign pilot for worker/configuration
  selection because the retained compact USalign JSON lacks explicit dual
  normalized scores and the old 946-labelled command provenance is ambiguous.
- Treat the one-query exact 946 pilot as a performance/alignment execution
  gate only; it cannot establish candidate quality or full-pipeline speed.

Evidence:

- `benchmark/prism_processed_results/prism_aligner_comparison_20260928/stage_ledger.md`
- `benchmark/prism_processed_results/prism_aligner_comparison_20260928/exact_panel_evidence/package_manifest.json`
- `benchmark/scripts/run_usalign_worker_sweep.py` and
  `benchmark/jobs/usalign_worker_sweep.sbatch`
- Slurm dry-run job 1708927 and bounded job 1708928; wrapper-only failed
  attempts 1708912 and 1708915 are preserved with explicit errors.

Consequence:

- USalign production may use default USalign with 16 measured workers, but
  `-fast` is not scientifically interchangeable on this pilot because its
  mappings, transforms, and scores differed materially. Production remains
  gated on resumable compact output and common downstream validation.

### USalign compact transformation/DockQ continuation - 2026-09-28

Decision:

- Run common transformation replay and compact candidate ledgers before
  removing raw USalign alignment JSON. Do not create an aligner-specific
  transformation pipeline.
- Score generated transformed pairs with the existing bijective evaluator,
  keeping complete GlobalDockQ distinct from requested cross-interface
  scores and preserving explicit score/failure statuses.
- Submit the corrected scoring array only as a dependency of compaction; each
  worker checkpoints rows and only then removes combined/scoring scratch and
  transformed halves.

Evidence:

- Focused tests: `33 passed` after adding the transformed scorer contract.
- `sbatch --test-only` accepted the corrected scoring array as job 1709044;
  production array 1709046 is queued after compaction job 1709007.
- Live production provider is job 1708992 with default USalign/16 workers;
  completion and scientific validation remain pending artifact evidence.

Consequence:

- No USalign quality or speed conclusion is promoted until provider
  completion markers, compact counts, transformed scoring outputs, and cleanup
  manifests validate. The exact-panel lane remains separate from the
  historical 19,855-template lane.

### Corrected refinement aggregation and retention gate - 2026-09-28

### USalign batch aggregation and refinement handoff - 2026-09-28

Decision:

- Add a read-only, dependency-gated batch aggregator for the 26 USalign
  production outputs. It must validate every per-batch compact status and
  TSV before emitting an aggregate; scheduler completion alone is not enough.
- Prepare the common-refinement manifest only from generated candidates with
  explicit refinable score states and existing transformed/native inputs.
  Keep all rejected rows and reasons in the compact ledger.
- Preserve USalign transformed halves through common refinement and paired
  DockQ validation; no cleanup is performed by the aggregation step.

Evidence:

- `benchmark/scripts/aggregate_usalign_batches.py`
- `benchmark/jobs/aggregate_usalign_batches.sbatch`
- `benchmark/scripts/prepare_usalign_refinement_manifest.py`
- Focused tests pass for both additions; selected validation suite is
  `46 passed`.
- Slurm test-only job `1709191` was accepted and production job `1709192` was
  submitted `afterany:1709046`.

Consequence:

- USalign remains `SUBMITTED` until provider, compaction, transformed DockQ,
  batch aggregation, and downstream refinement artifacts are independently
  validated. No run-scoped USalign deletion is currently permitted by the
  stage ledger.

### Final matched reducer contract - 2026-09-28

Decision:

- Join transformed and refined tables by the frozen candidate identity
  (case, template, orientation, query pair, and chain pair), with a short
  identity fallback only when chain labels are absent in both records.
- Compute `Delta_DockQ` only when both sides contain numeric GlobalDockQ for
  the same candidate. Use external Rosetta as the primary common-refinement
  arm and retain FiberDock as an explicitly labelled fallback/experimental
  arm.
- Report a normal-approximation 95% interval for paired deltas, with the
  method and small-sample rule recorded in the manifest; never use unpaired
  marginal means as a refinement effect.

Evidence:

- `benchmark/scripts/aggregate_final_matched_comparison.py`
- `tests/test_aggregate_final_matched_comparison.py`
- Focused test passes; source hash is retained in the staged provenance
  snapshot.

Decision:

- Treat checkpoint `dockq`/`dockq_sum` fields in the active common-refinement
  run as diagnostics only when raw DockQ JSON is available; use raw
  `GlobalDockQ` for complete-complex quality and explicitly aggregate only
  requested receptor-ligand interfaces.
- Preserve unresolved transformed structures after scoring failures or
  non-scoreable states. Delete transformed halves only after a compact score
  record reaches a validated score state, and record the retained unresolved
  count in the cleanup manifest.
- Join USalign compact candidates to the batch input manifest so case-level
  comparison uses durable BM5.5 pair identities rather than filename parsing.

Evidence:

- `benchmark/scripts/aggregate_corrected_refinement.py` and
  `tests/test_aggregate_corrected_refinement.py`.
- `benchmark/scripts/score_transformed_usalign_batch.py` and
  `tests/test_score_transformed_usalign_batch.py`.
- `benchmark/scripts/replay_compact_usalign_batch.py` and
  `tests/test_replay_compact_usalign_batch.py`.
- Existing GTalign corrected transformed score summary: 15,440 rows, 14,300
  scored, 1,140 score_failed across 216 cases; active refinement source
  example demonstrates `GlobalDockQ` can differ materially from `best_dockq`.

Consequence:

- Historical GTalign transformed scoring is reusable but panel-labelled and
  unrefined. Active KUACC refinement must finish before corrected T/M
  aggregation or any deletion. USalign production remains execution-pending
  validation.

### USalign refinement handoff gate - 2026-09-28

Decision:

- Submit a single dependency-gated preparation job after compact USalign
  aggregation. The job must verify `validated_compacted`, write selected and
  rejected manifests, and stop before array submission so the selected count
  and live resources can be reviewed.

Evidence:

- `benchmark/jobs/prepare_usalign_refinement.sbatch`
- `tests/test_prepare_usalign_refinement_job.py`
- Slurm test-only job `1709232` accepted; production job `1709233` submitted
  after `1709192`.

Consequence:

- The common USalign refinement stage is resumable and dependency-gated, but
  remains `NOT_SUBMITTED` until its validated manifest exists.

### USalign cleanup and timing provenance gate - 2026-09-28

Decision:

- Keep cleanup as a separate fail-closed operation after all downstream
  consumers validate. Require aggregate path/hash checks, refinement-handoff
  status, common-refinement cleanup eligibility, and final-package validation.
- Persist the dry-run cleanup manifest before applying deletion, reject
  symlinked or boundary-escaping targets, and retain compact aggregates,
  checkpoints, manifests, and provenance.
- Aggregate timing/resource ledgers separately from scientific quality tables;
  leave missing CPU/GPU fields empty rather than inferring resource use.

Evidence:

- `benchmark/scripts/cleanup_usalign_run.py`
- `tests/test_cleanup_usalign_run.py`
- `benchmark/scripts/aggregate_timing_resources.py`
- Comparison-focused regression set: `73 passed`.

Consequence:

- No USalign cleanup has been applied because provider, transformed-DockQ,
  common-refinement, and final matched consumers remain incomplete.

### Preserve nested common-refinement failures - 2026-09-28

Decision:

- Treat a common-refinement candidate with a terminal worker failure as an
  auditable result, not as a missing row or DockQ zero. Compact aggregation
  must carry failed stage names and their nested error/reason text.
- Do not synthesize a second predicted chain when the native mapping requests
  multiple ligand chains but the retained assembled model contains one.

Evidence:

- GTalign checkpoint
  `v2-chain-normalized-gtalign-medium_1wq1_045-3c6e1989873c8003.json`:
  input normalization failed because native ligand `G*` requires two chain
  segments while predicted model chain `D` provides one.
- `benchmark/scripts/aggregate_corrected_refinement.py`
- `tests/test_aggregate_corrected_refinement.py` (`5 passed` with cleanup and
  adapter regression tests).

Consequence:

- The failure will remain in the final failure/rejection summary and will be
  excluded from score denominators unless a scientifically justified mapping
  is independently established.

### Preserve active refinement ownership — 2026-09-28

Decision:

- Keep the existing VALAR GTalign array and KUACC TMalign/MultiProt refinement
  waves as the sole owners of their selected candidates while they are active.
- Do not submit overlapping replacement arrays or delete their inputs; advance
  only after checkpoint/exit markers and compact downstream records validate.

Evidence:

- GTalign job `1709167` remains active with fresh completed checkpoints.
- KUACC manifests and submission events identify the TMalign/MultiProt lane,
  while shard tasks continue progressing under the controller dependency.

### Use an opt-in common match/coverage contract for aligner comparison — 2026-09-28

Decision:

- Preserve provider-native acceptance as the default, but add
  `alignment_gate_mode=common_match_coverage` for controlled TMalign/USalign/
  GTalign/MultiProt comparisons.
- In common mode, use the shared minimum matched-residue and interface-coverage
  thresholds with an inclusive boundary for every aligner. Do not apply the
  TM-score threshold; retain scores only as post-hoc diagnostics.

Rationale:

- MultiProt's `tm_score` is an RMSD-derived proxy and cannot be made equivalent
  to a TMalign TM-score by assigning the same numeric cutoff. A common
  count/coverage contract is explicit, reproducible, and avoids a hidden
  algorithm-specific advantage.

Evidence:

- `src/transformation_config.py`
- `src/transformation.py`
- `src/stepwise_analysis.py`
- `prism.py`
- `docs/exec-plans/20260928-comparable-alignment-contract.md`
- 50 focused regression tests passed.
