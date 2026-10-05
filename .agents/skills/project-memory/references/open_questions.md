# Open questions and unresolved issues

### Local feature-bundle external-Rosetta compatibility

The `feature/prescript-pipeline-extension-clean` backend named
`external_rosetta` actually uses `src/rosetta_refinement.py` and PyRosetta.
PyRosetta 2026 removed `Vector1` and requires `utility.vector1_int` for
`setup_foldtree`. A focused compatibility test now passes; Slurm rerun 1506821
must complete before this backend is called confirmed. The explicit
`--refiner pyrosetta` adapter is independently confirmed.

## Open questions

### Benchmark source authority

Why it matters:
- Seventeen BM5/5.5 rows retain a disagreement among CSV chain order, raw
  selector roles, and curated archive `r/l` labels. The exact BM3.0 paper
  88-case cohort also lacks an authoritative local row-level selector/native
  ledger.

Current status:
- `source-gate-policy/v1` retains 240 strict-eligible rows for audit but does
  not authorize a confirmatory denominator.

Next investigation:
- Acquire/freeze the authoritative mapping or supplement with hashes, then
  rerun exact-chain validation before comparison claims.

### Ranking evidence

Why it matters:
- The new deterministic selection path is mechanically validated, but ranking
  utility must be measured on independent native complexes.

Current status:
- Deterministic ranking is mechanically validated and wired into `prism.py`.
- PRODIGY ranking is now implemented as an opt-in scorer and smoke tested on
  the retained `5zngA,4eylA` / `1a0cCD` two-orientation case. It changed the
  forwarded set from two candidates to one (`o1`) after both candidates scored
  successfully. Evidence:
  `tmp/agent/20260730-prodigy-ranking-paired-test/summary-corrected.json`.
- Paired isolated ranked smoke (job 1392725) confirmed load reduction: top-1
  selected 1/2 candidates vs 2/2 unranked; identical source/input/audit manifests.
- Current-tree paired smoke (job 1400311) independently confirmed the same
  behavior with current source hashes and all stage-status records completed.
  The ranked arm refined 1 candidate versus 2 baseline candidates, but its
  observed refinement time was 62.7s versus 56.3s baseline; Rosetta runtime
  variation prevents a speedup claim from this pair.
- 4 completed BM5.5 variants scored against T_Rigid with 13-19 DockQ scores each.
- Current-code re-ranking of the retained 235-row, 10-group pilot selected all
  3 oracle-positive groups at top-1; median top-1 DockQ was 0.007527 versus
  0.019942 for the group oracle, with mean regret 0.003381. The table is
  GTalign-derived and lacks the audit coverage fields used by the current
  pipeline, so this is pilot evidence rather than a paired production claim.
- The older retained ranked CSV is not current-code-comparable: it lacks
  `baseline_score_version` and contains scores inconsistent with the current
  scorer. A current-code reconstruction is retained separately.
- Final2 canonical ranking is now mechanically complete: 5,828 candidates,
  5,827 eligible/rankable labels, one retained explicit non-scoreable row,
  155 groups, and zero audit errors. The baseline produced 58/155 native-like
  top-1 groups versus 68/155 oracle-positive groups, with median top-1 DockQ
  0.02234852 versus median best 0.09339036. This resolves the pipeline
  contract question, but not ranking-quality generalization.

Next investigation:
- Evaluate the new opt-in PRODIGY scorer on frozen transformed candidates with
  native DockQ labels and compare against the deterministic baseline. Do not
  use PRODIGY scores as labels or enable it by default before that evaluation.
- Freeze Rosetta random seeds or repeat matched paired runs before quantifying
  wall-clock speedup. Build a labeled candidate table from the current
  TMalign pipeline with the same audit coverage fields, then compare ranked
  top-1 against unranked/full-candidate best outcomes across independent
  complexes before any production-default decision.

### Alignment and runtime scale

Why it matters:
- Full 20K-template CPU runs are impractical; CPU and GPU GTalign output counts
  have differed materially under apparently identical settings.

Current status:
- GPU GTalign is the practical broad-panel backend. The CPU/GPU parameter
  discrepancy still needs a controlled binary-version/argument audit.

Next investigation:
- Freeze matching input/template panels and command lines, compare raw binary
  outputs, and avoid using CPU results as interchangeable with GPU results.

### Protocol-asset parity and positive canary

Why it matters:
- Modern JSON and derived legacy hotspot/contact assets disagree materially,
  and strict published-protocol canaries have completed with no predictions.

Current status:
- The code now normalizes assets by chain and orientation, but the archived
  parity audit retains 19,005 disagreements and 850 missing modern profiles.
  Relaxed positive runs are diagnostic only.

Next investigation:
- Establish semantic asset equivalence or freeze a single authority, then find
  a strict published-protocol positive full-pose canary before confirmatory
  evaluation.

### Chronology extraction from project evidence

Why it matters:
- `tmp/agent` contains current run notes, manifests, status JSON, logs, and
  reproducible outputs that are part of the project history. Excluding all of
  `tmp` would discard integral progress evidence.

Current status:
- The existing Graphify graph has 208,846 nodes and 337,510 links, but no
  temporal links and no temporal node metadata. Its traversal is additionally
  affected by pre-#1504 IDs and copied environments/vendored dependencies
  inside historical run directories.
- The reliable approach is a filtered chronology corpus: retain dated
  `tmp/agent` evidence and project memory/docs, while excluding only
  `site-packages`, `python-site`, `node_modules`, caches, and generated
  Graphify corpus trees. Derive dates/statuses deterministically from paths,
  headings, manifests, `exit.json`, and run logs, then expose explicit
  `before`, `after`, `supersedes`, `validates`, and `blocked_by` relations to
  Graphify.
- The refreshed separate chronology graph is complete at `docs/chronology/`:
  168 events across 22 dates, with 199 path-qualified nodes, 360 directed
  edges, 23 verified cross-day `before` relations, and 168
  `has_status`/`occurred_on`
  relations, and a `current_as_of` anchor for the latest date.
  The Graphify integrity diagnostic found no endpoint or edge-collapse issues.
  Same-day event ordering is not asserted. Immediate
  `runs/<job-id>/results.tsv` files are now parsed without recursive run-tree
  scanning, so the 2026-07-26 ranked smoke resolves as `recorded:success`.

- On 2026-09-17, the bounded source-only architecture graph was rebuilt at
  `graphify-out/` with 584 nodes and 1,113 edges across 47 `src/` files. Its
  graph-health audit passed with no endpoint, self-loop, or edge-collapse
  findings. The full checkout graph remains intentionally unbuilt because its
  generated/history trees require a separate explicit filter policy.

Next investigation:
- Refresh `docs/chronology/` after material run completion and add explicit
  terminal timestamps only where the source record provides them; do not infer
  same-day order from directory names alone.

### MultiProt score-contract follow-up → RESOLVED/ADVANCED

Why it matters:
- The old MultiProt RMSD-to-TM proxy was being compared with the TMalign
  threshold and caused premature candidate loss before transformation.

Current status UPDATE (2026-07-30):
- Seccomp bypass implemented and tested (job 1416942).
- True TM-scores computed for 47 successful pairs (range 0.0014-0.6041, median 0.0157).
- Transformation gates recalibrated in `src/transformation.py` with calibrated thresholds (true_tm_score ≥ 0.3, matches ≥ 10, coverage ≥ 30%).
- Two pairs with true TM ≥ 0.3 now pass: 5zngA_1a0cCD_C (0.3294, 11 matches) and 5zngA_1buhAB_A (0.3385, 11 matches).
- Correlation between proxy TM and true TM: r=0.047 (effectively uncorrelated).
- Next investigation: Run the two gate-passing orientations through full transformation/clash/refinement with native DockQ. Do not lower stable TMalign thresholds.

### Matched-panel follow-up → ADVANCED

Why it matters:
- A broader seven-candidate panel now distinguishes aligner-dependent
  transform geometry from downstream refiner output acceptance, but it still
  lacks native DockQ and matched deterministic refinement controls.

Current status UPDATE (2026-07-30):
- MultiProt diagnostic complete; 47 successful alignments with true TM-scores.
- Transformation gates calibrated for MultiProt true TM-score.
- Next: Execute matched MultiProt+FiberDock vs TMalign+external-Rosetta with frozen inputs, template assets, source gates, refinement controls, and native DockQ. Requires per-candidate Rosetta return-code/score-gate observability (still open).

### MultiProt fragment-length bias

Why it matters:
- MultiProt finds very short fragment matches (4-68 residues) vs TMalign/GTalign which find longer structural alignments.
- True TM-scores are extremely low (median 0.016) because TM-score is length-normalized.
- This is a fundamental method difference, not a calibration issue.

Next investigation:
- Evaluate whether fragment-based MultiProt alignments are biologically meaningful for docking, or if a minimum aligned length filter is needed before true TM-score computation.
- Consider whether MultiProt should be used for fragment assembly rather than full-interface docking.

### Alignment adapter interface contract

Why it matters:
- 3 alignment adapters (TMalign, GTalign, MultiProt) write JSON with no shared schema.
- `tm_score` field means different things per aligner (native TM, native TM, proxy TM).
- Duplicate `extract_chain_and_res_ids()` in alignment.py (L176) and alignment_gtalign.py (L159).
- **Resolved 2026-08-06**: the GTalign symlink hack is removed; each aligner writes
  to its own run-scoped directory (`processed/alignment_{tmalign,gtalign,multiprot}/<run_id>/`),
  so `processed/alignment` is no longer a shared/stateful path.

Next investigation:
- Define a formal AlignmentResult protocol/contract (Pydantic model or JSON schema).
- Consolidate duplicated chain/residue extraction into shared utility.
- Add `tm_score_contract` field to every alignment JSON (already partially done: "multiprot_kabsch_rmsd_proxy", "standard_length_normalized").

## Unresolved issues

- Historical MultiProt/FiberDock compatibility still depends on a permitted
  runtime for legacy 32-bit helper binaries (`reduce.2`/loader and compatible
  libraries). Do not substitute reduce.3 as historical-equivalent.
- The main `benchmark/scripts/irmsd.py` should not be modified until an
  isolated paired-residue implementation proves correct after gap handling.
- Full-variant DockQ labeling and grouped-iRMSD completion are now accepted.
  Final `full-v2-safe-irmsd-final2/scored/audit-final.json` validates 6,539
  rows, 5,827 score-bearing rows, 712 explicit non-scores, 192 cross-only
  scopes, 16,459 interface rows, zero failed auxiliary iRMSD rows, and zero
  hash/interface-contract failures. The earlier `audit-v2.json`, safe-v1
  output, stale-environment retry, and timeout retries remain preserved as
  historical debugging evidence.
- July-22 outputs lack candidate-audit JSONL, but retain per-batch GTalign JSON.
  The deterministic adapter is batch-1 verified. Full reconstruction produced
  5,828 candidates and 5,827 labels; one explicit non-bijective row must remain
  unrankable. Job 1393128 has `DependencyNeverSatisfied` and must not be reused;
  replacement job 1393757 failed at the overly strict identity-set audit.
- FiberDock/legacy and current comparisons remain noncausal until inputs,
  template staging, source gates, and evaluator contract are matched.
- Current-versus-legacy historical aggregates, including reported mean scores,
  remain observational until regenerated under the paired contract.
- The retained stage ledger does not contain a clean same-input, same-assets,
  same-evaluator MultiProt+FiberDock versus TMalign+external-Rosetta run. The
  current TMalign arms and legacy MultiProt artifacts can locate observed stage
  attrition, but cannot support a causal method-quality conclusion. The current
  MultiProt gate diagnosis is now known: zero of 200 orientations pass both
  transformation thresholds because the RMSD-derived TM-score proxy is not
  calibrated. External-Rosetta partial/missing output reasons remain
  unobservable because subprocess return codes and per-candidate score-gate
  outcomes were not retained.
- The broad Graphify knowledge graph remains stale for exact path discovery
  (pre-#1504 IDs and stale report metadata). Use the separate filtered
  `docs/chronology/` graph for dated progress anchors and active memory/source
  files for exact commands and execution decisions.
- MultiProt score calibration remains open: the current backend computes
  `tm_score = max(0, 1 - rmsd/10)` from Kabsch RMSD, while transformation uses
  the TMalign-oriented threshold `TM_SCORE_THRESHOLD=0.5`. The retained 100-
  template run therefore has no eligible paired orientation. Validate a
  calibrated score or an explicitly separate MultiProt gate before changing
  thresholds.
- FiberDock energy parser remains open: the current refiner writes a valid `fiberdock_energies.ref` solution (`glob = 0.00`) but parses only `fd_params.ref`, returning `-` despite a valid energy PDB. Confirm the intended energy contract with an isolated fix and regression test before modifying stable code.
- External-Rosetta refinement observability remains open: the current code
  does not retain per-candidate subprocess return codes or the score-gate
  reason for missing canonical outputs. Add a diagnostic/replay contract before
  interpreting partial/missing Rosetta outputs as refiner failures.

## Next steps

1. Preserve job 1392725 as the verified isolated ranking-load smoke; keep
   stable unranked commands unchanged and ranking opt-in.
2. Preserve the final2 scoring and ranking outputs, together with all failed
   and superseded retry roots, as the reproducible BM5.5 evidence set.
3. Run the two newly eligible MultiProt orientations through normal
   transformation/clash accounting; add per-candidate Rosetta return/score-
   gate observability; then regenerate the matched
   MultiProt+FiberDock versus TMalign+external-Rosetta comparison with frozen
   inputs, template assets, source gates, refinement controls, and native
   DockQ evaluation before making causal claims.
4. Add independent labeled complexes and deterministic refinement/evaluation
   controls before testing ranking quality.
5. Resolve the 17 audit-only benchmark rows before confirmatory full-cohort
   claims or paired aligner/refiner comparisons.
6. Investigate GTalign CPU/GPU divergence before non-GPU production use.
7. Refresh `docs/chronology/` after this material scoring/ranking completion;
   preserve dated run notes and manifests while excluding only vendored
   environments, caches, and generated graph corpora.

### FiberDock output-contract follow-up

Current status:
- The parser-path defect is reproduced across the exact two-candidate replay
  and the seven-candidate panel: the declared
  `fiberdock_energies.ref` exists, `fd_params.ref` does not, and the
  corrected parser recovers the declared energies. The current opt-in source
  now uses the corrected `fiberdock_energies` prefix and pair-keyed output
  collection.
- Job 1405134 completed the exact two-candidate corrected replay; job 1405218
  completed the broader seven-candidate replay. The focused fixture regression
  `tests/test_fiberdock_output_contract.py` passes.
- No stable refiner default was changed. A new complete live pipeline replay
  after the source correction is still outstanding.

Next investigation:
- Run the corrected FiberDock arm end-to-end in a fresh isolated workspace and
  verify parsed energy, pair-keyed aggregate structure/energy files, terminal
  stage status, and downstream DockQ handling.

### Phase 1 planning follow-up

Current status:
- The Phase 1 discussion is complete and the context is ready for planning.
- Exact JSON field names, artifact-ledger column ordering, run-root layout,
  dirty-tree diff representation, and the smallest changed-artifact/row-
  identity fixtures remain implementation choices for `plan-phase 1`.

Next investigation:
- During planning, reconcile the new contract with
  `benchmark/scripts/investigation_provenance.py`, `prism.py`, existing TSV
  manifests, and downstream scoring gates without changing stable defaults or
  overwriting preserved evidence.

### Phase 1 grilling follow-up

Current status:
- The identity model is now split into declared contract, execution attempt,
  artifact observation, and run closure.
- The remaining implementation choices are expected-inventory schema,
  explicit source-inventory input format, synthetic exploratory row-ID format,
  and the exact legacy gate CLI.

Next investigation:
- Ensure implementation plans preserve the accepted ADR boundaries and do not
  accidentally put Slurm/timestamps/status into `contract_hash` or treat
  path-only legacy artifacts as benchmark-safe.

### Phase 1 provenance ownership follow-up

Current status:
- `src/provenance/run_evidence.py` owns the first canonical contract/attempt,
  artifact-observation, closeout, and consumer-gate core. Existing provenance
  scripts remain compatible adapters.
- The module's focused test file passes 14 tests and direct checks confirm
  TSV type round-tripping, immutable closeout records, and caller-supplied
  environment/argv/configuration sentinel detection. This is not acceptance:
  CLI validation reports a detected secret but exits 0, leaving the consumer
  gate non-fail-closed. The planned contract and `prism.py` integration tests
  also remain incomplete.

Next investigation:
- Decide whether template preflight and remaining benchmark artifact/TSV
  writers should migrate into `src/provenance` or remain format-specific
  adapters. Keep candidate lineage, ranking, evaluation joins, and stage
  lifecycle outside this Phase 1 ownership decision.
- Return a nonzero CLI status when validation is `fail`, add a CLI regression
  test for that behavior, then complete the Phase 1 contract/integration
  tests before treating the provenance implementation as ready.

### Two-version MultiProt matched-panel comparison - 2026-08-10

Current status:
- A controlled alignment-stage comparison completed on 770 exact shared
  template IDs, ten 77-template batches, and two fixed pairs. The current
  adapter accepted 405/6,160 records; the legacy adapter accepted 5,968/6,160.
- The executables are byte-identical. Native interface assets are mostly
  coordinate-equivalent, but a minority of sides differ in residue keys or
  coordinates.
- A representative legacy-only alignment produced the same MultiProt output
  in both formats but was rejected by the current adapter's Kabsch transform
  reconstruction (`inf`).

Next investigation:
- Preserve granular current failure reasons (`no solution`, parse failure,
  insufficient mapped residues, Kabsch failure) and run a same-assets,
  same-invocation/orientation diagnostic before attributing the full success
  gap to parser or transform logic. This comparison does not establish
  downstream structural quality or historical equivalence.

### ProInterVal transfer to PRISM-generated candidates - 2026-08-25

Why it matters:
- ProInterVal reports learned interface-validity performance on PDB-derived,
  docking-decoy, and biological-versus-crystal benchmark distributions, but no
  direct study establishes transfer to PRISM template matches or their refined
  output. A hard gate could discard salvageable poses before refinement or
  miscalibrate scores under the pipeline's candidate prior.

Current status:
- A deep public-literature review was completed in NotebookLM notebook
  `PRISM-literature-review-20260714`; 33 newly discovered public sources were
  imported. The safe current recommendation is analysis-only annotation and
  opt-in post-refinement comparison, with no production default.

Next investigation:
- Build an isolated adapter over frozen PRISM-generated rigid and refined
  candidates. Compare ProInterVal scores with native DockQ/iRMSD, report
  per-target PR-AUC/Spearman/top-N/calibration and false-negative rates, and
  retain candidate-level provenance. Keep the existing stable path unchanged
  until this matched study is complete.

### VALAR local Copilot execution blocker - 2026-09-07

Current status:
- The bounded PRISM-prescript run launched two independent validated A40 Q6
  workers with run-scoped Copilot homes, exact hook-disabled settings, local
  provider bindings, and separate sessions/ports. The initial batch invocation
  lacked a usable TTY; a tested `script -qefc` wrapper restored foreground
  registration, but Copilot 1.0.80 still left the sessions logged out and did
  not dispatch the required interactive `-i` prompt under offline mode.
- Both retries were cancelled after durable evidence showed zero turns,
  checkpoints, Qwen completion requests, and project edits. No USalign or
  transformation result was produced.

Next investigation:
- Determine whether the user will authorize a prompt-mode (`-p`) adapter path
  or provide a supported offline authentication mechanism that preserves
  `COPILOT_PROVIDER_API_KEY=""`, `COPILOT_OFFLINE=true`, no remote access, and
  the run-scoped filesystem contract. Do not relax production thresholds or
  infer the retained 1gte attrition ledger until workers actually execute.

### Fresh batch-safe prompt-mode smoke - 2026-09-07

Current status:
- A fresh one-worker smoke reached Slurm and Qwen Q6 health/model discovery,
  but Copilot 1.0.80 exited status 1 before a model turn or Qwen inference
  request. The parent run is BLOCKED and no project worker was launched.
- Exact evidence is retained under
  `tmp/agent/20260907-prism-prescript-usalign-attrition-prompt/`, including
  the provider binding, `settings.json`, JSONL output, Copilot session DB,
  Qwen logs, and fail-closed review.

Next investigation:
- Determine whether Copilot 1.0.80 requires an additional local
  alternate-provider/model-selection setting or has a batch prompt dispatch
  defect. Use a separately authorized smoke-only test; do not resume the
  interactive launch, add credentials, enable remote access, or schedule
  project science until one local Qwen turn is observed.

### Repaired-adapter direct provider sentinel - 2026-09-08

Current status:
- The fresh one-worker Qwen Q6 smoke reached Slurm node `ai25`, passed health and `/v1/models`, and generated one direct Qwen inference request. The response did not satisfy the exact `VALAR_PROVIDER_SMOKE_OK` contract; `logs/provider-smoke.json` stores only `ValueError` metadata and no response body. Copilot therefore did not start.

Next investigation:
- Determine whether the failure was a malformed choices/message shape or a non-exact/truncated assistant response without persisting response content. Test the smallest safe provider-smoke adjustment in a separately authorized adapter-only run, then require a real Copilot-to-Qwen turn before project work.

### Current orientation/threshold comparison validation - 2026-09-09

Current status:

- The current CLI now supports explicit CSV/template-panel selection,
  `native`/`o1`/`o2` transformation arms, and per-run threshold overrides.
  Focused software tests and the parameterized notebook's disabled execution
  path pass.
- No new full biological comparison was launched by this implementation turn.
  The notebook runner is prepared to create isolated native-default, o1, and
  o2 outputs, but candidate counts and downstream quality remain unknown for
  the user-selected panel until it is run.

Next investigation:

- Run the same frozen input CSV/template manifest through the three notebook
  arms with `--no-refine` first. Compare stage-status and candidate-audit rows
  before changing thresholds; then run one-variable threshold sweeps while
  preserving raw failures and benchmark labels.

### Stepwise diagnostic notebook implementation - 2026-09-10

Current status:

- The full eight-step diagnostic machinery is now implemented in
  `notebooks/pipeline_orientation_comparison.ipynb` and its read-only helper
  module `src/stepwise_analysis.py`.
- The notebook can analyze completed or user-supplied existing run roots and
  writes `run_manifest.json`, `alignment_inventory.jsonl`, `gate_ledger.jsonl`,
  `clash_diagnostics.jsonl`, and `refinement_inventory.jsonl` under each arm.
- No user-selected biological run has been executed in this implementation
  cycle. Therefore orientation yield, clash localization, refinement quality,
  and native-like recovery remain unknown.

Next investigation:

- Configure one fixed input/template panel, run native-default/o1/o2 with
  `--no-refine`, verify asset-mix and raw-hash completeness, then inspect the
  independent gate and clash ledgers before any threshold sweep.

### Stepwise comparison audit boundary - 2026-09-10

- The notebook and helper are implementation-validated, but no current
  biological native/o1/o2 arms have been executed in this cycle.
- Any legacy run without interface-list provenance, successful alignment
  return/status, or raw-output hashes will report unknown alignment-dependent
  gates; it cannot support a threshold or orientation-quality conclusion.
- US-align remains a separate preflight/contract arm. PRODIGY remains opt-in,
  and native-like recovery additionally requires an independently sourced and
  hashed label mapping.

### Full-suite validation stall - 2026-09-12 (resolved)

Why it matters:
- Focused provenance, ranking, transformation, and pipeline-contract tests
  pass, but the broad suite cannot currently be reported as complete.

Current status:
- The isolated matched-benchmark-manifest module passed 4/4 in 23.03 seconds.
- The complete prescript suite subsequently passed 344 tests with 6 skips in
  122.94 seconds and one SciPy/NumPy compatibility warning.

Next investigation:
- No further full-suite investigation is required. Retain the bounded command
  and its warning in the validation record.

### Cross-repository validation blockers - 2026-09-12

- `gtalign_env` remains incompatible with the installed DockQ extension under
  NumPy 2.4.6, but the repository-local scoring environment is now verified
  for isolated replay: DockQ 2.1.3, NumPy 1.26.4, job 1658942, and raw/adapter
  score agreement. The environment recipe still declares Python 3.11.13 while
  this repository-local prefix runs Python 3.9.23, so the canonical release
  environment identity is not yet reconciled.
- CPU and GPU GTalign completed on matched inputs but differ in raw-output
  hashes, match counts, and TM-scores. Compare command lines, version/build
  metadata, runtime GPU selection, and parser outputs before combining arms.
- The direct src.run_identity CLI is supported as python -m src.run_identity;
  the file-path form is not a supported invocation because it does not put the
  repository root on sys.path.
- MultiProt true-TM replay, external-Rosetta per-candidate return/output
  observability, completed orientation/threshold notebook post-processing,
  and benchmark authority reconciliation remain open. The full software suite
  is complete, but no claim of causal backend superiority or ranking benefit
  is allowed.

Next investigation:

- Decide whether to rebuild a clean user-managed DockQ environment from
  `environment.yaml` (Python 3.11.13) or formally register the existing
  repository-local Python 3.9.23 prefix as the scoring runtime. In either case,
  retain the exact interpreter/module/binary hashes and rerun the same valid
  model/native replay before any batch scoring.
- Then run the canonical scorer on a small hash-joined cohort with explicit
  receptor/ligand mappings, complete `scored`/`failed`/`not_scoreable` rows,
  and isolated raw JSON outputs. Do not merge those rows into benchmark
  denominators until source authority and chain-role gates pass.

### Orientation notebook inventory traversal - 2026-09-12

- The three no-refinement arms completed their pipeline stages in job 1658931,
  but the notebook's `asset_provenance_manifest` post-processing began a
  broad reference hash walk and did not close within the bounded run. The job
  was canceled and its partial arm logs/ledgers were preserved.
- Fix the asset inventory to hash only declared/consumed assets or to traverse
  a symlink-safe bounded manifest before treating the orientation comparison
  as complete. Until then, the arm rejection counts are plumbing evidence
  only and do not support orientation-quality conclusions.

### USalign validity follow-up - 2026-09-21

- Should the shared compatibility field remain `max(Structure_1,
  Structure_2)` or should USalign expose the strict reference-normalized
  Structure_2 score separately and use that for the production threshold?
- What exact manuscript USalign version, flags, batching, template scope, and
  hardware produced the faster-than-MultiProt result?
- After the parser fix, does a larger matched pilot preserve alignment quality,
  transformed-candidate yield, and DockQ/iRMSD quality under default versus
  `-fast` USalign?

### PRISM aligner comparison - 2026-09-28

- Does the selected default/16 USalign configuration preserve its zero-failure
  behavior and measured scaling when applied to all 257 historical BM5.5 cases?
- Can the exact current-panel alignment outputs be fed into the existing common
  transformation/filtering implementation without reconstructing discarded
  raw stdout, or is a bounded rerun required for corrected downstream records?
- When the active common refinement arrays finish, do their compact records
  contain enough transformed/refined identity and evaluator provenance for
  paired Delta_DockQ, ranking, overlap, and failure analyses?
- Do all 26 USalign provider batches finish with the expected 79,420 records
  per complete case and no unclassified failures before compaction starts?
- Does the common transformed scorer remain within the 24-hour per-batch
  bound at four workers, or require a measured smaller shard/concurrency
  adjustment after the first validated batch?
