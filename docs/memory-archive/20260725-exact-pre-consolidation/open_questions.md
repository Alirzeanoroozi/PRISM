# Organizer Planning Addendum

## Question

What is the fairest comparison protocol for current PRISM vs old PRISM on `T_rigid`, `T_medium`, and `T_hard`?

Why it matters:
- Without a stable protocol, conclusions about regression, improvement, or missing capability will remain noisy.

Current hypothesis:
- The comparison should explicitly control benchmark inputs, matching logic, scoring outputs, aggregation rules, and failure accounting.

Next investigation:
- Write the comparison matrix before expanding implementation work elsewhere.

## Question

How should `iRMSD` and `DockQ` be reported for single-chain and multichain contexts across the old and current paths?

Why it matters:
- The scoring layer has to support fair comparison rather than mixing incompatible contexts.

Current hypothesis:
- The benchmark flow needs clearer separation of contexts and clearer output conventions.

Next investigation:
- Redesign the scoring/reporting layer after the comparison protocol makes the critical contexts explicit.

## Question

Which benchmark slices should be used first to expose the most important current-vs-old differences?

Why it matters:
- A smaller, deliberate initial slice can surface the important gaps faster than a broad undifferentiated sweep.

Current hypothesis:
- Representative slices from `T_rigid`, `T_medium`, and `T_hard` can establish the main pattern before full expansion.

Next investigation:
- Select the first comparison subset and make the reasons for that subset explicit in the benchmark notes.

# Old Questions - form PRISM 

- Can a stricter rigid-core extraction step, such as RANSAC or iterative pruning, turn SoftAlign's correspondences into a PRISM-usable transform?
- If SoftAlign remains in the pipeline, should it be used as a prefilter only, rather than as the final transform source?
- Do SoftAlign-derived scores need separate acceptance thresholds from TMalign/GTalign because its match-count behavior and TM-like score scale differently on PRISM pairs?

# Active Questions

- Should the root `templates/` directory in `PRISM-prescript` be aliased to `new_template/template/` to reduce confusion for future compat-root runs?
- Should an isolated `benchmark/scripts/irmsd_paired_testcase.py` be implemented and validated before any change to the main `benchmark/scripts/irmsd.py`?
- What is the canonical scoring environment path across `PRISM`, `PRISM-old`, and `PRISM-prescript`?
- For multichain Fiberdock-like benchmark cases, should future storage keep grouped receptor/ligand `iRMSD`, total multichain `DockQ`, or both?
- ~~Do the current GTalign outputs need threshold tuning or a different filter path, given that the smoke test completed end-to-end but still produced `0` passed pairs on the current one-pair input?~~
  Resolution (2026-07-19): The 0 transformations with the expanded 19,855-template
  panel are explained by low structural similarity of most added templates. TMalign
  confirms this: only ~5.5% of 714,780 hits have TM-score ≥ 0.4. The `--pre-score`
  double-filtering has been fixed (`PRISM_GTALIGN_PRE_SCORE` env var), but the
  fundamental hit rate remains low for these query pairs against a general PDB
  template library. This is expected scientific behavior, not a bug.
- Should the NACCESS reinstall step be automated in the GTalign smoke-test harness, or kept as a manual prep step for isolated temp copies?
- AlphaFold 3 model parameter/weights location on VALAR: is there a centrally staged path (similar to `/datasets/alphafold3`), or must users download and stage weights themselves and point the runner via `PRISM_AF3_MODEL_DIR` / `--af3-model-root`?
- Safe cleanup: once AF3 is switched to `/datasets/alphafold3`, can the redundant `/home/rshadi25/public_databases` tree be removed to reclaim space, and are there any other jobs depending on it?
- Should future benchmark comparisons involving PRISM-main SoftAlign use their own pass/fail thresholds until a rigid-core extraction rule is defined?
- Should PRISM-prescript document the sibling ASA-replacement benchmark separately so future benchmark workflows can reuse `FreeSASA` or `RustSASA` instead of `Naccess` when appropriate?
- Should the current Python 3 downloader gain a BeEM-backed mmCIF fallback, and how should BeEM bundle chain-ID mappings be propagated into five-character target identities?

**GTalign CPU/GPU discrepancy** (2026-07-20):
- `gtalign_cpu` with `--pre-score=0.4` produces 27 vs 19,561 alignment JSONs from `gtalign_gpu` with identical parameters
- Needs investigation: is this a parameter interpretation difference or a bug in the CPU binary?
- Blocks reliable CPU-based pipeline runs for non-GPU partitions

## Question

Can CPU-based aligners (MultiProt/TMalign) with 20K templates finish a 26-batch BM5.5 run in practical time?

Why it matters:
- MultiProt CPU reached ~2% of 714K alignments in 12 hours
- TMalign CPU reached ~1% in 5 hours
- At this rate, full 26-batch runs would take weeks

Current hypothesis:
- GTalign GPU is the only practical alignment backend for 20K-template pipelines
- CPU aligners may need template subset filtering or GPU porting

## Question

Do the batch_0001 smoke-test DockQ scores (~0.78-0.84 mean) generalize to the full 26-batch results?

Why it matters:
- Only 10/257 pairs were tested in the smoke. Full-benchmark results are pending.

## Risks or Blockers

## Question

Can the exact BM3.0 88-case paper cohort be acquired as an authoritative raw-selector and native-role list with hashes?

Why it matters:
- The repository BM5/5.5 table cannot reconstruct the paper cohort; filtering by version yields 109 rows rather than 88.

Current hypothesis:
- The local paper/reference Markdown identifies the cohort and template regimes but does not provide the complete row-level source ledger.

Next investigation:
- Locate the original Benchmark 3.0 release or paper supplement, freeze all role files, and validate the 88-row manifest before exact-paper execution.

## Question

Can `working_version/multiprot` be run through a non-destructive local compatibility adapter?

Why it matters:
- The checked-in arm imports database modules, calls HTML/mail writers, cleans job directories, requires Python 2, and lacks compatible legacy template/assets.

Current hypothesis:
- A controller-boundary adapter plus a complete legacy-format template/tool staging root may make a one-pair/one-template smoke possible.

Next investigation:
- Stage one pair and one legacy-format template in a fresh directory, inject the existing local MySQL shim, disable writers/cleanup, and record executable/asset failures without modifying the pristine arm.

## Question

What fraction of all 257 BM5/5.5 rows pass strict source and exact-chain validation?

Why it matters:
- The final ten-row KUTEM smoke validates task isolation and provenance, not complete benchmark coverage.

Current hypothesis:
- Curated archive native roles will pass more rows than downloaded/full-PDB native caches, while some qualified local inputs need explicit chain materialization.

Next investigation:
- Execute the full 257-row source ledger in isolated KUTEM batches, aggregate `exit.json` and structure-validation statuses, and retain every unresolved row.

## Question

Which controlled candidate/filter/refiner arm changes downstream all-pair utility after provenance and evaluation are frozen?

Why it matters:
- The existing current-vs-legacy outputs are unpaired and cannot rank TM-align, MultiProt, Rosetta, or FiberDock causally.

Current hypothesis:
- Unequal template staging and measurement drift explain much of the July discrepancy; aligner/refiner effects remain unresolved.

Next investigation:
- Freeze the 946-template diagnostic panel, replay raw candidates through common/native filters, then run the byte-identical refiner crossover before the 255-pair confirmatory benchmark.

## Question

Should the standard score CSV expose per-interface rows directly, or should all interface metrics remain in the separate `scores_interfaces.tsv` artifact?

Why it matters:
- Existing downstream reports consume one row per model, while multichain evaluation requires one row per component interface without overwriting global values.

Current hypothesis:
- Keep model-level CSV compatibility and use `write_dockq_tsv()` for the explicit global/interface artifacts.

Next investigation:
- Integrate the separate score artifacts into the comparison batch collector and validate that every model row points to its raw DockQ JSON hash.

- `benchmark/scripts/irmsd.py` may be fragile when aligned residue lists do not stay index-compatible after gap handling.
- A safer future fix is to test an isolated variant that uses explicit paired residue lists instead of index-coupled residue lists.
- The cross-pipeline comparison may need a dedicated benchmark-mode reproduction harness if the local-style driver continues to diverge from historical PRISM-old behavior.

## Follow-up Items

- Do not change the main `benchmark/scripts/irmsd.py` until an isolated paired-residue variant is tested first.
- The TMalign JSON handoff mismatch and the Rosetta contact-helper signature mismatch are now fixed in `PRISM-prescript`; future current-pipeline smoke failures should be interpreted after those fixes, not as evidence that the old mismatch remains.
- The rigid-positive preprocessing question is resolved for `1rghB + 1a19B / 1b27AD`: chain-qualified files plus the `5.0` scaffold recover a final model; the otherwise matched `1.4` run produces zero passed pairs.

## TM-align Biological Ranking Open Questions

- Which additional independent current TM-align + Rosetta complexes should be staged next to enlarge the five-complex grouped table and support sequence-clustered held-out validation?
- Should current target handling be extended for multi-chain identifiers such as `1V8Z_AB`, or remain single-chain for the first biological-ranking benchmark?
- KUTEM job `1343750` is recorded as cancelled; job `1343808` remains pending by resources/priority and still needs final log collection if it eventually runs.
- Should the relaxed clash diagnostic (`PRISM_MAX_CLASHING_COUNT=100`) remain a rescue experiment, or should a biologically calibrated interface-clash metric replace the raw CA-clash count?
- The first current-pipeline positive rescue is now verified for `3lqcAB -> 1BPB/3K77` (DockQ `0.770`, iRMSD `0.922`), but more independent complexes are still needed before any production learned-ranking claim.
- Can the current template inventory be expanded with independently sourced interface templates like `3lqcAB` without introducing native-template leakage or uncontrolled template-family correlation?
- The tabular reranker has been evaluated on five grouped native complexes and did not improve native-like top-1 success; a larger labeled table is required before revisiting contact-GNN training.
- Determine whether intermittent Slurm failures on `ai01` are node routing, environment isolation, or cluster policy, and whether Codex should remain local/login-side with compute commands dispatched remotely.

## PRISM investigation open questions (2026-07-13)

- For the 17 final source-gate validation failures, which authority controls receptor/ligand orientation: the CSV `Complex` chain order, the raw selector roles, or the curated archive `r/l` labels? The evidence proves disagreement but not the intended correction.
- Should the 17 rows be repaired by an authoritative benchmark mapping, excluded from the primary denominator, or retained only in an audit cohort? No local full-PDB substitution is acceptable.
- MultiProt standalone compatibility is validated with Python 2.7.15 and two KUTEM binary smokes; the project-local NumPy dependency gate and separate native-tool gate now pass for the explicit current-NACCESS compatibility profile. The historical NACCESS `accall` and 32-bit helper path remain blocked.
- A project-local conda attempt for `python=2.7.15 numpy=1.16.6 pymysql=0.9.3` failed during solving, but a Python 2-compatible NumPy 1.16.6 wheel was subsequently staged into the derived project-local site and imported successfully.
- The exact PRISM paper BM3 88-case source list remains unavailable; the repository’s 257 BM5/5.5 rows must not be presented as a paper reproduction.

## Legacy toolchain status (2026-07-14)

- A positive diagnostic fixture is now identified: `T_Rigid.csv:47`, `1RGH_B + 1A19_B`, template `1b27AD` A/D. At the reference
  50% MultiProt filter it yields zero candidates; at a recorded 40% diagnostic threshold it yields one candidate per adapter.
- Historical `libgfortran.so.3` is now supplied in isolation and the historical NACCESS probe passes. The remaining unresolved
  runtime question is whether a permitted compatible 32-bit loader/`libstdc++.so.5` environment exists for the bundled FiberDock
  helpers; no substitution is currently allowed.
- The exact BM3 88-case source list and the 17 BM5/5.5 chain-contract authority decisions remain unresolved. The frozen policy is
  `source-gate-policy/v1`, with confirmatory execution blocked and 240 rows retained as strict-eligible audit candidates.

## Implementation follow-up (2026-07-14)

- Can a permitted container or host compatibility root provide the bundled FiberDock `reduce.2.21.030604` runtime, including its
  32-bit loader and `libstdc++.so.5`, without modifying the executable or mixing tool versions?
- Does reduce.3 produce equivalent hydrogenation, FiberDock energy, ranking, and evaluator behavior to the historical reduce.2 helper?
  A one-pair exploratory run completed only after explicit substitution and parent-directory initialization; it cannot answer equivalence.
- On the matched 12-row pilot, does GTalign meet the pre-registered margins (at most 0.02 loss in best GlobalDockQ@20, 5
  percentage-point success loss, and 0.5 Å grouped-iRMSD loss) relative to TM-align?
- Which authoritative source resolves the 17 curated chain-contract disagreements before the 257-row confirmatory denominator is
  opened?

## Output validity and repair options (2026-07-14)

- Where is the duplicated-chain output introduced: MultiProt staging, transformation, FiberDock/Rosetta writing, or the historical chain-fix postprocessor?
- Can the unavailable `rosetta_output_1_chainfixed` artifact be recovered from the prior benchmark workspace, or must both arms be regenerated with an explicit chain-preserving writer?
- Which repair is scientifically acceptable: native 32-bit reduce.2 runtime, native/source-built Reduce replacement, containerized legacy runtime, or exploratory reduce.3 substitution? Each requires matched hydrogenation, energy/ranking, and scoreability tests before primary use.
- Does FiberDock itself preserve partner identities when given chain-distinct inputs, or does its output writer require a postprocessing adapter? A repair must be tested on a two-chain fixture and the full positive diagnostic before benchmark arrays.
- The boundary-split probe is scoreable but not provenance-preserving by proof. Can a direct chain-distinct-input run reproduce its chain lengths, scores, and interface assignment without postprocessing?

## Final replay follow-up (2026-07-14)

- Can the 143 current retained models be regenerated or recovered with complete residue correspondence so the strict no-align evaluator yields valid structural metrics?
- Can the missing historical chain-fixed artifact be recovered, allowing the previous current/legacy report to be audited against the exact retained PDB rather than only its CSV row?
- Which candidate-producing distinct-chain fixture reaches FiberDock and permits a direct chain-preservation test? The current probe reached refinement with zero candidates and is inconclusive.
- After evaluator and source gates pass, do matched pipeline-generation runs preserve the observed legacy values and determine whether any remaining difference is from alignment, filtering, transformation, refinement, or ranking?

## PyRosetta and GTalign follow-up (2026-07-14)

- Can an authorized PyRosetta cp311 wheel be staged from a network-capable login host, and does `pyrosetta.init()` succeed without conflicting with the existing external-Rosetta environment?
- On identical canonical poses, how do PyRosetta and external Rosetta differ in chain preservation, total score, runtime, DockQ, iRMSD, interface size, and prediction/failure counts?
- Does GTalign CPU preserve transformation and refinement behavior on matched PRISM inputs, and does a Tesla V100 GPU run improve throughput enough to justify a separate GPU arm?
- Can a downstream paired GTalign/TM-align benchmark be completed under the source/evaluator gates so DockQ and iRMSD differences are measured rather than inferred from TM-score?

## PyRosetta installation update (2026-07-15)

- Resolved: PyRosetta is installed and initializes in `gtalign_env`; the one-pose adapter smoke succeeds.
- Remaining: run matched PyRosetta versus external-Rosetta benchmark shards and verify chain preservation, DockQ, iRMSD, interface size, runtime, prediction counts, and failure stages.
- Freeze PyRosetta random-seed options for benchmark arrays and decide whether to normalize temporary-path comments before comparing output hashes.

## Reusable setup follow-up (2026-07-15)

- The setup guide now freezes the standard array seed.  Should a structural-content hash (excluding PyRosetta energy-table path comments) be added to the adapter before benchmark comparison?
- After authoritative resolution of the 17 source-gate rows, which minimal matched shard should be used first for external-Rosetta versus PyRosetta and TM-align versus GTalign comparison?

## Matched pilot execution (2026-07-15)

- The planning manifest records candidate-generation, top-20-refinement, and strict-scoring task IDs, while the current batch launcher executes one 12-row arm-level pipeline. Should a follow-up adapter split and score those stage records per row, or should the pilot contract be revised to an arm-level execution unit before expansion?
- The submitted pilot arrays `1358990` and `1358991` are still in structural alignment; no completed model/score denominator is available yet.

## Template and PyRosetta diagnosis (2026-07-15)

- Should the project add a tested old-to-new template converter that copies `.int` files to `_int.pdb`, generates `interfaces_lists`/`contacts` JSON, and records coordinate-level equivalence, or should the old schema remain legacy-only?
- Which exact `1jxqA` TM-align output causes the 508 parser `IndexError` failures? The current temporary-directory implementation removes raw `out.tm` files after each call, so a retained failing-output mode is needed before changing the parser.
- How should generated transformation L/R halves be assembled and chain-mapped into the PyRosetta manifest? The PyRosetta adapter works on a supplied assembled model, but `prism.py` does not invoke it.
## Benchmark 5.5 execution follow-up (2026-07-15)

- Will all 26 current TM-align/external-Rosetta batches complete within the 72-hour array wall time, and how many outputs pass strict residue/chain validation?
- After the production array finishes, can the collector map every model to exactly one Benchmark 5.5 pair without filename ambiguity and produce strict DockQ/iRMSD rows?
- Should a separate full PyRosetta array be launched after external-Rosetta collection, or is the validated bounded PyRosetta arm sufficient until the strict-score denominator is available?
- Can FiberDock full refinement be enabled only after a compatible native 32-bit helper/runtime is verified, without treating the reduce.3 exploratory substitution as historical-equivalent?

## Verification-gate blockers (2026-07-18)

- Which authoritative source resolves the frozen `blocked_source_authority` decision so the 240 strict rows can be opened? The source-matrix implementation is complete and deterministic, but authorization cannot be inferred from similarity eligibility.
- Can modern JSON hotspot/contact assets be reconciled semantically with the derived legacy assets? The parity audit records 19,005 disagreements and 850 missing modern profiles; no confirmatory filter denominator should mix these sources.
- Can a strict published-protocol positive full-pose canary be produced under the current both-normalized-score GTAlign contract? Normal-threshold AI/COSBI canaries currently complete with no predictions, while relaxed thresholds are diagnostic only.

## Benchmark 5.5 scoring follow-up (2026-07-16)

- Does corrected array `1361649` finish all 26 batches with explicit completed/no-prediction/failed records, and does dependent scorer `1361656` retain a per-model error row for every failed or unmatched prediction?
- The full-contract smoke has a high DockQ (`0.8913`) but high analyzer iRMSD (`29.599 Å`) for the same generated model. The numbers are valid under the requested PRISM-main scripts, but their interface-selection semantics should be reported transparently and reviewed after a multi-pair sample; do not infer model quality from either value alone.
- Should the PRISM-main analyzer's sibling-path defect be repaired upstream in PRISM-main? This repository currently avoids modifying that sibling checkout and uses a documented disposable runner layout.
- How should the canonical benchmark scorer represent multichain partner groups? It must enforce a complete chain-bijective model→native assignment (for 1AHW: L→A, H→B, A→C), report `GlobalDockQ` and/or explicitly selected cross-partner interfaces, and calculate grouped iRMSD with equal-length partner chain lists. The current analyzer's `HLA:AC`/`HLA:BC` calls are invalid for this purpose.
## Legacy CLI compatible runtime

Why it matters:
- The FiberDock/MultiProt smoke cannot produce a scientific model on the current host despite Python 2 syntax and local assets being available.

Current hypothesis:
- A runtime/container providing `libgfortran.so.3` and allowing the bundled 32-bit MultiProt/NMA binaries is required; ABI substitution is unsafe.

Next investigation:
- Identify an authorized compatible library/container or cluster node profile, then rerun the same `codex_smoke_20260717` input under Slurm and compare stage outputs.
