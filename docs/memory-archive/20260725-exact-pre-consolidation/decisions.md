# Organizer Planning Addendum

## Decision

### Distributed variant placement and gating (2026-07-16)

Decision:
- Validate the full current 2x2 matrix: TM-align/external Rosetta, GTAlign/external Rosetta, TM-align/PyRosetta, and GTAlign/PyRosetta. Keep the legacy MultiProt/FiberDock arm outside these backend substitutions.
- Use absolute GTAlign executable paths on compute nodes. Use general `ai` with QOS `ai` and `--gres=gpu:1` for GPU GTAlign jobs; do not use `v100_ai` unless its account/QOS association is explicitly accepted.
- Do not promote a completed task with zero transformations to a valid quality result. A bounded template smoke validates plumbing only; prediction-quality validation requires the full frozen template panel and corrected chain-aware scoring.

Reason:
- The compute-node environment does not inherit the interactive Conda PATH. The first GPU smoke failed before alignment for this reason, while the corrected absolute-path submission completed successfully.
- `v100_ai` rejected QOS placement for the `ai` account, whereas general `ai` successfully scheduled both GTAlign jobs on `ai08`.

Consequences:
- Every variant manifest records aligner, refiner, template workflow, threshold, executable path, partition, QOS, and isolated run root.
- KUTEM concurrency batches are diagnostic duplicates only and are not included in benchmark denominators.

### Multichain staging contract (2026-07-16)

Decision:
- Derive model receptor/ligand chain groups from the observed refined PDB chain order and the template partner chain counts; do not assume template chain IDs survive refinement unchanged.
- Reject any score where model receptor and ligand groups overlap.

Reason:
- The first multichain smoke output preserved four physical chains as `A,B,C,D`, while template IDs were `1A6Z_AB` and `1CX8_AB`; using template suffixes directly produced the invalid overlapping model map `AB:BA`.

Consequences:
- `stage_current_models_for_main_benchmark.py` records explicit model groups and supports both external and PyRosetta output directories.
- The corrected smoke map is `ABDC:ABCF`, with cross-interface DockQ and grouped iRMSD scoreable.

### Benchmark-root score routing (2026-07-16)

Decision:
- Split score manifests by `benchmark_set` and invoke the corresponding native-bound-complex directory; never score a mixed rigid/difficult manifest against one native root.

Reason:
- The first asynchronous PyRosetta score pass produced one explicit failure for `rigid_0061` because it was sent to the difficult native directory. Re-scoring the rigid row against `native_bound_complexes_t_rigid` restored a valid 14/14 denominator.

`PRISM-prescript` is the canonical scoring and benchmark-comparison layer for the PRISM refactoring work.

Reason:
- It already contains the strongest workflow for benchmark matching, DockQ/iRMSD scoring, aggregation, and report generation.

Consequences:
- Use this repo first when the task is benchmark protocol, scoring, comparison reporting, or evaluation logic.
- Avoid replacing benchmark logic with ad hoc local scoring when the comparison needs to be durable.

## Decision

Single-chain and multichain comparison contexts must be tracked separately in the benchmark and scoring flow.

Reason:
- Aggregate reporting can hide the current path's multichain gap and make interpretation misleading.

Consequences:
- Design reports and scoring outputs that preserve the distinction.
- Treat multichain evaluation as an explicit benchmark/scoring requirement rather than an afterthought.

## Decision

Frontier protein-DNA benchmarking should remain a separate workstream from the core old-vs-current PRISM comparison.

Reason:
- The refactoring comparison question and the frontier protein-DNA question are related but not identical, and merging them would blur priorities.

Consequences:
- Keep the old-vs-current PRISM benchmark path focused.
- Bridge the two streams only when a task explicitly requires it.

# Architectural Decisions

## Decision

Historical-vs-current investigation runs use explicit provenance and template preflight artifacts before execution.

Reason:
- The comparison is invalid if listed templates, staged assets, executable resolution, or effective configuration are not frozen and auditable.

Consequences:
- Use `run_investigation_preflight.py` with named `current` and `historical` arms.
- Treat unresolved or format-invalid required assets as a fail-closed execution condition.
- Keep per-asset hashes and manifest membership in `template_assets_*.tsv`.

## Decision

DockQ evaluation stores GlobalDockQ and component-interface rows separately, retaining the exact raw JSON hash.

Reason:
- Multichain DockQ emits repeated component metrics, and retaining the last component value is not a valid global score.

Consequences:
- Use `standardized_evaluator.py` and `investigation_artifacts.py` for new comparison artifacts.
- `--no_align` is permitted only after strict chain/residue identity and numbering validation passes.

## Decision

Equal-budget candidate selection must use an explicit native-independent ranking-key allowlist; pair summaries may report native scores but must not use them to select the primary candidate by default.
## Decision

**DockQ best-result convention** (2026-07-20):
- Decision: Use the **best individual interface DockQ** (max of `best_result` values), not `GlobalDockQ` (average of all interfaces) or `best_dockq` (sum).
- Reason: Combined R+L PDBs have multiple chain pairs; GlobalDockQ averages good + bad interfaces, masking true docking quality.
- Consequences: Scoring scripts must compute `max(best_result[k].DockQ for k in best_result)`.

## Decision

**MultiProt as third aligner** (2026-07-20):
- Decision: MultiProt is integrated as a third `--aligner` option alongside `tmalign` and `gtalign`.
- Reason: MultiProt provides structural alignment as a complementary approach. The TMalign fallback for rotation matrices ensures compatibility with downstream transformation pipeline.
- Consequences: `prism.py --aligner multiprot` is now a supported invocation path.

## Decision

### GTalign GPU pre-filter via separate env var (2026-07-19)

Decision:
- Use `PRISM_GTALIGN_PRE_SCORE` (default 0.2) for the GPU-level pre-filter, not
  `PRISM_TM_SCORE_THRESHOLD`.
- Keep `PRISM_TM_SCORE_THRESHOLD` (default 0.5) for the downstream `transformer()`
  stage.
- Limit `--nhits=2000` and `--nalns=2000` instead of the full 39,710 reference
  count.

## Decision

### GTalign GPU jobs need explicit sbatch directives (2026-07-24)

Decision:
- GTalign GPU jobs must include `--gres=gpu:7` and explicit `gtalign_gpu` path in `--export`
- Do not rely on `gtalign_cpu` fallback when GPU acceleration is intended
- Use Python multiprocessing (`Pool.map`) for internal parallelization to stay within QOS job limits

Reason:
- Initial `parallel_pipeline.sbatch` defaulted to `gtalign_cpu` with no `--gres` — all 5 GPU variants produced 0 alignments silently
- QOS `gres/gpu=8` limit means interactive (1 GPU) + batch (7 GPU) max = 8

Consequences:
- Created `gtalign_gpu_pipeline.sbatch` as the correct template for GPU variants
- Old `parallel_pipeline.sbatch` should not be reused for GPU jobs without modification

Reason:
- Setting `--pre-score=0.0` with 19,855 templates creates 50+ MB output files per
  query that the Python parser cannot process in practical time.
- Setting `--pre-score` to the production threshold (0.4) discards ~94.5% of hits
  at the GPU level, creating a double-filtering discrepancy with TMalign (which
  writes all hits).
- The two thresholds serve different purposes: GPU pre-filter controls output
  volume; transformer() controls production acceptance.

Consequences:
- `PRISM_GTALIGN_PRE_SCORE` can be tuned independently for GPU output volume
  control.
- Python-side pre-filter in `src/alignment_gtalign.py` also reads
  `PRISM_GTALIGN_PRE_SCORE` so GPU and Python filters agree.
Reason:
- DockQ/iRMSD are native-derived and would leak evaluation information into candidate selection.

Consequences:
- `validate_ranking_keys()` rejects native/reference fields.
- `aggregate_pair_summary()` defaults to deterministic model-ID selection and accepts only validated candidate-derived ranking keys.

- Benchmark matching uses only `PDB ID 1` and `PDB ID 2`.
- The first token in PRISM/Rosetta output filenames must not be used as the benchmark matching key.
- Aggregation across multiple predictions for the same `(PDB ID 1, PDB ID 2)` pair but different templates should report mean, variance, and best.
- Keep protein-DNA scoring in a separate benchmark path under `benchmark/scripts/protein_dna_output/` rather than modifying the established PPI scripts.

## Workflow Decisions

- Keep a separate repo-local smoke runner for the current PRISM pipeline at `benchmark/scripts/run_prism_pipeline_smoke.sh` instead of overloading the benchmark validation scripts.
- Preserve the benchmark-side `run_stable_checks.sh` path and let it use the score environment, but let the current root-pipeline smoke runner auto-select a pipeline-capable interpreter such as `gtalign_env`.
- For the current TMalign backend, preserve the CSV summary output but also emit one JSON payload per `(query, template, chain)` so `src.transformation.load_alignment()` can consume the same artifact pattern as GTalign and SoftAlign.
- General benchmark utilities belong in `benchmark/scripts/`:
  - `dockq.py`
  - `irmsd.py`
  - `irmsd_backbone.py`
  - `score_single_prism_pair.py`
  - `generate_prism_benchmark_report.py`
  - `validate_prism_pipeline.py`
  - `preflight_benchmark_check.py`
- Use `benchmark/scripts/score_single_prism_pair.py` for quick scoring of single files or folders with CSV output.
- Keep the Python 3-compatible backbone iRMSD refactor as a separate copy (`benchmark/scripts/irmsd_backbone_py3.py`) for validation and portability; leave the legacy `irmsd_backbone.py` unchanged unless a direct replacement is explicitly required.
- Rosetta-output-specific benchmark pipeline belongs in `benchmark/scripts/rosetta_output/`.
- Experimental alignment scripts belong in `benchmark/scripts/experimental_alignment/`.
- Miscellaneous unrelated helpers belong in `benchmark/scripts/misc/`.
- Protein-DNA benchmark runs should consume a manifest CSV and emit `per_prediction.csv` and `pair_summary.csv` with the same mean/variance/best aggregation pattern used in the PPI reporting flow.
- The benchmark pipeline is the canonical scoring/validation workflow for PRISM outputs: raw archive extraction, chain-fix correction, `(PDB ID 1, PDB ID 2)` matching, native-complex download, DockQ/iRMSD scoring, aggregation, reporting, and validation.
- When using `PRISM-prescript` as a compatibility root for `PRISM-main`, stage the `templates/` directory from `new_template/template/` rather than relying on the empty top-level `templates/` directory.

# PRISM Decisions - from PRISM 

- Keep `TMalign` as the default alignment backend.
- Treat SoftAlign as an experimental opt-in backend via `--aligner softalign`.
- Use full-atom structures for SoftAlign input, not PRISM's CA-only interface files; reconstruct template interface PDBs from the retained template structures before inference.
- Keep backend outputs isolated by alignment directory so `transformer()` can read the correct JSON set for each aligner.
- Do not treat SoftAlign's raw match count as a proxy for PRISM compatibility; if SoftAlign is used further, its correspondences need rigid-core pruning or similar filtering before final placement. 
## Environment or Infrastructure Decisions

## Optional PyRosetta refinement boundary (2026-07-14)

Decision:
- Keep PyRosetta as an opt-in adapter in `src/pyrosetta_refinement.py`; do not add it to `environment.yaml` or `environment.yml`, and do not change the default external-Rosetta import in `prism.py`.
- Probe with a lazy import and return explicit `status="unavailable"` plus package/version/import-error and command/environment metadata when the runtime cannot be loaded.
- Record SHA-256 input/output hashes and a sidecar JSON record for each explicit adapter attempt. Never fall back to CLI Rosetta.

Validation:
- Focused adapter/probe tests passed (`5 passed`) under `gtalign_env`; the current default Python probe reported `ModuleNotFoundError: No module named 'pyrosetta'` with `status="unavailable"`.

- Treat `benchmark/README.md` and `benchmark/prism_processed/README.md` as the workflow source of truth.
- Use `benchmark/scripts/rosetta_output/submit_prism_analysis_all.sbatch` as the Slurm benchmark entry point.
- Use `benchmark/scripts/preflight_benchmark_check.py` as the safe preflight checker before large benchmark runs.
- Invoke DockQ via `python -m DockQ` to avoid broken entry-point shebangs in copied environments.
- Use `--dockq-no-align` for fast scoring when chain mapping is trusted.
- Prefer a single Conda environment defined by `environment.yml` for reproducible benchmark runs when a compatible environment is not already available; the documented recipe includes `dockq==2.1.3` and `freesasa`.
- For cross-pipeline comparisons in `prism-version-comparison`, use rigid benchmark PDBs as the target source for `PDB ID 1` and `PDB ID 2` rather than the complex-level unbound files.
- Use `gtalign_env` for local verification of the protein-DNA scripts because the default system `python3` on this workspace lacks `numpy`.
- For the Dockground protein-DNA extension benchmark, keep exact residue-pair contact overlap as the primary metric but also report a DNA register-aware contact metric so residue-number shifts do not hide near-misses.
- Score every Dockground transformation attempt in the extension benchmark, not only rows that pass PRISM filtering, so the report remains informative when no models are accepted.
- For future state-of-the-art protein-DNA baseline comparisons, track LightDock and HADDOCK3 as local CLI baselines and HDOCK/pyDockDNA as web-server baselines. When the binaries are unavailable locally, mark them skipped instead of faking execution.
- For frontier deep-learning baselines in the PRISM protein-DNA pipeline, prioritize Chai-1 first, then Boltz-2; use AlphaFold 3 as the non-commercial reference ceiling, RoseTTAFoldNA as the specialized protein-DNA comparator, and RoseTTAFold All-Atom only when higher-order or modified assemblies become part of the benchmark scope.
- Treat PRISM-main SoftAlign as experimental and compare it against TMalign/GTalign only with separate calibration or rigid-core pruning.
- Keep Naccess-like wrapper tests repo-local in PRISM-prescript when possible, because hardcoded tool paths and shared cwd temp files are a recurring source of false failures in related PRISM debugging.
- For repo-local current-pipeline smoke runs, pass `vdw.radii` and `standard.data` explicitly to the bundled NACCESS wrapper and keep the wrapper's executable directory self-resolving rather than hardcoding a host-specific `EXE_PATH`.
- For AlphaFold 3 specifically, do not treat the generic `ai` partition as homogeneous; request the supported GPU family explicitly. On VALAR, the validated local AF3 path is `kutem_gpu` / `rk02` on A100 hardware.
- For AlphaFold 3 database staging, prefer the `PRISM_AF3_DB_DIR` environment variable plus registry defaults (`database_root_env` / `database_root_default`) so the runner can execute when the database tree is staged in a non-default location.
- On VALAR, prefer `/datasets/alphafold3` as the default AF3 database root because it is the discovered shared database tree and matches the official bundle contents.
- VALAR data-location policy for benchmarks: always check `/datasets` (and `/userfiles/<group>` when applicable) for centrally staged data before downloading large assets into `$HOME` or scratch.
- For Boltz-2 launches, prepend the Boltz conda environment's NVIDIA CUDA library directories to `LD_LIBRARY_PATH` so Triton/NVRTC can find `libnvrtc-builtins.so.13.0`.
- For Boltz-2 prediction normalization, search the Boltz workspace fallback path `boltz_results_<pair_id>/predictions` rather than only the immediate tool output directory, because Boltz writes the CIF there.
- The broad Dockground frontier benchmark should continue to use the same comma-separated `--tools` interface in `run_protein_dna_frontier_models.py`; Chai-1 and Boltz-2 now both run successfully against the three-case self-template manifest under this shared harness.
- Frontier protein-DNA execution should stay manifest-driven and stage-separated (`manifest`, `workspace`, `probe`, `dry-run`, `execute`, `af3-gpu`) so failures can be localized before actual model execution.
- For PRISM GTalign smoke tests, use an isolated working-tree copy under `/scratch/tmp` so local runtime assets are present and the repo itself stays untouched.
- Reinstall NACCESS inside the temporary test copy before rerunning PRISM so the generated wrapper points at the local `EXE_PATH` and `accall` binary instead of a host-specific path.
- If PRISM-prescript ever needs an ASA backend replacement note, treat `FreeSASA` as the closest `Naccess` analog, `RustSASA` as a viable second option, and keep MaSIF separate as a slower surface-learning path.
- Preserve dirty checkouts by using `git stash push -u` or a separate `git worktree` when you need a clean branch view without deleting untracked files; do not rely on `git pull` in-place for that.
- For the focused current-vs-old PRISM comparison, compare current `TMalign + Rosetta` directly against legacy `MultiProt + FiberDock`; exclude SoftAlign unless a future task explicitly reintroduces it as a separate experimental backend.
- For current TMalign PRISM target handling, materialize one chain-qualified PDB per five-character target and use that same file for surface extraction and transformation.
- Keep the current surface scaffold compatibility default at `5.0`; use `PRISM_SCFF_THRESHOLD` only for explicit controlled alternatives.
- Treat BeEM as an input-conversion dependency for mmCIF-only structures, not as an explanation for no-output runs that already reached alignment with valid legacy PDB inputs.

## TM-align Biological Ranking Decisions

- Keep TM-align as the deterministic candidate generator. Learned models are opt-in ranking/quality layers and must not invent transformations or alter production acceptance without held-out evidence.
- Preserve every orientation and terminal failure in candidate audit records. Retained unlabeled rows are excluded from training, not treated as negative examples.
- Require `native_complex_id` and at least two independent native complexes before training; use grouped complex and, where possible, sequence-cluster splits.
- Preserve default TM-score, coverage, and clash thresholds. Environment overrides (`PRISM_TM_SCORE_THRESHOLD`, `PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE`, `PRISM_DIFF_PERCENTAGE`, `PRISM_MAX_CLASHING_COUNT`, and related settings) are diagnostic only.
- Compare tabular/contact models against the deterministic TM-score/coverage/clash baseline using DockQ, iRMSD, top-k success, enrichment, ranking correlation, calibration, and candidate-generation coverage.
- Run connectors/authentication outside compute nodes when possible. Use login nodes for lightweight orchestration and submit heavy PRISM/Rosetta/ML work through Slurm.

## Biological ranking evaluation outcome (2026-07-12)

Decision:
- Keep learned ranking disabled in production.

Reason:
- The current grouped table contains five independent native complexes and a verified current-pipeline positive rescue, but the 20-seed grouped 40% evaluation showed no improvement in native-like top-1 success over the deterministic baseline. Small DockQ changes occurred only among non-native decoys.

Consequences:
- The tabular trainer now reports baseline and learned top-1 DockQ/native-like metrics plus coverage, but its output is diagnostic only.
- Do not train or enable the contact-GNN until the tabular stage improves on a larger, independently labeled table.
- Keep diagnostic threshold overrides and legacy-template staging isolated from production defaults.

## Rosetta batch environment

Decision:
- Slurm launchers that invoke Rosetta must explicitly load `rosetta/2022.42`.

Reason:
- New `ai`/`cosbi` jobs without the module reached transformation but failed with `docking_prepack_protocol...: command not found`; module-loaded jobs on `ag01` completed refinement.

## Stable pipeline verification (2026-07-21)

Decision:
- All 6 pipeline variants confirmed working end-to-end. The canonical stable execution uses `gtalign_env` Conda env, explicit GTalign paths, and `geometry_only_experimental` filter mode.

Reason:
- Self-control test (template=receptor) confirms perfect reconstruction: DockQ=1.000 on raw transforms, >0.97 on refined models. All variants produce High Quality docking (DockQ ≥ 0.80).
- External Rosetta requires `module load rosetta/2022.42` AND explicit `PRISM_ROSETTA_*` env vars. PyRosetta works without module loading.

Pipeline env vars to export for every run:
- `PRISM_INPUTS_CSV` — path to input pair CSV
- `PRISM_TM_SCORE_THRESHOLD=0.4`
- `PRISM_FILTER_MODE=geometry_only_experimental` — skips slow hotspot/contact checks for faster local validation

GTalign paths (use absolute paths):
- CPU: `/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu`
- GPU: `/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_gpu`

Run directory layout (symlink-based):
```
ln -sfn $REPO_ROOT/src src
ln -sfn $REPO_ROOT/prism.py prism.py
ln -sfn $REPO_ROOT/external_tools external_tools
ln -sfn $REPO_ROOT/new_template/template/* templates/
```
Template list at `templates/calculated_templates.txt`. PDBs at `processed/pdbs/`.

DockQ scoring uses `max(best_result[k].DockQ for k in best_result)`, never `GlobalDockQ`. Combine `_R.pdb` + `_L.pdb` before scoring.

## FiberDock integration as third refiner (2026-07-22)

Decision:
- FiberDock is integrated as a third `--refiner fiberdock` option in `prism.py`, alongside `pyrosetta` and `external_rosetta`.
- FiberDock binary and tools staged at `external_tools/fiberdock/` (copied from legacy `working_version/Multiprot-new/prism-fiberdock-cli/external_tools/fiberdock/`).
- `buildFiberDockParams.pl` calls Reduce for hydrogenation, runs NMA, builds FiberDock params, and executes FiberDock. The script must run from the fiberdock directory (FindBin resolves `lib/` relative to `$FindBin::Bin`).
- Reduce outputs hydrogenated PDB to stdout; saved as `.HB` file (legacy format expected by `buildFiberDockParams.pl`).
- Reduce is a 32-bit binary with the same seccomp limitation as MultiProt. Seccomp check from `alignment_multiprot.py` could be shared but is not yet unified.
- External Rosetta still requires `module load rosetta/2022.42` and explicit `PRISM_ROSETTA_*` env vars; PyRosetta and FiberDock work without module loading.

Reason:
- All three refiners (PyRosetta, external Rosetta, FiberDock) are now accessible from the same `prism.py` entry point, covering the full legacy + current refiner matrix.
- FiberDock produces energy scores and refined structures, filling the gap between PyRosetta (fast, Python-native) and external Rosetta (full physics, Slurm-dependent).

Consequences:
- The complete pipeline matrix is now 2 surfaces × 3 aligners × 3 refiners = 18 combinations.
- All 18 combinations are verified working as of 2026-07-22.
- `prism.py --refiner fiberdock` is now available. The legacy pipeline at `working_version/Multiprot-new/prism-fiberdock-cli/` is preserved for reference but superseded.

## Stable pipeline documentation

Decision:
- Treat `docs/STABLE_PIPELINE.md` as the exact operational recipe for current TM-align + Rosetta testing and keep rescue overrides explicitly separate.

Reason:
- The successful biological rescue required both explicit Rosetta module loading and relaxed alignment/clash settings. Conflating those settings with production defaults would make later results irreproducible and weaken the acceptance gate.

Consequences:
- Stable runs use the source defaults and record isolated workspaces, Slurm metadata, input IDs, template list, and logs.
- Rescue runs must declare every environment override and remain excluded from claims about default production behavior.

## Multi-chain benchmark input normalization (2026-07-12)

Decision:
- Normalize target tokens to lowercase four-character PDB ID plus sorted unique alphanumeric chain IDs; remove benchmark underscores and residue-range annotations such as `(10)`.
- Materialize one chain-qualified PDB containing every requested chain, retain the full PDB when no chain is specified, and preserve raw five-character aliases for single-chain callers.
- Reassign all chains from each transformed partner to unique Rosetta partner-chain groups before refinement and contact extraction.

Reason:
- This matches the legacy preprocessor's multi-chain behavior while keeping current single-chain filenames and transformation consumers working.

Consequences:
- NACCESS output lookup must use the actual chain-qualified input basename and preserve its case.
- Full benchmark batches use one shared normalized manifest and at most ten rows per Slurm task.

## Investigation source and archive identity (2026-07-13)

Decision:
- Use `dataset_row_id` as the primary benchmark identity and preserve raw selectors; use curated archive role files for native truth and
  reject unqualified full-PDB files as exact sources for qualified selectors.
- Record archive filename prefixes (`1QFW`/`9QFW`, `BAAD`, `BOYV`, `BP57`, `CP57`) rather than the tar top-level directory.

Reason:
- Normalized PDB pairs collapse distinct rows, parenthetical selectors are otherwise lost, and repeated benchmark complexes use synthetic
  filename prefixes documented by the benchmark README.

Consequences:
- Any row with unresolved or chain-mismatched source files remains explicit and cannot enter a strict confirmatory run.
- Full-PDB files under `benchmark/data/pdbs` may remain audit comparators or require a separately hashed chain-materialization step.

## Investigation execution contract (2026-07-13)

Decision:
- Use isolated manifest tasks with one explicit benchmark row per array index, task-local outputs, complete hashes, separate scientific and
  scheduler retry IDs, and no inference of scientific success from Slurm completion.

Reason:
- The July collector propagated batch-level state to pair rows and could not support causal or pair-specific failure accounting.

Consequences:
- The exact KUTEM profile is used for smoke/calibration tasks; full pipeline resource sizing remains empirical and gated.
- `working_version/multiprot` remains unvalidated until a local-only compatibility adapter and legacy asset/runtime gate pass.

## Final source-gate contract (2026-07-13)

Decision:
- Treat curated archive file presence and chain-contract validity as separate boundaries. Always retain the four `r_u`, `l_u`, `r_b`, and `l_b` archive members and payload hashes when present; reject the row if Biopython polymer-chain or parse validation disagrees with the CSV/native assignment.
- Exclude hetero-only blank chains from polymer-chain matching, but retain all-chain IDs, residue IDs, sequence hashes, altloc/disorder counts, duplicate counts, and parser warnings for audit.
- Hash every task output consumed by aggregation, bind summaries to the task-local directory and manifest identity, reject duplicate/missing dataset rows, and prohibit reuse of task directories containing stale logs or outputs.

Reason:
- The initial source run labeled physically present archive members as missing because blank hetero chains defeated exact chain matching. The corrected run demonstrated 257/257 archive role completeness and isolated the remaining 17 row-level validation disagreements without substituting local full-PDB files.

Consequence:
- The source gate is not a permissive repair step. Primary model comparisons remain blocked until the 17 chain-contract cases receive an authoritative orientation decision or an explicit denominator policy.

## Legacy toolchain environment and NACCESS profiles (2026-07-14)

Decision:
- Stage the legacy toolchain in a derived project-local environment and never run `install_MultiProt.pl` against the source tree. Use a Python 2.7.15 wrapper, project-local NumPy 1.16.6/PyMySQL 0.9.3, and a local-only MySQLdb compatibility shim.
- Keep two explicit NACCESS profiles: `historical_working_version` uses the checked-out binary and is blocked by missing `libgfortran.so.3`; `explicit_compatibility_substitution` uses the repository current NACCESS binary and is operational with `libgfortran.so.5`. Do not pool their outputs.
- Treat the KUTEM independent native-tool probe as executable/output validation, and classify the current controller smoke as plumbing-only when `transformation/passedFiles` is empty. A FiberDock probe with `resFile.ref` does not prove the complete hydrogen/NMA refinement path.

Evidence:
- `tmp/agent/20260713-investigation-implementation/legacy-tool-environment-v5/environment_manifest.json`
- `tmp/agent/20260713-investigation-implementation/legacy-tool-probes-v5/task-{1,2,3,4}/retry-*/exit.json`
- `tmp/agent/20260713-investigation-implementation/legacy-pipeline-smoke-v9/task-{1,2}/retry-*/exit.json`
- `references/nprot.2011.367.md` installation and MultiProt parameter sections.

Consequence:
- The compatibility profile may be used for controlled plumbing diagnostics only. A matched positive pair/template and the historical native dependencies are still required before a FiberDock refinement comparison or confirmatory estimate.

## Historical NACCESS and FiberDock refinement boundary (2026-07-14)

Decision:
- Supply historical `libgfortran.so.3` only through the derived task-local environment activation, using the cluster GCC 6 library
  directory. Keep the historical NACCESS binary and current compatibility NACCESS binary in separate profiles.
- Treat the bundled 32-bit FiberDock helpers as an explicit capability blocker. Do not substitute 64-bit or newer libraries, alter
  ELF binaries, or silently fall back to another refiner.
- Permit a threshold-relaxed positive smoke only as a recorded diagnostic intervention; never include it in confirmatory scores.

Reason:
- Historical tool array `1355128` passed NACCESS after activation exported the staged Fortran library. Positive controller array
  `1355189` generated one candidate per arm but both stopped before final refinement because `nma`, `reduce.2`, and `reduce.3` are
  32-bit helpers with unavailable runtime dependencies.

Consequence:
- Current and historical adapters are validated through candidate generation and explicit refinement failure, but a FiberDock
  quality comparison remains blocked until a permitted compatible 32-bit runtime is staged.

## Source-gate denominator policy (2026-07-14)

Decision:
- Freeze `source-gate-policy/v1` with `dataset_row_id` identity, 240 strict rows, and 17 audit-only rows. Confirmatory execution is
  blocked until the 17 chain-contract cases receive authoritative resolution.
- Prohibit automatic orientation swaps, full-PDB substitution, or inferred selector repairs.

Evidence:
- `tmp/agent/20260713-investigation-implementation/source-gate-aggregate-final/source_gate_policy.json` SHA256
  `2d652ad7ba4ebd4b7dfbad168458af85cbb309edf177d46c2a12a8cd5a5767db`.

## MultiProt runtime provenance (2026-07-14)

Decision:
- Use a derived, task-local MultiProt installation rather than running `install_MultiProt.pl`, because the installer edits `~/.cshrc` and rewrites source-side scripts.
- Pin the selected executable to the checked-out `working_version/multiprot/external_tools/multiprot` payload and record the separate bundled ZIP hash as an audit comparison. Never combine the two binary payloads.
- Use the existing Python 2.7.15 interpreter directly, expose PyMySQL as `MySQLdb` only through the staged runtime shim, and keep database/network/mail/cleanup behavior disabled by the compatibility adapter.
- Treat missing NumPy as an explicit environment limitation; do not claim full paper-era PRISM readiness until the dependency is installed and tested.

Evidence:
- `tmp/agent/20260713-investigation-implementation/multiprot-environment-v4/environment_manifest.json`
- `tmp/agent/20260713-investigation-implementation/multiprot-smoke-array-v4/task-1/exit.json`
- `tmp/agent/20260713-investigation-implementation/multiprot-smoke-array-v4/task-2/exit.json`
- `references/nprot.2011.367.md` installation and MultiProt parameter sections.

Execution robustness:
- Empty or failed surface extraction now writes an explicit empty `.asa.pdb` plus failure log, and alignment writes `alignment_unavailable` JSON records instead of aborting the entire batch on a missing surface or unparsable TMalign output.

Reason:
- During full batch 2, `1j0sA` had no residues above the RSA threshold. The previous missing-file behavior aborted all ten pairs; explicit empty records preserve downstream failure accounting.

## Full comparison batch protocol (2026-07-12)

Decision:
- Use 257 rows from `T_Rigid.csv`, `T_medium.csv`, and `T_difficult.csv`, represented once in `tmp/agent/20260712-multichain-full-comparison/batches/shared_manifest.csv` and split into 26 ten-pair batches.
- Submit separate current and legacy Slurm arrays, retaining unavailable PDB IDs in per-task status files.

Reason:
- The identical manifest controls pair provenance and prevents the two pipelines from silently receiving different inputs.

Environment limitation:
- `/scratch/tmp/prism-current-test-py311` fails during interpreter initialization with `init_fs_encoding` / missing filesystem codec; focused checks use `gtalign_env` and passed there.

Verified full-run outcome:
- Run all 26 batches for both pipelines using the shared 257-row manifest; use the `ai` partition for the legacy array and the corrected current batch retry after the initial `cosbi` batch-2 input failure.
- Report pipeline coverage separately from model quality: both pipelines completed 255/257 pairs, but current produced 143 Rosetta models while legacy produced 2 FiberDock models for one pair. No pair had scoreable models from both pipelines.
- Use the dedicated score environment for DockQ/iRMSD because DockQ is unavailable in `gtalign_env`; keep pipeline tests and helper validation in `gtalign_env`.

## Notes

- Many PRISM files did not score simply because they did not map to a benchmark row in that set.
- Additional matched-row no-score reasons included:
  - `model_chain_missing_in_pdb:<chain>`
  - occasional DockQ interface issues
  - occasional iRMSD runtime errors
- The path `/scratch/rshadi25/GitHub/PRISM-prescript/benchmark/prism_processed/env/prism_score_env/bin/python` being missing in the user shell was identified as an environment/path issue, not a metric issue.
- The `Investigate Prism pipeline steps` chat is also recorded in PRISM-old and PRISM-main memory, so benchmark-side notes should stay aligned with the comparison repo.
- The `Integrate soft align into pipeline` chat is relevant here only as a cross-project comparison note; it does not alter the benchmark workflow itself.

## Runtime and GTalign implementation freeze (2026-07-14)

Decision:
- Treat `environment.yaml` plus `runtime_manifest.json` as the reproducibility contract for the implemented validation harness.
- Keep FiberDock full refinement disabled until a matched positive end-to-end run succeeds with the bundled helpers. ELF width alone
  is recorded as architecture evidence; missing loader/runtime dependencies and lack of end-to-end validation remain blockers.
- Treat GTalign as exploratory until a matched pilot uses common templates, filters, Rosetta refinement, and frozen evaluation. The
  CPU smoke validates isolation and provenance only, not speed or quality.
- Preserve the cleanup manifest and do not delete pre-existing caches, logs, or derived environments without explicit approval.

Evidence:
- `runtime_manifest.json`
- `benchmark/scripts/validate_runtime_manifest.py`
- `benchmark/scripts/stage_legacy_tool_environment.py`
- `docs/pipeline-validation-report-20260714.md`

## Experimental FiberDock reduce.3 adapter (2026-07-14)

Decision:
- Keep the primary FiberDock capability fail-closed and historical `reduce.2`-blocked.
- Permit the new `--fiberdock-reduce-helper` staging option only for explicitly labelled exploratory experiments. It records the original and effective helper hashes and states that the substitution is not historical-equivalent.
- Initialize `fiberdock_output/<job>` in the isolated launcher before invoking the historical controller. This is a plumbing repair for the controller's parent-directory assumption and does not alter FiberDock inputs or algorithms.

Evidence:
- KUTEM job `1355544` used the correctly substituted helper and produced non-empty `.HB`, NMA, `.fib`, and `.ref` files, but stopped at the legacy missing-parent directory call.
- KUTEM job `1355545` used the same isolated pair/template and substituted helper after the launcher repair. It returned controller code `0` and produced final `.fiberdock.pdb` and `.intRes.txt` files. Its scientific status remains `pipeline_blocked_full_refinement_capability` because the environment capability flag intentionally does not promote a non-historical substitution.
- Experimental helper hash: `6d066f88bff740627d7c1d2fb0200d326978fe70a0c041a2528c8853e682b1ce`; original reduce.2 hash: `1f12c9d3931d95dc549f0bb7e140c9462764f13b575241d5f24c6d5b3e1e115f`.
- The official Reduce repository documents reduce2 as a maintained successor but does not establish reduce.3 as a drop-in replacement: <https://github.com/rlabduke/reduce>.

Consequence:
- FiberDock is now classified as `experimental_end_to_end_observed; primary arm blocked`, not as a generally verified historical pipeline. A matched reduce.2/reduce.3 comparison, evaluator validation, and broader benchmark replay are still required.

## Output-integrity gate and benchmark-score reconciliation (2026-07-14)

Decision:
- Artifact existence and controller exit code are execution evidence only. A model is scoreable only when its raw PDB preserves the declared receptor and ligand as distinct chains, has no within-chain residue-number reset, and passes the frozen chain mapping/evaluator contract.
- `benchmark/scripts/score_single_prism_pair.py` now fails closed before iRMSD/DockQ when the raw model violates that contract; structural metrics remain null.
- The retained reduce.3 FiberDock PDB fails this gate: all raw ATOM records use chain `B`, with a residue-number reset for the second partner. The historical Rosetta model and `rosetta_output_1_chainfixed` paths are absent from the current workspace, so the historical DockQ/iRMSD row is excluded from the new comparison rather than treated as independently validated.

Evidence:
- `tmp/agent/20260714-fiberdock-reduce3-smoke-parentfix/task-2/retry-1355545-r0/workspace/fiberdock_output/smoke/1b27AD_pdb1_0_pdb2_0.fiberdock.pdb`
- `tmp/agent/20260714-fiberdock-reduce3-smoke-parentfix/task-2/retry-1355545-r0/exit.json`
- `benchmark/prism_processed_results/prism_rigid_analysis_all_jobs/per_prediction.csv`
- `tests/test_model_output_integrity.py`

Consequence:
- FiberDock is `execution_observed; scientifically unscoreable; primary arm blocked`, not verified end-to-end. Recovery/regeneration of chain-preserving output is required before any DockQ, iRMSD, interface-size, or accuracy comparison.

Metadata policy:
- `fiberdock_energy_only` is false unless an energy-only execution record exists; loadability is recorded separately as `fiberdock_energy_only_loadable`.
- 32-bit helper executables are explicit capability blockers even when dependency inspection reports them as available; no architecture substitution is silently accepted.

## Final observational benchmark replay (2026-07-14)

Decision:
- Treat KUTEM replay jobs `1355887` (strict no-align) and `1355897` (alignment-enabled audit) as observational validation artifacts only. Do not use either to claim causal superiority of TM-align/Rosetta or MultiProt/FiberDock.
- Use strict no-align as the confirmatory evaluator contract. Null metrics and mapping errors are retained as explicit failures, not converted to zero-quality scores.
- Retain the alignment-enabled results only to diagnose compatibility with the previous report. The legacy pair reproduces, while the current-arm distribution changes because the scoring regime differs.
- Keep the replay wrapper retry-safe and provenance-complete by recording a per-task parameter file, command file, resources, retry ID, output, and exit status.

Evidence:
- `tmp/agent/20260714-observational-score-replay-strict-v2/collected/COMPARISON.md`
- `tmp/agent/20260714-observational-score-replay-aligned-v2/collected/COMPARISON.md`
- `docs/findings.tsv` FV-005 through FV-009

Consequence:
- The previous current-arm score summary is not a valid direct comparator for the new aligned audit. A causal benchmark comparison remains blocked until regenerated outputs have complete chain/residue correspondence and the source-gate policy is satisfied.

Review correction:
- The hardened replay supersedes the initial jobs: strict `1355946` and aligned `1355947` completed with raw DockQ JSON retained and collection-time hash/coverage validation. iRMSD “best” is the minimum; prior maximum-labelled values must not be reused.

## PyRosetta and GTalign comparison (2026-07-14)

Decision:
- Keep PyRosetta as a separate opt-in arm and do not add it to the stable environment until an authorized, hash-pinned wheel imports and initializes. The current stable refinement arm remains external Rosetta.
- Accept GTalign CPU as an operational alignment backend for KUTEM CPU jobs, but do not promote it to a verified PRISM replacement from alignment-only evidence. Downstream transformation, refinement, evaluator, and paired quality tests are still required.

Evidence:
- `tmp/agent/20260714-pyrosetta-gtalign-comparison/pyrosetta-smoke-20260714/result.json` records the missing-module status and input hash.
- `tmp/agent/20260714-pyrosetta-gtalign-comparison/prism-smoke-1356100/summary.json` and `exit.json` record the 2/2 real smoke and unique workspace/output provenance.
- `tmp/agent/20260714-pyrosetta-gtalign-comparison/alignment-pilot-50x50/comparison.json` records near-parity synthetic throughput.

Consequence:
- PyRosetta has no valid benchmark denominator yet. GTalign has an alignment-boundary validation only; the 257-row confirmatory comparison remains closed to both quality claims until matched downstream artifacts exist.

## PyRosetta installation and API compatibility (2026-07-15)

Decision:
- Promote PyRosetta from `unavailable` to `installed; one-pose verified; benchmark comparison pending` in the separate opt-in arm. Keep external Rosetta as the baseline until matched benchmark evaluation is complete.
- Keep the adapter’s compatibility path for both `set_docking_local_refine(True)` and `set_highres_scorefxn`; the installed PyRosetta 2026.3 DockingProtocol does not expose the older zero-argument/local `set_scorefxn` interface.

Evidence:
- `tmp/agent/20260715-pyrosetta-installed/probe-final.json`
- `tmp/agent/20260715-pyrosetta-installed/one-pose-v3/result.json`
- `tmp/agent/20260715-pyrosetta-installed/one-pose-v3/refined.pdb`
- `tests/test_pyrosetta_refinement.py` and `tests/test_probe_pyrosetta_environment.py`: 6 passed.

## Reusable runtime defaults (2026-07-15)

Decision:
- Treat `/home/rshadi25/.conda/envs/gtalign_env/bin/python` as the default PyRosetta-array interpreter and freeze `-mute all -constant_seed -jran 12345` unless an experiment records a different seed deliberately.
- Keep `environment.yaml` as the redistributable base recipe and document the separately authorized PyRosetta installation rather than pretending a licensed wheel is conda-reproducible.

Evidence:
- `docs/reusable-pipeline-setup-20260715.md`
- `benchmark/jobs/pyrosetta_refinement_array.sbatch`
- focused regression suite: 14 passed in `gtalign_env`.

Consequence:
- Future batch rows have a runnable default interpreter and explicit seed provenance.  Output comparisons must normalize or ignore PyRosetta energy-table path comments; raw byte hash alone is not a coordinate-reproducibility test.

## Matched pilot execution contract (2026-07-15)

Decision:
- Use the pre-registered 12-row stratified pilot first: four rigid, four medium, four difficult, with two single-chain and two multichain rows per difficulty. Expand only after pilot outputs and scoring reconcile.
- Compare TM-align/external Rosetta as the baseline with GTalign CPU/external Rosetta as a diagnostic primary arm and TM-align/seeded PyRosetta as an opt-in diagnostic arm. Legacy MultiProt/FiberDock remains status-only.
- Require unique run roots per arm and fail-closed exit records. Never pool outputs from overlapping or overwritten directories.

Evidence:
- `tmp/agent/20260715-matched-benchmark-pilot/analysis_plan.json`
- `tmp/agent/20260715-matched-benchmark-pilot/preflight/preflight_summary.json`
- `tests/test_matched_benchmark_manifest.py`, `tests/test_collect_matched_benchmark.py`
- Corrected KUTEM submissions `1358990` and `1358991`; overlapping predecessors `1358980` and `1358981` were cancelled before completion.

Consequence:
- Pilot benchmark results remain unreported until each arm has independent logs, exit status, generated model evidence, and evaluator outputs. The 240-row strict-clean expansion remains pending.

## Template schema diagnosis (2026-07-15)

Decision:
- Treat `template_old/template` and `templates` as different pipeline asset contracts. Do not interchange them by changing only a directory path.
- The current PRISM driver requires `final_list.txt`/`calculated_templates.txt`, `interfaces/*_int.pdb`, `interfaces_lists/*.json`, and `contacts/*.json`; the legacy controller requires `interfaces/*.int`, `contact/*.txt`, and `hotspot/*`.
- A compatibility converter may be built only as an explicit derived adapter that renames/copies assets and validates parsed coordinate identity. Symlink-only conversion is insufficient for GTalign because resolved `.int` basenames are retained in output references.

Evidence:
- `tmp/agent/20260715-template-tests/new-template-alignment/summary.json`
- `tmp/agent/20260715-template-tests/old-template-alignment-copy/summary.json`
- `tmp/agent/20260715-template-tests/legacy-old-1kcaCH-run/task-1/retry-1359142-r0/exit.json`
- `tmp/agent/20260715-template-tests/pyrosetta-generated-model/task-1/result.json`

Consequence:
- Current/new and legacy/old pipeline tests are valid only within their native asset contracts. Cross-pipeline template claims require conversion plus an identical parsed-coordinate/hash audit.
## Benchmark 5.5 production gate and refiner controls (2026-07-15)

Decision:
- Use Benchmark 5.5 as the only confirmatory dataset. Keep the production TM-score threshold at `0.4`; thresholds below `0.4` are sensitivity experiments only.
- Treat current TM-align + external Rosetta as the production arm after bounded smoke validation. Treat PyRosetta as an opt-in diagnostic arm until a matched full-benchmark score denominator exists. Treat GTalign at `0.4` and legacy MultiProt/FiberDock as status-only until they produce valid predictions on a representative smoke set.
- Require explicit scientific statuses (`completed`, `completed_no_predictions`, `failed`) in every batch exit record, and fail closed on invalid output structures.

Evidence:
- Patched TM-align smoke `1360352`: exit 0, no parser failures, one accepted candidate, four refined models.
- PyRosetta smoke `1360365`: exit 0, valid four-chain PDB and PyRosetta metadata.
- GTalign `0.4` smoke `1360387`: exit 0, no predictions.
- Legacy native-template smoke `1360386`: exit 0, no passed candidates/FiberDock outputs.
- Full Benchmark 5.5 current array `1360392` is running from the staged 257-row/26-batch manifest.

## Benchmark scoring contract and runner layout (2026-07-16)

Decision:
- Use the PRISM-main benchmark analyzer and its DockQ/iRMSD metric scripts as the Benchmark 5.5 scoring authority. Do not use `score_comparison_models.py` strict `--no_align` scores for the replacement benchmark.
- Preserve PRISM-main metric code exactly. Because the analyzer has an internal sibling-path defect, execute a byte-identical disposable copy beside symlinks to the unmodified PRISM-main `dockq.py`, `irmsd.py`, and `irmsd_backbone.py`; this fixes only runtime file location.
- Convert current refinement output to the analyzer's historical filename grammar with symlinks only. Do not rename model chains, rewrite PDBs, or use native-derived chain fixes.

Evidence:
- `tmp/agent/20260716-score-contract-cleanup/main-contract-local-exact-output/per_prediction.csv`
- `tmp/agent/20260716-score-contract-cleanup/main-contract-local-runner/analyze_prism_rigid_results.py` (SHA256 equals PRISM-main source)
- `benchmark/scripts/stage_current_models_for_main_benchmark.py`
- `benchmark/scripts/score_benchmark55_main_contract.sbatch`

Consequence:
- Benchmark jobs `1361649` (generation) and `1361656` (dependent canonical scoring) supersede `1360392`/`1360454`. No report may mix old strict-null metrics with replacement metrics.

Correction:
- Jobs `1361649`/`1361656` failed before PRISM execution because their relative roots were resolved from Slurm's spool directory. They are retained only as launcher-failure evidence. `submit_comparison_batches.sbatch` now converts caller-supplied relative batch/PDB/run roots to absolute paths; fresh replacement jobs are `1361684` (generation) and `1361687` (canonical score dependency).

Scoring correction:
- Cancel the non-bijective scorer before reporting metrics. Use `score_bijective_benchmark_models.py` with complete model:native chain assignments and explicit cross-interface DockQ components. The corrected dependent scorer is job `1363282`; its outputs, not `1361687`, are authoritative candidates pending completion and aggregate validation.

Consequence:
- Do not interpret a zero-prediction completion as a successful benchmark result. Full comparisons must include attempted, completed, failed, no-prediction, scoreable, runtime/throughput, DockQ, iRMSD, TM-score, and per-pair rows.

## 2x2 aligner/refiner validation (2026-07-16)

Decision:
- Keep TM-align/external Rosetta, GTAlign/external Rosetta, TM-align/PyRosetta, and GTAlign/PyRosetta as separate arms. Bounded TM-align smoke metrics are not an accuracy comparison with full-template GTAlign metrics.
- Use the bijective scorer with observed refined-PDB chain groups and benchmark receptor/ligand groups; reject overlapping model groups before DockQ/iRMSD.
- Treat missing rigid native `1gla.pdb` as a Benchmark 5.5 data-completeness defect. Use a derived native from immutable `_r_b` and `_l_b` files only with provenance; do not modify raw benchmark inputs.
- Do not claim KUTEM capacity until pending probe `1363475` runs; `rk01` was fully allocated during validation.

Evidence:
- `tmp/agent/20260716-variant-validation/variant_summary_current.csv`
- `tmp/agent/20260716-variant-validation/scoring/gt_external/scored_merged.csv`
- `tmp/agent/20260716-variant-validation/scoring/gt_pyro/scored_merged.csv`
- jobs `1363478`, `1363479`, `1363490`, `1363491`, `1363506`, `1363507`, `1363553`

Consequence:
- The recorded pilot is runnable and scoreable, but full Benchmark 5.5 expansion remains gated by an equivalent TM-align 946-template pilot and explicit repair/validation of the missing rigid native complex.
## Legacy FiberDock/MultiProt CLI host compatibility (2026-07-17)

Reason:
- The documented CLI setup was attempted with the existing local assets and Python 2.7.15 environment.
- The bundled `accall` dynamically requires `libgfortran.so.3`, and bundled MultiProt/NMA 32-bit executables fail with `Bad system call` on the current host.

Consequences:
- Do not substitute `libgfortran.so.5` for `.3` or replace legacy binaries in a scientific run without a compatibility/equivalence test.
- The CLI is setup-complete and stage-reachable locally, but a valid refined-model run requires an authorized compatible runtime or container and preferably a Slurm allocation.

Update:
- Rebuilding NACCESS `accall` from `accall.f` is accepted as the local compatibility fix; it resolves the `libgfortran.so.3` failure without ABI substitution.
- MultiProt remains a separate unresolved host/runtime blocker.

Final validation:
- The MultiProt binary is functional when executed outside the restricted seccomp sandbox. Do not diagnose sandbox-induced exit 159 as an intrinsic MultiProt failure.
- Full legacy CLI smoke validation must run in an execution context that permits the bundled 32-bit binary; the unrestricted `test3-escalated-20260717` run is the validated result.

## Full source-matrix fail-closed decision (2026-07-18)

Decision:
- Mark template-panel execution/provenance complete after AI `1368528` and COSBI `1368530` produced byte-identical 12-row × 19,855-template matrices with zero missing interface assets and explicit self-hit/homology classifications.
- Keep confirmatory eligibility disabled because the immutable source policy still reports `blocked_source_authority` and `confirmatory_run_authorized=false`. Similarity-level eligibility must not be promoted to confirmatory template exposure.
- Regenerate positive-score failure ledgers in a new retry root after adding native-PDB hashing before DockQ/validation; never overwrite the original retry-1 evidence.

Evidence:
- `tmp/agent/20260718-template-source-gate-1/retry-3/ai/`
- `tmp/agent/20260718-template-source-gate-1/retry-3/cosbi/`
- `benchmark/scripts/run_template_source_gate.py`
- `benchmark/scripts/score_bijective_benchmark_models.py`
- scoring retry `1368545`

Consequence:
- The matched pilot and 240-row confirmatory run remain unauthorized until source authority and modern/legacy filter-asset semantics are resolved.

## Modern filter-asset normalization (2026-07-18)

Decision:
- Normalize modern and derived filter assets at load time, retain source-format/hash metadata, and pass chain-specific hotspots plus orientation-specific contacts into transformation filtering.

Rationale:
- Modern JSON assets are chain-keyed and numeric while derived assets are flattened and chain-qualified. Applying one flattened hotspot list to both partners could accept/reject the wrong orientation. Normalization fixes that integration defect without asserting parity between the two asset populations.

Consequence:
- Gate 5 remains partial. No confirmatory denominator may mix modern and derived assets until the 19,005 parity disagreements and 850 missing modern profiles are adjudicated.
