# Summary

## Summary

PRISM-prescript is the current Python 3 PRISM docking pipeline and its
benchmark/reproducibility tooling. The maintained pipeline is `prism.py`:
input download (PDB fetch) → surface extraction (NACCESS/FreeSASA) →
structural alignment (TMalign/GTalign/MultiProt) → transformation and
filtering (geometry + optional protocol hotspots) → optional candidate
ranking (top-K per pair) → selected refinement backend (external-Rosetta /
PyRosetta / FiberDock) → optional DockQ evaluation.

The legacy Python 2 MultiProt/FiberDock tree under
`working_version/Multiprot-new/prism-fiberdock-cli/` is retained only as a
reference/compatibility arm. Do not alter it when working on the current
pipeline.

## Learning workflow

- Ongoing concept learning uses `mentoring-juniors` for Socratic guidance and
  `teach` for persistent lessons, reference material, and learning records.
- The deleted learnship learning skill is not part of the user's teaching
  workflow. This preference does not change the project's separate task
  planning or routing instructions.

## Current working status

- The exact archived notes support nine confirmed pipeline combinations in the
  controlled 2026-07-21 validation. Later 18-variant BM5.5 text records
  submissions, pending GPU work, and cancelled CPU work—not a completed
  18-variant benchmark. Do not claim all 18 as verified without a new
  per-variant artifact ledger.
- The operational defaults remain NACCESS, TM-align, and external Rosetta.
- CLI parity work is captured on branch `feature/prism-cli-parity` in commits
  `004721156c2`, `56aae1fd860`, `886bea75cf3`, and `20ebf4c3cc3`. The current
  `prism.py` accepts `/PRISM` hyphenated and legacy prescript underscored
  option spellings, runtime input/backend paths, and explicit `--no-refine`;
  prescript's refinement-on default is preserved. Compatibility tests pass
  (`31 passed`); unrelated dirty benchmark files remain unstaged.
- Current pipeline interpreter:
  `/home/rshadi25/.conda/envs/gtalign_env/bin/python` (Python 3.11,
  GTalign, Biopython, PyRosetta). Use absolute GTalign paths:
  `/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu` or `gtalign_gpu`.
- Current verified DockQ scoring interpreter:
  `benchmark/prism_processed/env/prism_score_env/bin/python` with DockQ 2.1.3,
  Biopython 1.85, and NumPy 1.26.4; invoke DockQ as `python -m DockQ`.
  The historical `/scratch/tmp/prism-dockq-env/bin/python` entry now fails with
  `ModuleNotFoundError: No module named 'DockQ'` and must not be used. This
  supersedes its former active-environment status while preserving the failed
  launcher evidence.
- External Rosetta requires `module load rosetta/2022.42` plus the explicit
  `PRISM_ROSETTA_PREPACK`, `PRISM_ROSETTA_DOCK`, and `PRISM_ROSETTA_DB`
  settings documented in `docs/STABLE_PIPELINE.md`.
- The stable local smoke command is:
  `PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python bash benchmark/scripts/run_prism_pipeline_smoke.sh`.
  A completed zero-pair smoke is a health check, not a biological-positive result.
- Use isolated work directories for runs. Stable defaults are TM=0.5, minimum
  matches=15, match percentage=50, difference allowance=20, clash distance=3,
  maximum clashes=5, and scaffold threshold=5.0. Diagnostic overrides never
  become production defaults.
- Canonical current template assets are under `new_template/template/`; the
  root `templates/` directory in an isolated run is a staged/symlinked view.
  Inputs retain chain-suffixed IDs, while downloaded PDB files use four-letter
  IDs in `processed/pdbs/`.
- For a faithful historical-filter claim use `PRISM_FILTER_MODE=published_protocol`.
  `geometry_only_experimental` is appropriate for plumbing/diagnostic runs but
  is not a reproduction claim.
- Opt-in deterministic ranking is now validated end-to-end on an isolated
  paired current-pipeline smoke. Slurm job 1392725 completed `0:0` under stable
  thresholds: both arms had byte-identical inputs/source/audits and two
  transformed candidates; baseline refined two, while `--rank true --top-k 1`
  selected and refined one (`1h5bAB/o1`). This proves refinement-load
  reduction, not quality improvement. Evidence is under
  `tmp/agent/20260726-isolated-ranked-smoke/runs/1392725/`.
- A current-tree paired rerun (Slurm job 1400311) completed `0:0` with
  identical current-source/input/audit manifests: both arms generated two
  candidates, baseline refined two, and `--rank true --top-k 1` selected and
  refined one. All input, alignment, transformation, ranking (ranked arm),
  and refinement stage-status records completed successfully. In this single
  run, baseline refinement took 56.3s and ranked refinement 62.7s; therefore
  load reduction is confirmed, but wall-clock acceleration is not established
  under stochastic Rosetta. Evidence is under
  `tmp/agent/20260728-ranking-pipeline-evaluation/runs/1400311/`.
- Canonical BM5.5 scoring is now validated on the first complete July-22
  batch of `naccess_gt_external_rosetta`. Job 1392974 scored 235/235 unique
  final poses across 10 strict-clean `dataset_row_id` values in 2m54s on
  `ai08`; it produced 938 global/interface records with DockQ 2.1.3 and zero
  model/native/raw-JSON hash mismatches. Evidence is under
  `tmp/agent/20260726-bm55-canonical-scoring/scored-batch-0001-v2/`.
- Use requested receptor-ligand cross-interface DockQ for ranking labels and
  report GlobalDockQ separately. In multichain antibody-like cases,
  GlobalDockQ can be dominated by the receptor-internal interface and is not
  interchangeable with the requested cross-interface score.
- Full-variant preparation job 1393026 is accepted: 6,539 model records, zero
  blank `dataset_row_id`, 5,828 valid strict-clean poses, 359 explicit
  strict-clean model-contract failures, 282 staged audit-only poses, 70
  audit-only model-contract failures, and 156/156 represented strict-clean
  natives assembled. The eight native failures are all audit-only.
- Full strict scoring job 1393041 completed `0:0` on `ai05` with 6,539 model
  rows: 5,539 scored, 711 not scoreable, and 289 score failures. Audit job
  1393064 correctly failed `2:0`; preserve `full-v2/scored/` as failed evidence.
- The 289 failures are classified: 192 DockQ 2.1.3 empty-array crashes for
  `difficult:000016`, 96 legacy `irmsd.py` symmetric-chain `KeyError`s, and
  one genuine non-bijective partner-cardinality case (`rigid:000023`). The
  corrected full expectation is 5,827 score-bearing and 712 not-scoreable rows.
- Real `ai` diagnostic job 1393691 completed `0:0`: complete mapping reproduced
  the DockQ Cython crash, pair `BA` scored, `CA` was explicitly
  `no_native_interface`, GlobalDockQ remained unavailable, and safe grouped
  iRMSD returned 24.502. Its independent audit passed with no hash or interface
  errors under `diagnostic-cross-fallback-v1/`.
- Retry job 1393707 completed `0:0` on `ai08` with eight internal workers and
  produced 289 model rows plus 1,260 interface records. The merged table under
  `full-v2-repaired-v1/scored/` has 6,539 models, 5,827 score-bearing rows,
  712 not scoreable, 192 cross-interface-only scores, and zero hash mismatches.
  `audit-v2.json` passes DockQ/mapping/hash/interface contracts, but 96 grouped
  iRMSD values remain missing because the retry launcher mistakenly invoked
  legacy `benchmark/scripts/irmsd.py` instead of `irmsd_grouped_safe.py`.
- The safe grouped-iRMSD repair is complete without modifying the legacy
  `benchmark/scripts/irmsd.py`. Job `1394212` produced the 96-row safe-v1
  intermediate, but its original merge job `1394260` failed because 22
  auxiliary iRMSD rows still timed out and the dependent ranking job `1394261`
  became `DependencyNeverSatisfied` and was later canceled after the
  replacement completed. Preserved retries diagnosed the stale
  DockQ environment in `1400294`, 900-second symmetry timeouts in `1400320`,
  and the first cache optimization's remaining timeouts in `1403992`.
  The final pair-mask/cached implementation in
  `benchmark/scripts/irmsd_grouped_safe.py` completed all 22 residual rows in
  job `1404093` using the repository-local DockQ environment. The canonical
  merged output is `tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/`.
- Final canonical BM5.5 audit passed: 6,539 model rows, 5,827 score-bearing,
  712 explicit `not_scoreable`, 16,459 interface rows, 192 requested
  cross-interface-only scopes, zero failed auxiliary iRMSD rows, zero hash or
  interface-contract failures, and DockQ 2.1.3 throughout scored rows.
- Final ranking audit passed from the canonical output: 5,828 candidate rows,
  5,827 eligible/rankable labeled rows, one retained explicit non-scoreable
  candidate, 155 ranking groups, and zero label/hash/alignment/provenance
  errors. The deterministic baseline evaluated 58/155 native-like top-1
  groups versus 68/155 oracle-positive groups; median top-1 DockQ was
  0.02234852 versus median best DockQ 0.09339036. This is an evaluation result,
  not authorization to enable ranking by default.
- Ranking audit eligibility is repaired: explicit `score_status=not_scoreable`
  candidates remain visible for provenance but are excluded from rankable
  score-identity equality; unexpected unlabeled rows still fail. The focused
  regression suite passes 23 tests.
- Stage comparison artifacts are retained under
  `tmp/agent/20260727-pipeline-comparison/`, including
  `pipeline_stage_ledger.csv`, `pipeline_stage_ledger.json`,
  `pipeline_stage_comparison.md`, `matched_candidate_ledger.csv/json`, and
  `multiprot_gate_ledger.csv/json`. The ledgers use Bio.PDB integrity checks
  for current TMalign+external-Rosetta, current TMalign+FiberDock, current
  MultiProt experiments, and the legacy MultiProt+FiberDock workspace. No
  clean same-input/same-assets/same-evaluator comparison of the requested two
  pipelines is retained, so quality differences remain observational.
- The stage-level interpretation is recorded in
  `tmp/agent/20260727-pipeline-comparison/stage_diagnosis.md`: in the
  same-input current TMalign pilot, inputs, alignments, and transformations are
  byte-identical and the first candidate-level divergence is refinement output
  state. External Rosetta has five canonical final structures, one partial,
  and one missing; FiberDock has seven energy PDBs and one logged missing
  `zero-trial` input. The current MultiProt replay has 39 successful sides,
  only two paired successes, and zero paired orientations passing the current
  transformation thresholds because its RMSD-derived TM-score proxy is not
  calibrated to the TMalign threshold. The legacy workspace remains not
  row-pairable. The retained current MultiProt external-Rosetta and PyRosetta
  roots both stop before transformation, so their refiner backends were not
  exercised in this comparison.
- The compliant headless PyMOL renderer
  `tmp/agent/20260727-pipeline-comparison/render_1ahw_stage_comparison.py`
  completed on 2026-07-28 using `uv`/OSMesa. It loaded native, transformed,
  external-Rosetta, and FiberDock structures with 9,830, 1,672, 1,611, 6,450,
  and 4,577 atoms respectively, and produced
  `1ahw_stage_comparison.png` plus `1ahw_stage_comparison.pse`. The image is a
  qualitative integrity check, not a DockQ or quality comparison.
- The matched MultiProt/TMalign calibration is now complete in `tmp/agent/20260728-multiprot-tmalign-calibration/results/`. Slurm job `1404776` used the verified `gtalign_env` interpreter and one `ai` CPU job with eight internal TMalign workers on the exact retained panel: two queries, 100 templates, and 400 records. TMalign succeeded on 400/400 records; the retained MultiProt side remains 39 successful and 361 `alignment_unavailable`. Across the 39 matched successful records, score Pearson correlation is `-0.0009494593` and median absolute difference is `0.27306`. The current gate admits 0/200 MultiProt orientations and 1/200 TMalign orientations (`1ahwAF/o1`). This confirms the current Kabsch-RMSD proxy is not calibrated to the TMalign-oriented TM-score threshold; it does not authorize changing production thresholds or claim method quality.
- The earlier `tmp/agent/20260728-alignment-comparison/` experiment is failed evidence only: environment activation, symlink/idempotence, and import-path errors occurred, and its final JSON has zero successful alignments. Do not use it for calibration or pipeline-quality conclusions.
- The exact matched transform/refiner replay completed in Slurm job `1404934` under `tmp/agent/20260728-matched-align-refiner-replay/results-v2/`. Both arms used `1fgnHL + 1tfhA`, the same `1ahwAF/o1` template candidate, source PDBs, current transformation/clash code, copied FiberDock tools, and Rosetta 2022.42. TMalign generated the common transform; MultiProt wrote valid transform intermediates but was rejected by the clash filter. For the common candidate, MultiProt had 16 CA clashes with minimum distance 1.478 Å versus TMalign 0 clashes and minimum distance 6.015 Å; direct CA coordinate RMSDs between arm transforms were 59.8668 Å and 45.8653 Å across 428 and 202 common CA atoms. This locates the first matched deviation at transformation/clash filtering, not refinement.
- The same replay’s MultiProt-only `1ahwBC/o1` candidate generated valid FiberDock and raw Rosetta structures. FiberDock’s valid 6,765-atom PDB was paired with a blank parsed energy because the current parser looks for `fd_params.ref` while the actual file is `fiberdock_energies.ref` with `glob = 0.00`. Rosetta’s valid 9,685-atom raw structure had interaction score `0.0` and was therefore excluded by the existing `-5.0` canonical-output gate. The TMalign common candidate produced a valid canonical Rosetta structure. These are diagnostic output-contract findings, not quality scores.
- The prior matched replay job `1404921` is failed setup evidence only; it lacked the second arm’s transformation directory after module reuse. Job `1404934` fixed that isolated-workdir defect and is authoritative.
- Retained-feature reconstruction is accepted for batch 1: all 235 candidate
  rows were recovered from the batch's GTalign JSON, labeled by
  `dockq_cross_mean`, and passed the dedicated ranking audit with complete
  model/JSON/raw-output hash provenance. Full reconstruction and deterministic
  original ranking job 1393128 is permanently pending with
  `DependencyNeverSatisfied` after failed audit 1393064 and must not be reused.
- Replacement ranking job 1393757 reconstructed 5,828 candidates and attached
  5,827 canonical labels, then failed `1:0` at its audit. The sole unlabeled
  row is the correctly non-scoreable `rigid:000023` model; ranking audit must
  retain it for provenance while excluding it from rankable identity equality.
- Ranking reconstruction is accepted on batch 1: the stage-manifest adapter
  recovered 235/235 candidates from retained GTalign JSON and canonical labels
  joined 235/235 by `dataset_row_id` plus model SHA256. Ranking is now computed
  independently per durable benchmark row. Across the 10-row batch, three rows
  contain an oracle native-like candidate and the deterministic baseline chose
  a native-like top-1 for one; this is pilot evidence only. Outputs are under
  `tmp/agent/20260726-bm55-canonical-scoring/ranking-batch-0001/`.
- Re-ranking the retained 235-row labeled pilot with the current
  `biological-baseline/v2-real-coverage` implementation gives native-like
  top-1 on 3/10 groups, matching the 3/10 oracle-positive groups; median
  top-1 DockQ is 0.007527 versus 0.019942 for the per-group best candidate,
  with mean DockQ regret 0.003381. The retained older ranked CSV lacks the
  score-version column and used a different score (for example 0.66206 versus
  0.79010 for its first row), so current-code output is retained separately at
  `tmp/agent/20260728-ranking-pipeline-evaluation/`. This is a small,
  GTalign-derived pilot, not a paired current TMalign accuracy benchmark.
- Opt-in PRODIGY ranking is implemented as `--rank-method prodigy` and smoke
  tested on the retained `5zngA,4eylA` / `1a0cCD` two-orientation case. After
  fixing the PRODIGY argument order, both candidates scored successfully and
  top-1 selected `o1` (`-65.827 kcal/mol`) over `o2` (`-65.274 kcal/mol`),
  reducing the forwarded set from two to one. Evidence is under
  `tmp/agent/20260730-prodigy-ranking-paired-test/summary-corrected.json`.
  This is a candidate-selection/load observation, not a DockQ quality claim.
- BM5.5 cross-variant scoring against T_Rigid (2026-07-27): 4 completed GTalign
  GPU variants scored via `benchmark/scripts/score_bm55_models.py`. Rosetta-based
  refiners (NACCESS+ext.Rosetta: mean DockQ 0.600, 9/13 > 0.5; FreeSASA+PyRosetta:
  mean 0.558, 10/17 > 0.5) significantly outperform FiberDock (NACCESS: mean 0.167,
  2/11 > 0.5; FreeSASA: mean 0.207, 4/19 > 0.5). CPU variants (TMalign/MultiProt)
  produced zero models at 20K-template scale.

- The broader matched transformed-candidate panel is complete in
  `tmp/agent/20260728-multiprot-fiberdock-broader-replay/`: job 1405218
  replayed seven common retained candidates through MultiProt+FiberDock in
  one internally parallelized `ai` job, all seven returned successfully and
  produced valid declared-energy PDBs; job 1405237 compared them with the
  corresponding TMalign transforms. Input and source-PDB inventory hashes
  match, while alignment/surface inventories differ. Direct CA geometry
  differs for 3/7 left partners and 1/7 right partners; total 3-Angstrom
  cross-partner CA clashes are 454 (MultiProt) versus 385 (TMalign).
  MultiProt FiberDock has seven valid candidate PDBs; TMalign external
  Rosetta has five canonical, one partial, and one missing. This remains
  diagnostic evidence without native DockQ or matched refiner seeds.


- The isolated FiberDock output-contract diagnostic is complete in `benchmark/scripts/diagnose_fiberdock_output_contract.py`. On the retained matched replay, both work directories declare `fiberdock_energies` in `fd_params.txt` and contain `fiberdock_energies.ref` plus valid 6,765-atom `*.ref.pdb` outputs; neither contains the current parser's `fd_params.ref` guess. Correctly parsed global energies are 0.00 (MultiProt `1ahwBC/o1`) and 0.51 (TMalign `1ahwAF/o1`). This proves an output-path contract bug without modifying stable source. Evidence: `tmp/agent/20260728-matched-align-refiner-replay/fiberdock_output_contract.json`.
- A retained six-pair current MultiProt root exists at `tmp/agent/20260721-full-validation/naccess_multiprot_rosetta` with the same six input rows and 144 alignment records as the `tm_external` arm. It has 24 MultiProt-success records, 120 unavailable records, 14 transformed PDBs, and 7 Rosetta PDBs. It is useful observational stage evidence, but not a matched quality comparison: the aligner contracts differ and its surface/alignment/transformation artifacts are not byte-identical to the TMalign arm.

- The isolated corrected FiberDock replay is authoritative in `tmp/agent/20260728-fiberdock-corrected-replay`. Slurm job 1405134 used one `ai` allocation with two internal tasks and reran the declared FiberDock output contract for MultiProt `1ahwBC/o1` and TMalign `1ahwAF/o1`. Both returned code 0, produced valid 6,765-atom `LHA` complexes, and correctly parsed energies 0.00 and 0.51. The parameter stem was rewritten only into the fresh output root; `src/fiberdock_refinement.py` and shared tools were not modified. Headless PyMOL loaded native, transformed, and corrected FiberDock structures and wrote `corrected_fiberdock_replay.png` and `corrected_fiberdock_replay.pse`.

- The focused FiberDock output-contract fixture test `tests/test_fiberdock_output_contract.py` passes (1 test). The corrected replay and fixture test together verify the declared output-stem behavior without modifying stable source.

- The MultiProt score-contract repair is implemented in
  src/transformation.py: MultiProt records no longer use the incompatible
  TMalign TM-score threshold; they use native match-count/coverage gates,
  while TMalign/GTalign retain the stable TM-score gate of 0.5. The retained
  100-template diagnostic now reports 2/200 paired orientations eligible
  (previously 0/200), 37/39 successful sides passing, and the focused
  transformation suite passes 21 tests. This fixes premature score-gate
  attrition; downstream transform/clash behavior and native quality remain
  open.

## Required resources and path contract

Run from an isolated workspace that contains (normally as symlinks) `prism.py`,
`src/`, `external_tools/`, `templates/`, `inputs.csv`, and `processed/pdbs/`.
The canonical smoke launcher creates this layout automatically:
`benchmark/scripts/run_prism_pipeline_smoke.sh`.

| Need | Verified location or setting | Failure-safe guidance |
| --- | --- | --- |
| Pipeline runtime | `/home/rshadi25/.conda/envs/gtalign_env/bin/python` | Verify `import pandas, numpy, Bio` before launching. |
| Inputs and template assets | `inputs.csv`; `new_template/template/{pdbs,interfaces,interfaces_lists,contacts,hotspots,rsas}` | Stage/symlink them as `templates/`; do not assume a root `templates/` tree is populated. |
| Surface calculation | `external_tools/naccess/naccess`, `standard.data`, `vdw.radii`; optional `PRISM_NACCESS_EXECUTABLE` | Use `PRISM_SURFACE_BACKEND=freesasa` only with a working `PRISM_FREESASA_PYTHON`. |
| Aligners | `external_tools/TMalign`; `gtalign_env/bin/gtalign_cpu` or `gtalign_gpu`; `external_tools/multiprot.Linux` | Use absolute GTalign paths on compute nodes; never treat CPU/GPU GTalign output as interchangeable. |
| Refiners | Rosetta module + `PRISM_ROSETTA_*`; PyRosetta in `gtalign_env`; FiberDock at `external_tools/fiberdock/` | External Rosetta requires the module. FiberDock also requires `FiberDock`, `nma`, `Reduce`, `buildFiberDockParams.pl`, and `addHydrogens.pl`; historical equivalence remains unresolved. |
| Outputs and scoring | `processed/{alignment,alignment_gtalign,transformation,rosetta_refinement,pyrosetta_refinement,fiberdock_refinement}`; `benchmark/prism_processed/env/prism_score_env/bin/python` | Score assembled refined models only; preserve raw scorer JSON and mapping/hash provenance. The former `/scratch/tmp/prism-dockq-env/bin/python` is stale and unavailable for DockQ. |
| Visualization | `/home/rshadi25/.local/bin/uv`; `pymol-open-source-whl` 3.1.0.4 with NumPy 1.26.4 resolved by `uv run` | Headless OSMesa rendering works after dependency staging; preserve both PNG and `.pse` session. Do not assume PyMOL is installed in `gtalign_env`. |

For clean Slurm scoring, use the verified repository-local entry path
`benchmark/prism_processed/env/prism_score_env/bin/python`; verify `import
DockQ` before submission. The former `/scratch/tmp/prism-dockq-env/bin/python`
is retained only as failed historical evidence because it now lacks DockQ.

Important runtime controls: `PRISM_INPUTS_CSV`, `PRISM_FILTER_MODE`,
`PRISM_TMALIGN`, `PRISM_GTALIGN_PRE_SCORE`, `PRISM_MULTIPROT`,
`PRISM_FIBERDOCK_DIR`, `PRISM_SURFACE_BACKEND`, `PRISM_FREESASA_PYTHON`,
`PRISM_NACCESS_EXECUTABLE`, `PRISM_STAGE_STATUS_PATH`, and the documented
`PRISM_ROSETTA_*` settings. Ranking-only controls are `PRISM_RANK`,
`PRISM_TOP_K`, `PRISM_RANK_MIN_SCORE`, and `PRISM_CANDIDATE_AUDIT_PATH`.

## Knowledge-graph status

The broad `graphify-out/graph.json` received a code-only incremental update on
2026-07-25 (208,846 nodes, 337,510 links, 11,852 communities). Its exported
report/corpus metadata remains stale and its IDs predate Graphify #1504, so
treat it as navigation only—not authority for exact paths or execution
decisions. A chronology rebuild must retain the run notes, manifests, status
files, and logs under `tmp/agent`; exclude only copied environments, vendored
dependencies, caches, and self-generated Graphify corpora within those run
directories. Explicitly decide that filtered scan scope before rebuilding.

`docs/chronology/` is now the separate, deterministic chronology graph: 167
dated events across 22 calendar dates, with 199 nodes and 360 directed edges.
It preserves source paths in `manifest.tsv`, uses only cross-day
`before` relations, and provides a `current_as_of` anchor; same-day ordering
is intentionally not inferred. The builder also inspects bounded
`runs/<job-id>/results.tsv` files without recursing through copied run trees;
the 2026-07-26 ranked smoke is therefore recorded as `recorded:success` with
both job result tables retained as evidence. Rebuild
with the Graphify interpreter, not `gtalign_env`:

```bash
$(cat graphify-out/.graphify_python) \
  tools/build_project_chronology.py --root . --out docs/chronology
```

## Verification and debugging

Run the stable isolated health check before debugging a new pipeline failure:

```bash
PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  bash benchmark/scripts/run_prism_pipeline_smoke.sh
```

For ranking changes, run the focused regression suite with the verified
pipeline interpreter:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_candidate_selector.py tests/test_candidate_audit.py \
  tests/test_prism_cli.py tests/test_candidate_ranker.py \
  tests/test_ranking_data.py tests/test_train_reranker.py
```

The expanded focused suite covering selector, audit, ranker, tables, metrics,
transformation audit, and CLI passed 28 tests on 2026-07-26.

For any run, inspect evidence in this order: launcher and `run.log`; terminal
stage-status JSONL; ranked candidate-audit JSONL when ranking is enabled; raw
alignment output plus parsed alignment records; refined-pose integrity; then
raw evaluator JSON. Do not infer success from directory/file counts.

Legacy-only debugging: the 32-bit MultiProt/Reduce tools in
`working_version/Multiprot-new/prism-fiberdock-cli/` can exit 159 with `Bad
system call` when executed under a restricted seccomp sandbox. That is an
execution-context block, not proof that MultiProt is intrinsically broken.
Use an approved unrestricted/Slurm context and validate parsed alignment output
before diagnosing a legacy-binary failure or substituting tool versions.

## Confirmed fixes and bugs

- 2026-07-25 ranking integration: `--rank` now creates a fresh, run-scoped
  audit at `processed/candidate_audit/<timestamp>-<pid>.jsonl` unless
  `--candidate-audit-path` is explicit. The root cause was import-time capture
  of `PRISM_CANDIDATE_AUDIT_PATH` in `src/transformation.py`, before CLI
  parsing. Audit resolution is now runtime-configurable.
- 2026-07-25 ranking safety: a partial/stale audit cannot silently empty a
  nonempty transformation result; unmatched pairs are retained. `--top-k` must
  be positive. Focused ranking, audit, CLI, and table tests passed (21 tests).
- 2026-07-26 verification: the documented six-file focused ranking suite
  passed 24 tests with 1 expected skip under `gtalign_env`.
- 2026-07-26 isolated ranking smoke: job 1392708 verified top-2 ranking
  plumbing with two candidates in both arms; the evidence-driven top-1 job
  1392725 then reduced refinement from two candidates/structures to one with
  identical source, input, and candidate-audit manifests. All expected stage
  records completed, all pipeline returns were zero, and refined PDBs contained
  6,419 ATOM records on chains A/H. External-Rosetta accepted-energy rows were
  stochastic between arms, so this is resource-reduction evidence only.
- 2026-07-26 BM5.5 staging: external Rosetta leaves sequential
  `_rosetta.pdb`, `_rosetta_0001.pdb`, and `_rosetta_0001_0001.pdb`
  artifacts. Staging now selects exactly one most-refined artifact per pose,
  while retaining lower-suffix fallback for incomplete historical layouts.
- 2026-07-26 multichain staging: model partner groups are now split from query
  receptor/ligand chain counts, matching `src/rosetta_refinement.py`; the old
  template-chain-count split corrupted multichain mappings. Durable
  `dataset_row_id`, raw selectors, model hashes, and explicit chain-count
  failures are retained.
- 2026-07-26 canonical native/scoring: row-specific native PDBs are assembled
  only from hash-verified curated BM5.5 `r_b/l_b` role files. The strict scorer
  accepts explicit native paths, complete within-partner bijections, raw DockQ
  JSON, requested cross-interface metrics, grouped forward/reverse iRMSD, and
  source-gate audit-only exclusions. Fourteen focused tests pass.
- 2026-07-26 clean-batch DockQ fix: job 1392965 retained 235 explicit failures
  (`No module named DockQ`) because the shard launcher used `Path.resolve()`
  on the virtual-environment Python symlink. Preserving the entry path fixed
  the environment; replacement job 1392974 scored all 235 rows successfully.
- 2026-07-26 full-score repair diagnosis: complete-mapping DockQ 2.1.3 can
  crash on an unrelated empty multichain interface. Recovery is limited to the
  verified `Buffer has wrong number of dimensions` signature and scores only
  requested receptor-ligand interfaces. Every requested pair is retained as
  `scored` or `no_native_interface`; GlobalDockQ is explicitly unavailable.
- 2026-07-26--28 grouped-iRMSD repair: `getSymmetricChainsList()` deletes chain
  keys while later iterations still access them and combines independent
  permutation counts by addition rather than Cartesian product. The safe
  wrapper initially reproduced the expensive repeated interface scan: one
  six-chain `rigid:000024` row and twenty-one eight-chain `rigid:000112` rows
  timed out, with the latter expanding two 4-chain symmetry groups. The final
  implementation caches aligned chain pairs and reuses residue-index interface
  masks across orders, preserving legacy coordinate ordering while avoiding
  repeated scans. Jobs `1400320`, `1403992`, and `1404093` document the
  progression; all 22 residual rows scored in `1404093`. The legacy
  `benchmark/scripts/irmsd.py` remains unchanged.
- 2026-07-28 DockQ environment correction: retry `1400294` used the former
  `/scratch/tmp/prism-dockq-env/bin/python` path and failed all 22 rows because
  DockQ was absent. The repository-local
  `benchmark/prism_processed/env/prism_score_env/bin/python` was verified with
  DockQ 2.1.3 and used for the successful repair and final audit.
- 2026-07-26 ranking-label contract: `attach_native_labels.py` now reads the
  strict TSV schema, requires `dataset_row_id` plus model SHA256, excludes
  audit-only/failure rows, and labels from `dockq_cross_mean` rather than
  GlobalDockQ. Candidate tables preserve durable row identity, final-model
  hash, and most-refined artifact selection. The combined ranking/scoring/
  staging/label regression passed 44 tests with one expected skip.
- 2026-07-26 per-row ranking fix: `rank_candidate_table.py` previously ranked
  unrelated benchmark complexes in one global list, and `train_reranker.py`
  evaluated one global top candidate across held-out complexes. Ranking now
  groups first by `dataset_row_id` (legacy fallbacks: `native_complex_id`, then
  query pair), and reranker quality averages per native-complex group.
  `evaluate_ranked_candidates.py` records per-group top-1, oracle, and DockQ
  regret metrics. The focused ranking suite passes 13 tests with one expected
  optional scikit-learn skip.
- 2026-07-26 coverage-proxy fix: production candidate audits had never passed
  template sizes, so `match_coverage_*` was always absent and the baseline
  substituted arbitrary `mean_match_count / 50` evidence. This biased ranking
  toward longer partners. Future audits now record real per-chain template
  coverage; historical rows without denominators use TM-score alone. On the
  10-row batch-1 diagnostic, native-like top-1 increased from 1/10 to 3/10,
  matching all three oracle-positive rows. This is encouraging diagnostic
  evidence, not independent validation or authorization to enable ranking by
  default.
  Ranked tables identify this policy as
  `biological-baseline/v2-real-coverage` for reproducibility.
- 2026-07-25 DockQ compare stage: `--compare` flag added to `prism.py` as
  opt-in final evaluation after refinement. Filename-convention adapter
  `compare_pairs_from_outputs()` parses `_L.pdb`/`_R.pdb` paths to reconstruct
  metadata. 4 compare tests pass.
- 2026-07-25 module restoration: `src/eval/` (DockQ, iRMSD, CAPRI),
  `src/eda/` (template sequence extraction), `src/compare.py`,
  `src/sasa_utils.py` restored from `origin/main` and adapted to prescript
  APIs. `src/eval/irmsd.py` import fixed (`from utills` → `from .utills`).
- 2026-07-25 test adaptation: 13 test files from origin/main adapted to
  prescript's refactored API. 44 unit tests pass, 3 skipped (need real PDBs).
  Tests for removed functions (`_center_of_mass`, `_parse_matrix`,
  `_extract_chain_ids`, `_parse_tmalign_output`) rewritten for current APIs.
- 2026-07-25 graphify update: code-only incremental update (23 changed Python
  files → 179 new AST nodes). Graph: 208,846 nodes, 337,510 edges,
  11,852 communities. No LLM cost (code-only).
- 2026-07-26 chronology status extraction: the deterministic builder now reads
  immediate `runs/<job-id>/results.tsv` aggregate return codes. Two focused
  tests pass; it does not recurse through copied per-arm source/environment
  trees.
- 2026-07-27 BM5.5 model scoring: Created `benchmark/scripts/score_bm55_models.py`
  for scoring GTalign GPU pipeline outputs against benchmark CSVs. Handles 3
  refinement backends (external Rosetta, PyRosetta, FiberDock) with TER-based
  chain inference from combined model PDBs. Scored 4 completed variants:
  NACCESS+ext.Rosetta (13 DockQ, mean 0.600), FreeSASA+PyRosetta (17, 0.558),
  NACCESS+FiberDock (11, 0.167), FreeSASA+FiberDock (19, 0.207). Chain-mapping
  remains the limiting factor (~20% of rows successfully scored).
- 2026-07-27 ranking module creation: `src/candidate_selector.py` connects the
  candidate audit JSONL trail to `candidate_ranker.biological_baseline_score()`.
  `select_top_candidates()` groups candidates by (receptor, ligand), filters by
  status and optional min-score, ranks, and returns top-K. Partial/stale audits
  preserve unmatched pairs rather than silently dropping them.
- 2026-07-27 `prism.py` ranking wiring: `--rank` flag (env: `PRISM_RANK`),
  `--top-k` (env: `PRISM_TOP_K`, default 5), `--rank-min-score` (env:
  `PRISM_RANK_MIN_SCORE`), `--candidate-audit-path` (env:
  `PRISM_CANDIDATE_AUDIT_PATH`). When `--rank` is set without explicit audit
  path, a run-scoped path at `processed/candidate_audit/<ts>-<pid>.jsonl` is
  auto-created. Audit path is forwarded to `transformer()` for recording.
- The prior GTalign GPU pre-score double filter is resolved: use the separate
  `PRISM_GTALIGN_PRE_SCORE` control. Low hit rates against the 19,855-template
  panel are expected for weakly similar query/template pairs, not a pipeline
  failure.
- Current MultiProt uses native match extraction, one-letter→three-letter
  residue conversion, molecule-ID-independent coordinate matching, and Kabsch
  transforms; it is no longer a TM-align fallback.
- Rosetta failures caused by missing executables are configuration errors:
  batch jobs must load `rosetta/2022.42` and export `PRISM_ROSETTA_*` paths.
- Pipeline status is not inferred from directory/file counts: cancelled jobs
  can leave partial artifacts. Require a normal pipeline return, terminal stage
  records, refined-pose integrity, and valid evaluator output.
- 2026-07-28 stage-debugging correction: the broad refinement counts in
  `pipeline_stage_ledger.csv` count intermediate PDBs and are not candidate
  yields. Use `matched_candidate_ledger.json` for the matched TMalign arms.
  It verifies identical inputs, alignment JSON, and transformed PDBs, and
  localizes the first candidate-level divergence to refinement output state.
- 2026-07-28 MultiProt gate diagnosis: `multiprot_gate_diagnosis.json`
  replays the current transformation thresholds over 200 orientations from
  the retained 100-template run. Only two orientations have both successful
  sides and none pass both thresholds; the current MultiProt Kabsch-RMSD
  TM-score proxy therefore causes the pre-transformation drop. Do not lower
  the TM threshold or call this a quality result without a calibrated score
  contract.

## Pipeline layout

### Working directories
- `processed/` — runtime output: `alignment/`, `alignment_gtalign/`,
  `surface_extraction/`, `transformation/`, `pdbs/`, `rosetta_refinement/`,
  `fiberdock_refinement/`, `candidate_audit/`
- `templates/` — template assets: `pdbs/`, `interfaces/`, `interfaces_lists/`,
  `hotspots/`, `contacts/`, `rsas/`
- `external_tools/` — binaries: `TMalign`, `multiprot.Linux`, `naccess/`,
  `fiberdock/`, `SoftAlign/`
- `new_template/template/` — canonical current template assets (symlinked
  into `templates/` for runs)
- `src/` — current Python pipeline, ranking, refinement, utility, and scoring
  modules
- `tests/` — 76 test files (pytest)
- `benchmark/scripts/` — benchmark validation, scoring, and analysis scripts
- `fixed_pipeline/` — stable staging copy (do not modify)

### Required input files
- `inputs.csv` — Receptor/Ligand pair list (CSV with `Receptor,Ligand` columns)
- `templates/checked_templates.txt` — template manifest (one ID per line)
- `templates/calculated_templates.txt` — filtered template output list
- PDB files in `templates/pdbs/` (downloaded or staged)

### Environment variables (all configurable)
- `PRISM_INPUTS_CSV` — pair list path (default: `inputs.csv`)
- `PRISM_SURFACE_BACKEND` — `naccess` or `freesasa` (default: `naccess`)
- `PRISM_FREESASA_PYTHON` — Python interpreter with FreeSASA
- `PRISM_TMALIGN` — TMalign binary path (default: `external_tools/TMalign`)
- `PRISM_TMALIGN_WORKERS` — parallel alignment workers
- `PRISM_GTALIGN_PRE_SCORE` — optional GTalign prefilter-score override
- `PRISM_MULTIPROT` — MultiProt executable path (default:
  `external_tools/multiprot.Linux`)
- `PRISM_MULTIPROT_WORKERS` — MultiProt parallel workers (default: 8)
- `PRISM_TM_SCORE_THRESHOLD` — TM-score minimum (default: 0.5)
- `PRISM_MINIMUM_RESIDUE_MATCH_COUNT` — min match count (default: 15)
- `PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE` — min match % (default: 50.0)
- `PRISM_DIFF_PERCENTAGE` — difference allowance (default: 20.0)
- `PRISM_CLASHING_DISTANCE` — clash threshold in Å (default: 3.0)
- `PRISM_MAX_CLASHING_COUNT` — max clashes tolerated (default: 5)
- `PRISM_FILTER_MODE` — `geometry_only_experimental` or `published_protocol`
- `PRISM_FILTER_ASSET_ROOT` — protocol filter asset directory
- `PRISM_NACCESS_EXECUTABLE` — NACCESS executable path (default:
  `external_tools/naccess/naccess`)
- `PRISM_FIBERDOCK_DIR` — FiberDock tool directory (default:
  `external_tools/fiberdock`)
- `PRISM_REFINER` — default refinement backend (default: `external_rosetta`)
- `PRISM_CANDIDATE_AUDIT_PATH` — candidate audit JSONL path
- `PRISM_STAGE_STATUS_PATH` — stage event JSONL path (structured logging)
- `PRISM_RUN_ID` — run identifier (auto-generated from timestamp+pid)
- `PRISM_NATIVE_COMPLEX_ID` — native complex ID for audit metadata
- `PRISM_COMPARE` — enable DockQ evaluation (default: off)
- `PRISM_COMPARE_JOBS` — parallel compare workers (default: 1)
- `PRISM_RANK` — enable pre-refinement ranking (default: off)
- `PRISM_TOP_K` — top candidates per pair after ranking (default: 5)
- `PRISM_RANK_MIN_SCORE` — minimum baseline score (default: 0.0)
- `PRISM_ROSETTA_PREPACK`, `PRISM_ROSETTA_DOCK`, `PRISM_ROSETTA_DB` — Rosetta
  binary/database paths (for `--refiner external_rosetta`)

### Run commands

Examples use the verified absolute pipeline interpreter:
`/home/rshadi25/.conda/envs/gtalign_env/bin/python`.

**Basic (defaults: NACCESS + TM-align + external Rosetta):**
```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py
```

This command requires the Rosetta module and `PRISM_ROSETTA_*` paths above.
Use `--refiner pyrosetta` explicitly for the module-free PyRosetta backend.

**With optional stages (ranking + DockQ evaluation):**
```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --compare True --rank True --top-k 5
```

**Verified final2 ranking audit/evaluation:**
```bash
SCORE_ENV=benchmark/prism_processed/env/prism_score_env/bin/python
$SCORE_ENV benchmark/scripts/audit_ranking_table.py \
  tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/ranking/labeled.csv \
  tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/scored/scores_models.tsv \
  tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/ranking/audit.json
$SCORE_ENV benchmark/scripts/rank_candidate_table.py \
  tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/ranking/labeled.csv \
  tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/ranking/baseline-ranked.csv
$SCORE_ENV benchmark/scripts/evaluate_ranked_candidates.py \
  tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/ranking/baseline-ranked.csv \
  tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/ranking/baseline-evaluation.json
```
These commands were run successfully after the final score audit; the
candidate-builder and label-attachment inputs are retained beside them in the
same `ranking/` directory.

**FiberDock refinement:**
```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --refiner fiberdock
```

**FreeSASA surface:**
```bash
PRISM_SURFACE_BACKEND=freesasa PRISM_FREESASA_PYTHON=/path/to/python \
  /home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py
```

**GTalign GPU (full template panel):**
```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --aligner gtalign \
  --gtalign_path /home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_gpu
```

**GTalign CPU:**
```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --aligner gtalign \
  --gtalign_path /home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu
```

**MultiProt:**
```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --aligner multiprot
```

**External Rosetta (requires module load):**
```bash
module load rosetta/2022.42
export PRISM_ROSETTA_PREPACK="/opt/.../docking_prepack_protocol.static.linuxgccrelease"
export PRISM_ROSETTA_DOCK="/opt/.../docking_protocol.static.linuxgccrelease"
export PRISM_ROSETTA_DB="/opt/.../database/"
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --refiner external_rosetta
```

**Structured logging (stage events + candidate audit):**
```bash
PRISM_CANDIDATE_AUDIT_PATH=audit.jsonl PRISM_STAGE_STATUS_PATH=stages.jsonl \
  /home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --rank True
```

**Slurm batch (see `run.slurm`):**
```bash
sbatch run.slurm
```

The `run.slurm` defaults to `cosbi` partition, 8 CPUs, 40GB, 72h.

## Next steps

1. Keep deterministic ranking opt-in and use the verified paired launcher at
   `tmp/agent/20260726-isolated-ranked-smoke/run_paired_ranked_smoke.sbatch`
   when rechecking load reduction; set and record `PRISM_SMOKE_TOP_K` only for
   that experiment launcher.
2. Treat `full-v2-safe-irmsd-final2/` and its passed `audit-final.json` as the
   canonical BM5.5 scoring root. Preserve the failed `1394260`/`1394261` chain
   and retries as diagnostic evidence, not as current status.
3. Preserve the corrected ranking audit contract: explicit non-scoreable rows
   remain visible but are excluded from rankable score-identity equality. The
   final ranking audit and per-row evaluation are already complete under the
   canonical final2 ranking directory.
4. Do not enable learned reranking: the five-complex grouped evaluation did
   not improve native-like top-1 success over the deterministic baseline.
5. Add independent labeled complexes and a frozen random-seed/evaluator
   contract before making ranking-quality or model-quality claims.
6. Resolve the 17 audit-only source rows before opening a confirmatory BM5/5.5
   denominator or making full-cohort comparison claims.
7. Consult `docs/memory-archive/20260725-exact-pre-consolidation/` for the
   exact historical notes; it is reference-only, while these active files are
   the concise operational authority.
8. Phase 1 planning context is captured in
   `.planning/phases/01-run-identity-and-manifest/01-CONTEXT.md`; use its
   run-identity, artifact-ledger, fail-closed validation, and runtime-capture
   decisions before implementing provenance changes.
9. The 2026-07-29 grilling session sharpened the model: `contract_hash`
   identifies the declared source/tool/config/input contract, `run_id`
   identifies an execution attempt, and artifact/closure digests identify
   observed evidence. ADRs are recorded under `docs/adr/0001-*.md` and
   `docs/adr/0002-*.md`; the vocabulary is in `docs/glossary.md`.
10. On 2026-07-30/31, Phase 1 evidence-ledger deepening added the
    standard-library module `src/provenance/run_evidence.py` for contract and
    attempt records, row-aware artifact observations, symlink-safe hashing,
    and separate closeout views. The focused provenance test file passes 14
    tests, but this is not Phase 1 acceptance: the CLI detects a secret
    sentinel while returning exit code 0, so its consumer gate is not yet
    fail-closed. The full-suite invocation has no captured pytest summary and
    is not treated as a complete-suite claim.

- **MultiProt true TM-score contract and calibrated gates (2026-07-29)**:
  Diagnostic SBATCH job 1416942 (8 CPUs, ai partition, 1 QoS slot, PRISM_MULTIPROT_FORCE=1)
  ran with seccomp bypass on 56 pairs (2 queries × 14 templates × 2 orientations).
  47/56 pairs succeeded (was 2/28 without bypass). True TM-scores computed from
  MultiProt match_dict via Kabsch alignment of matched CA pairs using standard
  length-normalized formula: TM = (1/L_t) Σ 1/(1+(d_i/d_0)²), d_0 = 1.24(L-15)^(1/3)-1.8.
  True TM-scores range 0.0014–0.6041 (median 0.0157, mean 0.0538). Proxy TM-score
  (1-RMSD/10) is essentially uncorrelated with true TM (r=0.047). Transformation
  gates calibrated in src/transformation.py: MultiProt uses true_tm_score ≥ 0.3,
  match_count ≥ 10, coverage ≥ 30%. Only 2/47 pairs pass true TM ≥ 0.3:
  5zngA_1a0cCD_C (TM=0.3294, 11 matches) and 5zngA_1buhAB_A (TM=0.3385, 11 matches).
  TMalign/GTalign retain shared TM_SCORE_THRESHOLD=0.5 on their native TM-scores.
  Alignment JSON now includes both tm_score (proxy) and true_tm_score fields with
  tm_score_contract metadata ("multiprot_kabsch_rmsd_proxy" vs "standard_length_normalized").

- **Alignment adapter interface contract** (identified during architecture exploration):
  Three alignment adapters (TMalign, GTalign, MultiProt) write JSON to processed/alignment/
  with no shared schema. `tm_score` field means different things per aligner (native TM,
  native TM, proxy TM). Duplicate `extract_chain_and_res_ids()` in alignment.py (L176)
  and alignment_gtalign.py (L159). GTalign symlink hack makes processed/alignment stateful
  (points to last run). Next: define formal AlignmentResult protocol (Pydantic/JSON schema),
  consolidate chain/residue extraction into shared utility, remove GTalign symlink statefulness,
  add tm_score_contract field to all alignment outputs.

## Chats


### Graphify
- Main work: Graphify codebase analysis (462 nodes, 956 edges, 26 communities), implemented ranking/progression system with composite confidence scoring, and integrated full PRISM-prescript features into PRISM-main-archive on feature branches while keeping main clean.
- Last bold steps: **Graphify analysis + god node identification**; **Full prescript integration (8 modules + CLI)**
- Durable updates: decisions.md (Graphify findings, feature branch workflow, prescript-main integration); open_questions.md (no new questions; existing questions cover validation work); summary.md (this entry)
- Key files or outputs: `graphify-out/`, `src/ranking.py`, `src/candidate_audit.py`, `src/candidate_selector.py`, `src/template_filtering.py`, `src/alignment_multiprot.py`, `src/pyrosetta_refinement.py`, `src/fiberdock_refinement.py`, `prism.py` (full CLI), commits 33b0551 (PRISM-prescript) and 2571180 (PRISM-main-archive)
### PR - BM5.5 scoring repair and ranking audit
- Main work: Repaired canonical BM5.5 DockQ scoring, isolated grouped-iRMSD and ranking-audit blockers, and preserved stable pipeline versions.
- Last bold steps: none explicitly marked
- Durable updates: repaired-table status in `summary.md`; fallback/ranking contracts in `decisions.md`; remaining work in `open_questions.md`.
- Key files or outputs: `tmp/agent/20260726-bm55-canonical-scoring/full-v2-safe-irmsd-final2/`, `tmp/agent/20260727-pipeline-comparison/pipeline_stage_ledger.csv`, `benchmark/scripts/irmsd_grouped_safe.py`.

### implementing ranking option - 2026/07/27
- Main work: Wired opt-in deterministic ranking and optional DockQ comparison into `prism.py`, restored evaluation modules, adapted tests, and scored four BM5.5 variants.
- Last bold steps: **Wired candidate ranking between transformer() and refiner()**; **Wired --compare DockQ stage**
- Durable updates: ranking and DockQ decisions in `summary.md`, `decisions.md`, and `open_questions.md`.
- Key files: `src/candidate_selector.py`, `src/candidate_ranker.py`, `src/eval/`, `src/compare.py`, `benchmark/scripts/score_bm55_models.py`.

### Update and test ranking option - 2026-07-28
- Main work: Verified opt-in ranking behavior and consolidated aligner/refiner stage-localization evidence.
- Last bold steps: none explicitly marked
- Durable updates: ranking remains opt-in; matched-panel and observability questions remain recorded in project memory.
- Key files or outputs: `tmp/agent/20260728-ranking-pipeline-evaluation/`, `tmp/agent/20260728-multiprot-fiberdock-broader-replay/`.

### MultiProt diagnostic run & pipeline divergence analysis
- Main work: Completed MultiProt diagnostic SBATCH (job 1416942), computed true TM-scores for 47 successful alignments, recalibrated transformation gates, and confirmed pipeline divergence originates at alignment stage.
- Last bold steps: **Diagnostic SBATCH submission & completion**; **True TM-score computation & gate calibration**
- Durable updates: summary.md (MultiProt true TM-score contract, alignment adapter contract); decisions.md (MultiProt true TM-score contract decision); open_questions.md (MultiProt score-contract RESOLVED/ADVANCED, matched-panel ADVANCED, new: MultiProt fragment-length bias, Alignment adapter interface contract)
- Key files or outputs: `docs/alignment_comparison_report.md`, `docs/multiprot_tmscore_analysis_report.md`, `benchmark/scripts/diagnostic_multiprot.sbatch`, `src/compute_true_tmscore.py`, `src/alignment_multiprot.py`, `src/transformation.py`

### Mentoring and teaching tools
- Main work: Recorded the user's preferred learning workflow for ongoing PRISM concepts.
- Last bold steps: none explicitly marked
- Durable updates: Learning workflow section in `summary.md`; workflow decision in `decisions.md`.
- Key files or outputs: `/home/rshadi25/.agents/skills/mentoring-juniors/SKILL.md`, `/home/rshadi25/.agents/skills/teach/SKILL.md`

### codebase - Phase 1 evidence-ledger deepening
- Main work: Deepened Phase 1 run evidence into `src/provenance`, preserving compatibility adapters and adding the first artifact-ledger validation core; CLI fail-closed behavior remains incomplete.
- Last bold steps: none explicitly marked
- Durable updates: Phase 1 implementation boundary in `decisions.md`; remaining provenance ownership question in `open_questions.md`; implementation summary in `summary.md`.
- Key files or outputs: `src/provenance/run_evidence.py`, `tests/test_run_evidence.py`, `docs/exec-plans/20260730-deepen-phase1-evidence-ledger.md`, `tmp/architecture-review-20260730-154500.html`

### Planning phase 1 debugging
- Main work: Reviewed schema-alignment claims for the Phase 1 provenance prototype and isolated the remaining fail-closed CLI defect.
- Last bold steps: none explicitly marked
- Durable updates: corrected Phase 1 status in `summary.md`; acceptance blocker in `open_questions.md`.
- Key files or outputs: `src/provenance/run_evidence.py`, `src/provenance/__init__.py`, `tests/test_run_evidence.py`.

### Add CLI instructions
- Main work: Added reversible `/PRISM` CLI option compatibility to prescript,
  including input CSV, aliases, MultiProt/PyRosetta/FiberDock settings, and
  explicit refinement gating while preserving prescript defaults.
- Last bold steps: none explicitly marked
- Durable updates: CLI compatibility policy and branch/checkpoint state are
  recorded in `decisions.md` and this summary.
- Key files or outputs: `prism.py`, `src/pdb_download.py`,
  `src/alignment_multiprot.py`, `tests/test_prism_cli_parity.py`,
  commit `20ebf4c3cc3`.
