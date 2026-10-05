# Summary

## Benchmark presentation update — 2026-10-01

- `docs/PRISM-benchmark.pptx` now has 34 slides. The original 21-slide deck is preserved as `docs/PRISM-benchmark.pre-20261001.pptx`.
- Group pages distinguish 257 BM5.5 complexes, 161 pair-summary rows, and 565 DockQ-scored model rows. The September extension separates identity controls, the six-case filter diagnostic, the 257-case candidate-overlap ledger, and historical results.
- Seven new 1s78/1e6j PyMOL PNGs and matching `.pse` sessions are in `tmp/agent/20261001-prism-benchmark-slides/figures/`; the figure generator is in the framework run workspace. The 1s78 control includes TM-align self scores A=1.000 and D=0.951, DockQ=0.874; the within-case 1e6j comparison shows lower PyRosetta total score with much lower DockQ.
- The presentation was extended to 39 slides on 2026-10-01. It now reports the earlier 88-complex prediction counts alongside September candidate-row counts and includes a six-case DockQ/energy ranking-reversal table. Examples include TM-align baseline 1e6j H–P (+371.03 / DockQ 0.0054 versus +405.50 / 0.9776), TM-align relaxed 1ahw A–F (+1842.21 / 0.0165 versus +2306.29 / 0.5344), and GT-align relaxed 2igs A–D (−82.26 / 0.0059 versus +220620.67 / 0.1269). These are observational ranking mismatches between different metrics.
- A 40th slide now isolates low-energy, poor-DockQ poses: TM-align baseline 1e6j H–P (+371.03 / 0.0054), TM-align relaxed 2i25 N–O (−305.60 / 0.0301), GT-align relaxed 2igs A–D (−82.26 / 0.0059), plus 1s78, 1e6j, and 2fd6 examples. It includes available 2i25 PyMOL figures from the September figure set.

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

## Transformation attrition diagnosis (2026-08-29)

- Question: why do many structural-alignment files get removed after the
  transformation step? Evidence from the retained `1gte` GTalign+PyRosetta run
  (`tmp/agent/20260823-1gte-variant-runs/current/`, job 1593611, 19,062
  templates, queries 1gteA/1gteB):
  - Alignment stage writes only hits that already pass its own pre-filter
    (`src/alignment_gtalign.py`): `tm_score >= PRISM_GTALIGN_PRE_SCORE` (default
    0.0) AND `match_count >= 15`. Result: 2,997 JSON files, ALL `status=success`,
    ALL `tm_score >= 0.501`. So the alignment stage is the first real filter.
  - Transformation threshold gate (`alignment_passes_thresholds` in
    `src/transformation.py`: match_count>=15, tm>=0.5, match_pct>50/30) removes
    ZERO of the present files — replay over all 574 present sides: 0 failures.
  - The dominant removal is the CLASH filter (`pair_has_acceptable_clashes`,
    `CLASHING_DISTANCE=3`, `MAX_CLASHING_COUNT=5`): of 560 transformed pairs
    (o1 287 + o2 273), 489 (87%) are rejected, median 99 CA clashes, min inter-CA
    distance as low as 0.34 A. Only 71 pairs pass (matches run.log "Passed pairs
    71").
  - 1,877 of 2,997 alignment files (62.6%) are ORPHANS: a hit on one query side
    with no passing hit on the partner side, so they never form a pair. 1,540
    orphans have the other query hitting the same template+chain (orientation
    mismatch), 337 have no partner hit at all.
- Interpretation: the "removal after transformation" is mostly (a) the alignment
  pre-filter deciding which hits exist, and (b) the clash filter rejecting
  geometrically overlapping placements. The transformation threshold gate is not
  the bottleneck for TMalign/GTalign. This is a diagnostic observation, not a
  threshold-change authorization.

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

### Demo notebook and hotspot checking

- The available isolated feature-bundle demo is
  `/scratch/tmp/prism-prescript-pipeline-extension-clean/notebooks/prism_pipeline_demo.ipynb`.
  It uses the verified `gtalign_env` interpreter, defaults to read-only analysis
  of retained runs and artifacts, and makes expensive pipeline execution opt-in.
- In the current transformation code, `hotspot_analysis()` delegates to
  `evaluate_hotspots()` in `published_protocol` mode and requires at least one
  matched hotspot. `geometry_only_experimental` mode bypasses the hotspot gate.
  Hotspot/contact asset parity remains an unresolved protocol-validation issue.

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

The current `graphify-out/` navigation graph was rebuilt on 2026-09-17 from
the maintained `src/` tree only: 47 code files, 584 nodes, 1,113 edges, and
26 communities. The graph-health audit found zero missing or dangling
endpoints, self-loops, duplicate-edge collapse, or unresolved endpoint groups.
No LLM backend was configured, so community names remain deterministic hub or
`Community N` labels; this is an architecture-navigation artifact, not
scientific or execution authority.

The former broad graph (208,846 nodes, 337,510 links, 11,852 communities)
is historical and stale for exact path discovery. The full checkout was not
used for this rebuild because generated benchmark/history trees and copied
environments make an unfiltered graph noisy and impractical. A future
evidence/chronology graph must retain dated run notes, manifests, status files,
and logs under `tmp/agent` while excluding copied environments, vendored
dependencies, caches, and self-generated Graphify corpora; its filtered scope
must be decided explicitly before rebuilding.

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
  Three alignment adapters (TMalign, GTalign, MultiProt) write JSON with no shared
  schema. `tm_score` field means different things per aligner (native TM,
  native TM, proxy TM). Duplicate `extract_chain_and_res_ids()` in alignment.py (L176)
  and alignment_gtalign.py (L159). Next: define formal AlignmentResult protocol (Pydantic/JSON
  schema), consolidate chain/residue extraction into a shared utility, and add a
  `tm_score_contract` field to all alignment outputs.
  **Resolved 2026-08-06**: the GTalign symlink hack is removed; each aligner now
  writes to its own run-scoped directory
  (`processed/alignment_{tmalign,gtalign,multiprot}/<run_id>/`), so
  `processed/alignment` is no longer a shared/stateful path. See "Pipeline run commands"
  below.

## Pipeline run commands

Run from an isolated workspace that contains (normally as symlinks) `prism.py`,
`src/`, `external_tools/`, `templates/`, `inputs.csv`, and `processed/pdbs/`.
The canonical smoke launcher creates this layout automatically. Verified pipeline
interpreter (Python 3.11, Biopython, NumPy, PyRosetta, GTalign):

```bash
PY=/home/rshadi25/.conda/envs/gtalign_env/bin/python
```

### Smoke / health check
A completed zero-pair smoke is a health check, not a biological-positive result.

```bash
PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  bash benchmark/scripts/run_prism_pipeline_smoke.sh
```

### Core run commands (verified 2026-08-06)
Defaults: `--aligner tmalign`, `--refiner external_rosetta`, surface `naccess`,
refine on. Each aligner writes to its own run-scoped directory; no shared
`processed/alignment` and no GTalign symlink hack.

```bash
# Baseline: TMalign + PyRosetta (no external module; plumbing/diagnostic)
$PY prism.py --aligner tmalign --refiner pyrosetta --no-refine --template-limit 10

# External Rosetta (requires module + PRISM_ROSETTA_* env before launch)
module load rosetta/2022.42
export PRISM_ROSETTA_PREPACK="/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/bin/docking_prepack_protocol.static.linuxgccrelease"
export PRISM_ROSETTA_DOCK="/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/bin/docking_protocol.static.linuxgccrelease"
export PRISM_ROSETTA_DB="/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/"
$PY prism.py --aligner tmalign --refiner external_rosetta

# FiberDock (tangled external binary; runs pair-by-pair, slow)
$PY prism.py --aligner tmalign --refiner fiberdock

# MultiProt (only PyRosetta supported; 32-bit binary, needs unrestricted runtime)
$PY prism.py --aligner multiprot --refiner pyrosetta

# GTalign CPU / GPU (absolute paths; never treat CPU/GPU output as interchangeable)
$PY prism.py --aligner gtalign --gtalign_path /home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu
$PY prism.py --aligner gtalign --gtalign_path /home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_gpu

# FreeSASA surface backend (drop-in for NACCESS)
PRISM_SURFACE_BACKEND=freesasa $PY prism.py --aligner tmalign --refiner pyrosetta

# Optional stages: candidate ranking before refinement; DockQ eval after refinement
$PY prism.py --aligner tmalign --refiner pyrosetta --rank --top-k 1 --rank-method baseline
$PY prism.py --aligner tmalign --refiner pyrosetta --compare --compare-jobs 4
```

### Relaxed thresholds for smoke / diagnostic runs (never production defaults)
GTalign TM-scores are systematically lower than TMalign, so it needs lower
thresholds to admit the same candidates.

```bash
export PRISM_TM_SCORE_THRESHOLD=0.2                # TMalign; use 0.1 for GTalign
export PRISM_MINIMUM_RESIDUE_MATCH_COUNT=4
export PRISM_CLASHING_DISTANCE=2.0
export PRISM_MAX_CLASHING_COUNT=10
export PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE=20   # TMalign; use 10 for GTalign
export PRISM_GTALIGN_PRE_SCORE=0.05                # lower GTalign pre-filter for more hits
```

Stable production defaults (documented in `docs/STABLE_PIPELINE.md`) are
TM=0.5, minimum matches=15, match percentage=50, difference allowance=20,
clash distance=3, maximum clashes=5.

### Key runtime env controls
`PRISM_INPUTS_CSV`, `PRISM_FILTER_MODE`, `PRISM_TMALIGN`, `PRISM_GTALIGN_PRE_SCORE`,
`PRISM_MULTIPROT`, `PRISM_FIBERDOCK_DIR`, `PRISM_SURFACE_BACKEND`,
`PRISM_FREESASA_PYTHON`, `PRISM_NACCESS_EXECUTABLE`, `PRISM_STAGE_STATUS_PATH`,
`PRISM_ROSETTA_*` (see above). Ranking-only: `PRISM_RANK`, `PRISM_TOP_K`,
`PRISM_RANK_MIN_SCORE`, `PRISM_CANDIDATE_AUDIT_PATH`.

### Isolated output directories (per-tool, since 2026-08-06)
- `processed/alignment_tmalign/<run_id>/`, `processed/alignment_gtalign/<run_id>/`,
  `processed/alignment_multiprot/<run_id>/`
- `processed/pdbs/` (download), `processed/surface_extraction/`,
  `processed/transformation/`, `processed/candidate_audit/<run_id>.jsonl`,
  `processed/ranking/prodigy/`
- Refiners each keep their own root: `processed/rosetta_refinement/`,
  `processed/pyrosetta_refinement/`, `processed/fiberdock_refinement/`
- `processed/compare/` (DockQ native cache)

## Chats

### PRISM-prescript architecture graph
- Main work: Rebuilt the bounded source architecture graph and verified its graph-health audit.
- Last bold steps: none explicitly marked
- Durable updates: updated graph scope in `summary.md`, `decisions.md`, and `open_questions.md`.
- Key files or outputs: `graphify-out/graph.json`, `graphify-out/GRAPH_REPORT.md`, `graphify-out/graph.html`

### notebook and hotspot checking - 2026-08-21
- Main work: Located the isolated PRISM feature-bundle demo notebook and documented the current transformation hotspot gate.
- Last bold steps: none explicitly marked
- Durable updates: summary.md and decisions.md; no new open question because existing hotspot/contact parity coverage is sufficient.
- Key files or outputs: `tmp/prism-prescript-pipeline-extension-clean/notebooks/prism_pipeline_demo.ipynb`, `src/transformation.py`, `src/template_filtering.py`

### Local runnable PRISM feature bundle - 2026-08-16
- Prepared `/scratch/tmp/prism-prescript-pipeline-extension-clean` as a local-only frozen-case bundle with TMalign, MultiProt, FiberDock, NACCESS, `1kcaCH` template assets, and `1FGNH`/`1TFHA` inputs.
- Added preflight, isolated-run, Slurm matrix, interactive guidance, asset hashes, and read-only `origin/main` comparison helpers.
- Validation: 77 tests pass after a PyRosetta vector compatibility regression fix. Jobs 1506810, 1506813, and 1506819 confirmed TMalign, GTalign CPU/GPU execution, MultiProt execution under `Seccomp: 0`, baseline ranking, explicit PyRosetta refinement, and FiberDock refinement. PRODIGY is absent. The legacy-named `external_rosetta` path exposed PyRosetta API drift; corrected rerun job 1506821 is pending.


### Pipeline run commands
- Main work: Captured durable per-tool run commands for all PRISM pipeline variants and removed the shared/GTalign-symlink `processed/alignment` path in favor of per-tool, per-run isolated output directories.
- Last bold steps: **Test TMalign/GTalign/MultiProt alignment isolation**; **Verify refinement tools functional**
- Durable updates: summary.md (new "Pipeline run commands" section + Chats entry; corrected the GTalign-symlink note to resolved); open_questions.md (marked alignment symlink statefulness resolved).
- Key files or outputs: `prism.py`, `src/alignment.py`, `src/alignment_multiprot.py`, `src/alignment_gtalign.py`, `src/pyrosetta_refinement.py`, `tests/test_prism_cli.py`.

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

### MultiProt legacy-compatible mode - 2026-08-11
- Main work: Added an opt-in compatibility path in the current MultiProt
  adapter that preserves legacy interface/query order, optional `params.txt`,
  `Reference Molecule`, `Trans`, and three solver solutions; downstream
  transformation now evaluates retained solution variants.
- Validation: Slurm job 1495420 on `rk01` (`Seccomp: 0`) produced 8/8 current
  adapter successes with 3 solutions each. Post-change default mode job
  1495424 retained its previous 1/8 Kabsch result. Downstream variant smoke
  job 1495430 produced six transformed PDBs from the retained solutions.
- Key files: `src/alignment_multiprot.py`, `src/transformation.py`, `prism.py`,
  `tests/test_multiprot_legacy_compat.py`, `CHANGELOG.md`.

### Matched MultiProt compatibility panel - 2026-08-12
- Main work: Ran the current `legacy_compatible` adapter and the legacy
  Python-2 adapter on the same normalized 770-template panel (10 batches of
  77) for pairs `1cew/2ghuD` and `2uwjG/2uwjE`. Current job 1496043 and the
  corrected legacy-only job 1496073 used identical MultiProt binary hashes,
  the exact legacy `params.txt`, and isolated batch workspaces.
- Validation: 6,160 records and keys per arm; 6,157 successful sides and
  three no-solution/unavailable sides per arm; 6,026 successful records had
  exact retained solution payloads. Focused compatibility/parsing tests passed
  (`11 passed`).
- Finding: native current and legacy interface files are not semantically
  identical across the panel (16 residue-key-set differences and 27 CA
  coordinate-different sides among 1,540). A six-template serial probe matched
  21/48 with native current assets but 48/48 when legacy interface contents
  were staged under current filenames. This supports an asset-content cause on
  the probe, not a general downstream equivalence claim.
- Evidence: `tmp/agent/20260811-multiprot-compat-panel/aggregate/`, with
  `comparison.md`, `summary.json`, `record_comparison.tsv`, and
  `asset_parity.json`; Terra verdict: PASS WITH CAVEATS.

### Full legacy-interface asset replay - 2026-08-12
- Main work: Ran the current `legacy_compatible` MultiProt adapter across the
  same 770-template, 10-batch panel using the legacy interface contents staged
  under the current `*_int.pdb` filenames. The run used isolated per-batch
  roots and Slurm job 1496281 on `ai01`; the redundant `kutem` copy 1496273
  was canceled while pending.
- Validation: 6,160 records and keys overlapped the corrected legacy run;
  6,157 successful solution payloads matched exactly across solution count,
  solution number, match count, reference molecule, `Trans`, and
  `match_dict`; max `Trans` difference was 0.0. The same three no-solution
  keys were retained, with current status `alignment_unavailable` versus
  legacy status `no_solution`. All current subprocess return codes were zero.
- Asset audit: 1,540 interface sides were audited file-wise. 1,305 had exact
  atom identities, 3 had coordinate differences (maximum 73.887 A), and all
  1,540 differed in B-factors; semantic exactness was therefore zero under the
  audit's identity-plus-numeric-field definition. Both MultiProt binaries and
  the exact legacy params file had matching SHA256 hashes.
- Evidence: `tmp/agent/20260812-multiprot-legacy-assets-full/aggregate/`,
  `asset_audit.json`, `manifest.json`; focused tests passed (`11 passed`).
  Terra verdict: PASS WITH CAVEATS. This establishes adapter-level
  equivalence conditional on legacy interface contents, not downstream
  transformation, FiberDock, DockQ, ranking, or scientific equivalence.

### Multiprot-fix & pipeline check - 2026-08-13
- Main work: Checked current executable, environment, pipeline, and benchmark
  readiness without launching a benchmark or modifying source.
- Validation: `gtalign_env` imports Python 3.11, Biopython, NumPy, pandas,
  PyRosetta, and FreeSASA; the repository-local scoring environment imports
  DockQ; TMalign, GTalign CPU/GPU, MultiProt, NACCESS, and FiberDock binaries
  are present; focused CLI/transformation/MultiProt/FiberDock checks passed
  (`12 passed`). Rosetta 2022.42 is available as a module and must be loaded
  explicitly in Slurm jobs.
- Readiness boundary: the stable benchmark baseline is current TMalign +
  NACCESS + external Rosetta. GTalign CPU/GPU equivalence, complete current
  MultiProt end-to-end validation, and causal current-versus-legacy FiberDock
  comparison remain separate validation arms.
- The working tree contains extensive generated and untracked benchmark state;
  full benchmarks require an explicit frozen source/assets/environment
  manifest. The local Slurm controller was unreachable from the active shell,
  so scheduler state was not inferred.
- Last bold steps: none explicitly marked
- Durable updates: summary.md, decisions.md, and open_questions.md
- Key files: `docs/STABLE_PIPELINE.md`, `prism.py`,
  `src/fiberdock_refinement.py`, `tests/test_fiberdock_output_contract.py`

### Transformation gate ablation notebook - 2026-08-23
- Added and executed `notebooks/transformation_gate_ablation.ipynb` as a
  diagnostic companion for the frozen
  `/scratch/tmp/prism-prescript-pipeline-extension-clean/runs/20260816T121957Z-tmalign`
  case with canonical `new_template/template` assets.
- The notebook preserves both candidate orientations in a nine-gate ledger,
  reports independent outcomes, cumulative diagnostic attrition, current
  production-order replay, leave-one-gate-out rescues, and a descriptive
  one-variable TM-score sensitivity table. It does not alter production code
  or defaults.
- In the retained case, both orientations pass data availability, transform
  fields, minimum matches, interface coverage, hotspot mapping, complementary
  contacts, and transformation materialization. Both fail the current
  TMalign TM-score threshold of 0.5; `o2` additionally has 29 C-alpha clashes
  under the current `<3 A`, reject-at-5 contract, while `o1` has zero.
- The notebook now includes a visible `gte_against_1h7x` selector for
  receptor `1gteA`, ligand `1gteB`, and template `1h7xCD` with chains C/D.
  Selecting that case activates an opt-in preparation cell that writes the
  pair manifest, stages canonical interface/filter assets, invokes the
  downloader and TMalign pipeline in an isolated run root, and points the
  ledger at the generated alignment directory. The earlier GTalign attempt is
  retained as failed evidence; the notebook does not silently treat it as a
  completed alignment run.

### Full 1gteA/1gteB three-arm execution - 2026-08-23
- Authoritative Slurm job `1593611` is running in
  `tmp/agent/20260823-1gte-variant-runs/` with full default template panels.
- The current GTalign/PyRosetta arm has completed input, alignment, and
  transformation filtering; PyRosetta refinement remains active. The legacy
  MultiProt/FiberDock and TMalign/Rosetta arms are sequentially pending in the
  same controller job.
- Jobs `1593591` and `1593599` are setup-diagnostic evidence only: the first
  exposed incorrect symlinked workspace paths; the second was canceled after
  exposing missing canonical current template assets. Neither is a biological
  result.

### Parameterized orientation comparison controls - 2026-09-09
- Added current-CLI template selection through `--templates` or
  `--template-list`, with six-character ID validation and post-selection
  `--template-limit` behavior. Explicit panels no longer require opening the
  default template manifest.
- Added `src/transformation_config.py` with typed, environment-compatible
  transformation thresholds. CLI overrides now reach alignment prefiltering
  where applicable and all transformation, contact, hotspot, coverage, TM,
  and C-alpha clash gates.
- `--inputs-csv` is now propagated to both input download and transformation;
  the transformation stage clears per-call accumulators so sequential notebook
  arms do not mix candidates.
- Extended `notebooks/orientation_orphan_diagnostic.ipynb` with an opt-in
  current-pipeline runner for the no-option native default and explicit `o1`
  / `o2` arms, isolated run roots, stage-status/audit paths, command previews,
  and threshold-sweep planning. The notebook was executed successfully with
  31 cells and execution disabled for pipeline arms.
- Added the focused `notebooks/pipeline_orientation_comparison.ipynb`, whose
  separate guarded cells invoke the current CLI for native-default, `o1`, and
  `o2` in isolated roots. Its 12-cell disabled path executed successfully.
- MultiProt match/coverage filtering now uses the chain-specific template
  interface size as the coverage denominator; the maintained function
  inventory was regenerated after the API changes.
- Focused validation passed: 64 tests; full-suite execution remains a
  separate check because the pre-existing broad suite previously stalled.

### Stepwise orientation reliability notebook - 2026-09-10
- Added `src/stepwise_analysis.py` and extended
  `notebooks/pipeline_orientation_comparison.ipynb` to implement the full
  recommended diagnostic order: frozen baseline, input/asset parity,
  alignment/pairing provenance, independent and leave-one-gate-out ledgers,
  clash distance/event replay, refinement evidence validation, same-set
  ranking metrics, and a separate US-align preflight/contract check.
- The notebook records explicit `unknown`, `deferred`, `mixed`, and
  `not_available` states rather than treating missing evidence as success.
  Pipeline and clash execution remain opt-in.
- Added raw-output hashes and subprocess return-code fields to TMalign and
  MultiProt alignment JSONs. Added the `--scaffold-threshold` override while
  retaining the stable 5.0 Å default.
- Validation: the 28-cell notebook executed successfully in a temporary
  `gtalign_env` kernel with `RUN_PIPELINES=False`; stepwise helper tests passed
  along with the focused pipeline suite. No biological comparison was run.

### Final audit hardening - 2026-09-10
- Added interface-list JSON provenance because the transformation stage uses
  `templates/interfaces_lists/<template>.json` for chain-specific coverage.
  Consumed assets are classified by resolved path and reference-byte hash;
  copied/renamed assets are detected and missing consumed assets remain
  unresolved rather than being inferred as current.
- Alignment inventory now reports contract status for each partner. Gate
  replay forces alignment-dependent gates to `unknown` when status, return
  code, required fields, or raw-output hashes are incomplete; a parseable JSON
  record cannot create a cumulative pass by itself.
- Ranking recovery is fail-closed unless the pre-ranking panel matches
  exactly and the native labels have source, source hash, and deterministic
  mapping-hash provenance. PRODIGY affinity remains a ranking input only.
- Final software validation: 74 focused tests passed and 1 skipped; Python
  compilation, notebook schema/AST checks, and the 28-cell notebook disabled
  execution passed (18 code cells, zero notebook errors). No biological arm
  was launched in this implementation cycle.

### Cross-repository workflow hardening - 2026-09-12

- Both dirty repositories were frozen under
  `tmp/agent/20260912-workflow-freeze-retry1/` with status, commit/branch,
  tracked diff, index diff, changed paths, benchmark roots, binary hashes, and
  cluster state. The prescript untracked-path manifest is retained as a gzip
  artifact because nested generated environments expanded it to 256 MB.
- The maintained boundary is recorded in
  `docs/adr/0003-prism-repository-boundary.md`: PRISM-prescript owns the
  maintained CLI, provenance, benchmark, and scientific evidence; PRISM is a
  reviewed experimental feature source; the legacy MultiProt/FiberDock tree
  remains reference-only.
- PRISM documentation consistency and transformation invariants were repaired.
  The full PRISM hermetic suite passed 69 tests with compilation and diff
  checks.
- Prescript provenance validation now emits controlled fail-closed JSON for
  duplicate ledgers, accepts `--output` as an alias for `--out`, and has
  regression coverage for mutated artifacts, duplicate primary keys, and
  secret sentinels. The focused prescript suite passed 63 tests.
- The broad prescript suite reached 109 passed and 2 skipped before being
  interrupted after 63 seconds in
  `benchmark/scripts/build_matched_benchmark_manifest.py`; it is not a full
  suite pass. The exact stall remains an open validation issue.

### Cross-repository smoke and provenance completion - 2026-09-12

- PRISM now has restored guidance documentation, an explicit backend contract,
  orientation-safe transformation regression coverage, and observable DockQ
  command metadata. Its full hermetic suite passes 70 tests.
- PRISM-prescript now accepts a declared expected artifact inventory in TSV or
  JSON, returns nonzero for warning/fail consumer gates, rejects malformed or
  duplicate ledgers in controlled JSON, and records terminal skipped
  refinement events for no-candidate and no-refinement runs.
- Run identity manifests record attempt ID, output root, command, bounded
  dirty-tree/source hashes, Slurm array metadata, safe environment identity,
  template file hashes when a template directory is explicitly supplied, and
  executable identity without serializing credentials.
- Focused prescript validation passed 67 tests; compilation and shell syntax
  checks passed. The broad suite remains incomplete after the known stall at
  benchmark/scripts/build_matched_benchmark_manifest.py:58.
- Stable TMalign smoke job 1658925 exited 0 with all stages terminal and zero
  candidates. GTalign CPU job 1658928 and GPU job 1658929 each exited 0 with
  four raw-hashed alignment records and zero candidates; their score/match
  differences remain unresolved.
- Optional backend job 1658926 succeeded for PyRosetta and FiberDock in an
  isolated copied FiberDock workspace. DockQ was attempted through
  DOCKQ_PYTHON but failed because the installed compiled extension is
  incompatible with NumPy 2.4.6.

### Cross-repository validation update - 2026-09-12

- The first full-suite attempt was conservatively interrupted while the
  matched-benchmark-manifest test was still running. That test was isolated
  and passed all four cases in 23.03 seconds.
- After repairing the dirty benchmark evaluator's explicit-chain, repository-
  root, and JSON-output contracts, the complete prescript suite passed
  `344 passed, 6 skipped` in 122.94 seconds. One warning records that the
  installed SciPy build expects NumPy below 2.3; this is separate from the
  DockQ compiled-extension failure.
- `tests/test_model_output_integrity.py` now passes 6/6. Invalid raw chain
  contracts fail closed before either scorer, explicit mappings are honored,
  and multi-interface DockQ JSON does not promote the first interface's
  metrics to model-level fields.
- The explicit DockQ replay is preserved under
  `tmp/agent/20260912-dockq-replay/` with job `1658930`, status `failed`, and
  return code `1`; the NumPy ABI issue remains an evaluator blocker.
- A bounded run-identity smoke under
  `tmp/agent/20260912-run-identity-smoke/` passed with a one-file template
  inventory. Recursive hashing of the repository's symlinked template tree
  was stopped before output and remains unsafe without an explicit bounded
  inventory.
- Runtime-manifest validation passed. The global prescript diff check remains
  blocked by six modified benchmark CSV first-line whitespace records and one
  trailing-space line in `src/rosetta_refinement.py`; generated evidence was
  not rewritten.
- Controlled orientation notebook job `1658931` ran native-default, `o1`, and
  `o2` with the same input CSV, one verified template, TMalign, stable
  thresholds, and `--no-refine`. All three input/alignment/transformation
  stages completed; native default had two alignment-threshold rejections and
  each fixed orientation had one. Refinement was terminally skipped. The job
  was canceled during notebook provenance post-processing because the asset
  hash index was traversing the large template/reference tree; partial arm
  evidence is preserved under `tmp/agent/20260912-orientation-study/`.

### Repository-local DockQ runtime replay - 2026-09-12

- The repository contains a usable scoring environment at
  `benchmark/prism_processed/env/prism_score_env/`. Its Python module entry
  point imports DockQ `2.1.3` with NumPy `1.26.4` and completed on Slurm job
  `1658942` using the `kutem` account/QOS on `rk01`.
- The standalone `bin/DockQ` file has a stale shebang pointing to the sibling
  `PRISM` tree. Use the environment's explicit `bin/python -m DockQ` invocation
  or the canonical scoring adapter instead; do not repair or mutate the copied
  environment in place.
- The isolated replay root is
  `tmp/agent/20260912-dockq-repo-env/`. Both raw DockQ and
  `benchmark/scripts/score_single_prism_pair.py` returned zero with mapping
  `OA:GF` and score `0.2116967149685021`. Raw metrics were fnat `0.2105263`,
  iRMSD `4.5240`, LRMSD `12.2353`, F1 `0.2807`, and zero clashes; the adapter's
  auxiliary iRMSD was `4.476`.
- The model/native hashes, raw JSON hash, commands, environment observation,
  Slurm logs, and terminal status are retained in that run root. This is
  evaluator/wiring evidence for one prior model/native pair, not evidence that
  the model is biologically successful or that ranking improves quality.
- `src/eval/dockq.py` now honors an explicit `DOCKQ_PYTHON` path by invoking
  that interpreter as a module and preserving the entry path. Focused tests
  passed (`11 passed` across DockQ runtime, compare, and model-output
  integrity). Pipeline integration job `1658943` ran the real `gtalign_env`
  process with the repository-local override and reproduced the same score.
- `benchmark/scripts/score_single_prism_pair.py` now exposes its existing
  `--dockq-json-dir` capability at the CLI. Job `1658944` retained one raw
  `model_1glcFG.dockq.json` file in a fresh run root with return code 0 and
  the same DockQ score; this closes the documentation/CLI mismatch without
  changing stable scoring defaults.

### USalign manual validity pilot - 2026-09-21

- The production USalign array was canceled before debugging; no replacement
  job was submitted. A bounded manual pilot used three single-chain interface
  pairs from `templates_test`.
- TMalign wall times were 0.017-0.024 s, USalign default 0.043-0.051 s,
  USalign `-fast` 0.042-0.046 s, and MultiProt 0.037-0.048 s. All calls
  returned zero; MultiProt produced `2_sol.res` and Largest Solution values
  37, 16, and 15. The local environment therefore does not reproduce a large
  USalign speed advantage over MultiProt.
- The parser regression now recognizes USalign `Structure_1/Structure_2`
  labels and accepts explicit `aligner_name`; the focused alignment module
  passes 8 tests. The compatibility `tm_score` remains the maximum of both
  normalized scores, so strict reference-normalized scoring remains open.

### PRISM aligner comparison recovery - 2026-09-28

- The selected comparison workstream is `/scratch/rshadi25/GitHub/PRISM-prescript`.
  The compact historical 19,855-template package is retained under
  `benchmark/prism_processed_results/prism_aligner_comparison_20260928/compact_historical_19855/`.
- Exact panel evidence is now copied compactly under
  `benchmark/prism_processed_results/prism_aligner_comparison_20260928/exact_panel_evidence/`.
  It freezes 19,948 checked, 19,062 calculated, 19,058 materialized, and
  historical 19,855 panel identities with hashes and exclusions.
- The shared parser now preserves `tm_score_query`, `tm_score_ref`, and an
  explicit reference-normalized Structure_2 contract. Focused parser,
  contract, and worker-sweep tests pass (`18 passed` in the latest combined
  run); identical real TMalign/USalign inputs also passed the mapping,
  transform, RMSD, and dual-score gate.
- Exact current-panel alignment-only artifacts are reusable but not quality
  evidence: jobs 1659356/1659358 (GTalign), 1659448/1659449 (TMalign), and
  1659360/1659361 (USalign). They have no common transformation/filtering,
  ranking, refinement, or DockQ outputs.
- Corrected USalign worker/configuration sweep job 1708928 completed and passed
  all 10 configurations (default/`-fast` × 1/2/4/8/16 workers), each with
  1,892/1,892 records and zero execution failures. Default/16 reached
  105.337 records/s; `-fast`/16 reached 110.677 records/s but changed
  1,385 mappings/scores and 1,306 transforms, so default/16 is selected.
  Compact evidence is under `benchmark/prism_processed_results/prism_aligner_comparison_20260928/usalign_pilot_946/`.
- The earlier wrapper-only attempt 1708915 is preserved with its explicit
  serial-dispatch invalidation reason; its raw records/scratch were cleaned
  after the failure status was recorded.
- USalign historical-panel production array 1708992 is running on VALAR
  `kutem` as a documented fallback while KUACC remains at the per-user
  association limit. Compaction array 1709007 is dependency-queued after the
  provider, and corrected transformed-DockQ array 1709046 is dependency-
  queued after compaction. Both use run-scoped paths only; active KUACC
  refinement jobs and all canonical/source assets remain untouched.
- The transformed scoring adapter is under
  `benchmark/scripts/score_transformed_usalign_batch.py` with job wrapper
  `benchmark/jobs/score_transformed_usalign_batches.sbatch`. It uses the
  existing bijective evaluator, preserves GlobalDockQ separately from
  requested cross-interface components, checkpoints each candidate, and
  deletes combined/scoring scratch only after the checkpoint is written.
- Before enabling cleanup, the compactor was tightened to retain per-case
  query/Structure_1 and reference/Structure_2 score summaries plus explicit
  dual-score provenance on generated candidate rows. The focused suite remains
  `34 passed`; raw alignment JSON is not eligible for deletion without these
  compact summaries.

### PRISM aligner comparison continuation - 2026-09-28

- The transformed-DockQ cleanup gate now deletes transformed halves only for
  `scored`, `scored_cross_only`, or `valid_unscored` rows; `score_failed` and
  unresolved non-scoreable candidates remain auditable. The production wrapper
  now passes `--retain-transformed-inputs`, because common refinement is still
  a downstream consumer; only combined/scoring scratch is deleted at DockQ.
  The focused retention, compaction, and corrected-refinement tests pass
  (`46 passed` in the latest selected run).
- USalign compaction now joins generated candidate rows to each batch's
  `inputs.csv`, retaining `pair_id`, `benchmark_set`, source row, and complex
  identity in the compact ledger. This preserves BM5.5 case grain for the
  final EDA without filename-derived identity inference.
- Existing GTalign full-BM55 transformed scoring was reconciled as reusable
  historical-panel evidence: 15,440 candidate rows, 14,300 scored and 1,140
  score failures across 216 cases. GlobalDockQ is separate from the
  diagnostic best/interface score. A compact requested-interface table was
  rebuilt from retained raw JSON plus the frozen dataset manifest (18,052
  interface rows; zero JSON parse failures). No exact-panel or refined
  GTalign claim is promoted from this lane.
- A read-only `aggregate_corrected_refinement.py` was added and staged into
  the active KUACC run. It reads raw DockQ JSON, preserves requested
  receptor-ligand components and GlobalDockQ, computes paired FiberDock versus
  external-Rosetta deltas only when both scores exist, and writes a cleanup
  eligibility gate without deleting files.
- The first compact case-wise overlap output is retained under
  `benchmark/prism_processed_results/prism_aligner_comparison_20260928/candidate_overlap/historical_tmalign_multiprot/`.
  With candidate identity defined as template, query pair, orientation, and
  chain pair, the historical lane has 82 shared TMalign/MultiProt candidates,
  2,986 TMalign-unique candidates, 66,745 MultiProt-unique candidates, and
  mean case-wise Jaccard `0.0014676` across 257 cases. This is an execution
  observation, not a quality or causal conclusion.
- The overlap analysis now includes compact GTalign transformed candidates.
  GTalign/TMalign has 210 shared candidates (mean case-wise Jaccard
  `0.0532620`) and GTalign/MultiProt has 19 shared (mean Jaccard
  `0.0003044`); candidate presence is 216, 195, and 257 cases for GTalign,
  TMalign, and MultiProt respectively. These identities are only comparable
  where query/chain/orientation contracts agree.
- Added `benchmark/scripts/aggregate_matched_comparison.py` with tests. It
  merges compact candidate/score tables, emits one method-by-case row plus
  candidate and top-k diagnostic tables, preserves explicit failure statuses,
  and never coerces absent DockQ to zero.
- GTalign common-refinement preparation validated 14,470 materialized models
  and 970 rejected rows with explicit reasons. Array 1709167 (`0-144%8`,
  100 candidates per task) is now running on VALAR `kutem`, using the staged
  validated FiberDock/external-Rosetta/corrected-DockQ worker. Early
  checkpoints show the worker contract is executing; refinement/scoring
  outcomes remain pending and no GTalign raw/source trees have been deleted.
- Read-only corrected-refinement aggregation job 1709178 (test-only 1709177)
  is dependency-linked after 1709167. It will require all 14,470 checkpoints,
  write compact GlobalDockQ/cross-interface and paired FiberDock/Rosetta
  records, and leave deletion to the separate hash-checked cleanup gate.
- GTalign refinement currently has 24 completed and 4 running checkpoints,
  with no stage-failure records; among those records external Rosetta has 10
  completed and 9 explicit `no_model` outcomes. This is execution evidence
  only. The KUACC
  TMalign/MultiProt common-refinement array still has one running task and
  remains untouched.
- Added and tested `aggregate_usalign_batches.py` plus its Slurm wrapper. The
  read-only aggregation job 1709192 (test-only 1709191) is queued after the
  corrected transformed-DockQ array 1709046 and will accept results only when
  all 26 per-batch status/TSV packages validate.
- Added and tested `prepare_usalign_refinement_manifest.py`; it references
  retained transformed halves in place, selects only explicit refinable score
  states, and records rejected rows/reasons. No USalign raw or transformed
  artifacts have been deleted.
- Added and tested `aggregate_final_matched_comparison.py`; it will join exact
  candidate identities across transformed/refined tables and emit compact
  case-level, candidate-level, ranking, and paired-delta outputs with explicit
  normal-approximation confidence intervals.
- A fail-closed handoff job 1709233 (dry-run 1709232) is now dependency-linked
  after USalign batch aggregation 1709192. It prepares the common-refinement
  manifest only after `validated_compacted` and does not submit nested Slurm
  work or delete transformed structures.
- Latest verified scheduler/artifact poll: USalign batches 1--4 remain active
  with roughly 398k/392k/406k/426k raw alignment records and no compact batch
  status files. GTalign has 24 completed and 4 running checkpoints; FiberDock
  has 24 scored records and external Rosetta has 12 completed plus 12 explicit
  `no_model` outcomes, with no stage failures observed.
- Added and tested `aggregate_timing_resources.py`; it emits compact stage
  timing/resource rows and preserves empty resource fields when a source does
  not record CPU/GPU allocation data.
- Added and tested `cleanup_usalign_run.py`; its dry-run/apply gate verifies
  retained aggregate hashes, the refinement handoff, common-refinement
  aggregate eligibility, and final-package validation before planning removal
  of only `current/batch_*`, refinement `results`, and `adapter_inputs`.
  No cleanup has been applied because the live consumers are incomplete.
- Re-ran the comparison-focused regression set after the cleanup and parser
  provenance updates: `74 passed`. The DockQ runtime fixture now explicitly
  verifies both mapping placement and the intentional `--n_cpu 1` argument.
- Latest live poll at 2026-09-28T05:54:41+03:00: USalign provider tasks 1--6
  are active or have just started, later array tasks remain pending, and no
  compact status files exist. GTalign common refinement has 92 checkpoint
  files: 87 top-level completed, 4 running, and 1 failed. The failed
  `medium_1wq1_045 / 1de4AC / o1` record is an explicit
  input-normalization mismatch (native ligand `G*` has two chains; predicted
  model chain `D` has one), retained without unsupported reconstruction.
  FiberDock has 90 completed stages and 86 DockQ scores, external Rosetta has
  43 completions and 43 explicit `no_model` outcomes.
  KUACC task `3132295_874` remains active; its documented GPU helper path was
  unavailable, so direct scheduler inspection was used. No cleanup or job
  intervention was performed.
- Follow-up scheduler poll at 2026-09-28T05:58:33+03:00: USalign provider
  tasks 1, 3, 5, 6, and 7 were active, later tasks pending, and no compact
  marker existed; all dependent stages remained held. GTalign stayed at 92
  checkpoint files (87 completed, 4 running, 1 explicit failure). No cleanup
  or job intervention was performed.
- At 2026-09-28T06:00:11+03:00, USalign task `1708992_1` exited `0:0` with
  `available=10 unavailable=0` and empty stderr. Tasks 3, 5, 6, and 7 remain
  active; later tasks and all dependency stages remain pending, so the single
  completed task is not promoted to a validated provider-stage result.
- At 2026-09-28T06:12:29+03:00, GTalign refinement advanced to 93 checkpoint
  files: 88 completed, 4 running, and 1 explicit failure. USalign tasks 3, 5,
  6, and 7 remain active with no compact marker; no dependency stage or cleanup
  was advanced.
- At 2026-09-28T06:13:42+03:00, GTalign refinement advanced again to 94
  checkpoint files: 89 completed, 4 running, and 1 explicit failure. USalign
  tasks 3, 5, 6, and 7 remain active; no dependency stage or cleanup was
  advanced.
- At 2026-09-28T06:15:23+03:00, GTalign refinement advanced again to 95
  checkpoint files: 90 completed, 4 running, and 1 explicit failure. USalign
  tasks 3, 5, 6, and 7 remain active; no dependency stage or cleanup was
  advanced.
- At 2026-09-28T06:16:48+03:00, one GTalign checkpoint transitioned to
  completed without a new file: the ledger is now 95 files with 91 completed,
  3 running, and 1 explicit failure. USalign tasks 3, 5, 6, and 7 remain
  active; no dependency stage or cleanup was advanced.
- At 2026-09-28T06:17:27+03:00, GTalign produced one additional checkpoint;
  the ledger is now 96 files with 91 completed, 4 running, and 1 explicit
  failure. USalign tasks 3, 5, 6, and 7 remain active; no dependency stage or
  cleanup was advanced.
- The corrected-refinement aggregator now preserves nested worker failure
  stage names and error/reason text. This is required for the observed GTalign
  chain-cardinality failure to remain auditable rather than becoming a missing
  row or score zero; the focused aggregation/cleanup/adapter tests pass `5/5`.
- A bounded audit of active USalign raw records found explicit
  `alignment_unavailable` rows alongside valid records with both
  `tm_score_query` and `tm_score_ref`. Updated
  `benchmark/scripts/replay_compact_usalign_batch.py` so compact numerical
  means use complete successful records only, while status/failure counts and
  valid-record counts remain explicit. Focused affected comparison tests pass
  `35/35`; the provider remains valid and does not require rerun. Provenance
  delta: `usalign_production_19855/source_snapshot_after_alignment_summary_contract_fix.json`.
- KUACC reconciliation resolved the prior `3132295_874` disappearance to child
  `3140177` (`COMPLETED`, exit `0:0`) and controller `3132296`
  (`COMPLETED`, exit `0:0`). It scheduled continuation arrays `3140365` and
  `3140401`--`3140407`; the same remote run has 16,651 checkpoint files and
  no aggregate yet. These are the authorized existing common-refinement run;
  no unrelated job was modified and no cleanup was applied.
- Strengthened `benchmark/scripts/aggregate_final_matched_comparison.py`
  before production tables exist: it now reports explicit no-ranking/all,
  deterministic, PRODIGY-when-available, top-1/3/5, and diagnostic-oracle
  strategies, with transformed/refined coverage, cross-interface quality, and
  method-by-split summary tables. Focused reducer tests pass `3/3`; provenance is recorded in the
  final-ranking source-snapshot delta. KUACC refinement reached 16,484
  checkpoints with no aggregate package yet.
- At 2026-09-28T06:22:15+03:00, the corrected GTalign checkpoint probe found
  97 files: 92 completed, 4 running, and 1 explicit input-normalization
  failure. USalign provider tasks `3`, `5`, `6`, and `7` remain active; all
  downstream compact/DockQ/aggregation/handoff jobs remain dependency-held.
  No cleanup or job intervention was performed.
- At 2026-09-28T06:25:50+03:00, KUACC Slurm showed authorized common-refinement
  continuation progress: array `3140365` reached running task `_919`, and
  array `3140401` began running tasks `_786/_787`; later shards remain under
  `AssocMaxJobsLimit`. A large-tree remote checkpoint count probe timed out,
  so no count or completion was inferred and no intervention was performed.
- At 2026-09-28T06:27:19+03:00, GTalign common refinement reached 99 checkpoint
  files: 94 completed, 4 running, and 1 explicit failure. KUACC continuation
  `3140401` reached running task `_862`, while `3140365` remained active at
  `_919`; local USalign and all dependent stages remained active/held.
- At 2026-09-28T06:28:26+03:00, GTalign common refinement reached 100
  checkpoint files: 95 completed, 4 running, and 1 explicit failure. KUACC
  continuation `3140401` reached task `_900`, while `3140365` remained active
  at `_919`; local USalign and downstream stages remained active/held.
- At 2026-09-28T06:29:06+03:00, KUACC continuation `3140401` advanced to
  running task `_919` and `3140365` remained active at `_919`. GTalign stayed
  at 100 checkpoints (95 completed, 4 running, 1 explicit failure); local
  USalign and dependent stages remained active/held.
- At 2026-09-28T06:30:32+03:00, KUACC continuation `3140401` reached running
  task `_964` and `3140365` remained active at `_919`. GTalign remained at
  100 checkpoints (95 completed, 4 running, 1 explicit failure); no local
  compact, DockQ, aggregation, or cleanup marker existed.
- At 2026-09-28T06:31:24+03:00, KUACC continuation `3140401` reached running
  task `_988` while `3140365` remained active at `_919`. Local USalign and
  GTalign workers remained active and downstream jobs remained dependency-held;
  a broad marker scan timed out and was not treated as completion evidence.
- At 2026-09-28T06:32:05+03:00, KUACC continuation `3140401` reached tasks
  `_997`–`_999` and continuation shard `3140402` began tasks `0`–`10`; later
  tasks remained association-limit pending. Local USalign/GTalign workers and
  their downstream dependency chain remained active/held.
- At 2026-09-28T06:32:41+03:00, KUACC refinement shard `3140402` advanced to
  running tasks through `_31`; tasks `_32–999` remained association-limit
  pending. Local USalign/GTalign workers and downstream stages remained
  active/held; GTalign remained at 95 completed, 4 running, 1 failure.
- At 2026-09-28T06:33:43+03:00, KUACC refinement shard `3140402` advanced to
  running task `_65`; tasks `_66–999` remained association-limit pending. Local
  USalign/GTalign workers and downstream stages remained active/held; GTalign
  remained at 95 completed, 4 running, 1 failure.
- At 2026-09-28T06:35:30+03:00, KUACC refinement shard `3140402` advanced to
  running task `_122`; tasks `_123–999` remained association-limit pending.
  Local provider/refinement arrays remained active and GTalign remained at 95
  completed, 4 running, 1 explicit failure; no dependency or cleanup gate
  opened.
- At 2026-09-28T06:34:29+03:00, KUACC refinement shard `3140402` advanced to
  running task `_93`; tasks `_94–999` remained association-limit pending. Local
  USalign/GTalign workers and downstream stages remained active/held; GTalign
  remained at 95 completed, 4 running, 1 failure.
- At 2026-09-28T06:36:08+03:00, KUACC refinement shard `3140402` advanced to
  running task `_146`; tasks `_147–999` remained association-limit pending.
  Local USalign/GTalign workers remained active; GTalign remained at 95
  completed, 4 running, 1 explicit failure.
- At 2026-09-28T06:37:59+03:00, KUACC refinement shard `3140402` advanced to
  running task `_208`; tasks `_209–999` remained association-limit pending.
  Local USalign/GTalign workers remained active; GTalign remained at 95
  completed, 4 running, 1 explicit failure. No downstream gate opened.
- At 2026-09-28T06:38:43+03:00, USalign provider task `1708992_8` started
  while tasks `9–26` remained array-limit pending. KUACC refinement shard
  `3140402` reached task `_234`; GTalign remained at 95 completed, 4 running,
  1 explicit failure and local downstream stages remained dependency-held.
- At 2026-09-28T06:39:34+03:00, USalign task `1708992_8` remained running with
  tasks `9–26` array-limit pending. KUACC refinement shard `3140402` reached
  task `_262`; later tasks remained association-limited. GTalign remained at
  95 completed, 4 running, 1 explicit failure.
- At 2026-09-28T06:40:31+03:00, USalign task `1708992_8` remained active with
  tasks `9–26` array-limit pending. KUACC refinement shard `3140402` reached
  task `_288`; later tasks remained association-limited. GTalign remained at
  95 completed, 4 running, 1 explicit failure.
- At 2026-09-28T06:42:01+03:00, KUACC refinement shard `3140402` advanced to
  running task `_341`; tasks `_342–999` remained association-limit pending.
  USalign task `1708992_8` and all four GTalign workers remained active;
  GTalign remained at 95 completed, 4 running, 1 explicit failure.
- At 2026-09-28T06:42:43+03:00, KUACC refinement shard `3140402` advanced to
  running task `_369`; tasks `_370–999` remained association-limit pending.
  USalign task `1708992_8` and all GTalign workers remained active; GTalign
  remained at 95 completed, 4 running, 1 explicit failure.
- At 2026-09-28T06:44:45+03:00, GTalign common refinement reached 101
  checkpoint files: 96 completed, 4 running, and 1 explicit failure. KUACC
  refinement shard `3140402` reached task `_440`; tasks `_441–999`
  remained association-limit pending and local downstream stages remained held.
- At 2026-09-28T06:43:21+03:00, KUACC refinement shard `3140402` advanced to
  running task `_390`; tasks `_391–999` remained association-limit pending.
  USalign task `1708992_8` and GTalign workers remained active; GTalign
  remained at 95 completed, 4 running, 1 explicit failure.
- At 2026-09-28T06:44:02+03:00, KUACC refinement shard `3140402` advanced to
  running task `_415`; tasks `_416–999` remained association-limit pending.
  USalign task `1708992_8` and GTalign workers remained active; GTalign
  remained at 95 completed, 4 running, 1 explicit failure.
- At 2026-09-28T06:45:40+03:00, KUACC refinement shard `3140402` advanced to
  running task `_471`; tasks `_472–999` remained association-limit pending.
  GTalign remained at 101 checkpoints (96 completed, 4 running, 1 explicit
  failure); USalign task `1708992_8` remained active and downstream local
  stages remained held.
- At 2026-09-28T06:46:49+03:00, KUACC refinement shard `3140402` advanced to
  running task `_511`; tasks `_512–999` remained association-limit pending.
  GTalign remained at 96 completed, 4 running, 1 explicit failure; USalign
  task `1708992_8` remained active and downstream local stages remained held.

- At 2026-09-28T06:58:49+03:00, KUACC refinement shard `3140402` advanced
  through running task `_924`; tasks `_925–999` remained association-limit
  pending, while shards `3140403–3140407` and controller `3140408` remained
  dependency/association-limited. GTalign advanced to 98 completed, 4 running,
  and 1 explicit input-normalization failure. USalign task `1708992_8` remained
  active with tasks `9–26` array-limit pending; downstream USalign stages stayed
  dependency-held.
- At 2026-09-28T06:51:13+03:00, KUACC shard `3140402` had running tasks
  through `_663`, with task `_636` in `COMPLETING`; tasks `_664–999`
  remained association-limited. GTalign remained at 96 completed, 4 running,
  1 explicit failure; USalign task `1708992_8` remained active and downstream
  local stages remained held.
- At 2026-09-28T06:51:57+03:00, KUACC shard `3140402` advanced to running
  task `_694`; tasks `_695–999` remained association-limit pending and the
  previously completing task was no longer listed by Slurm. GTalign remained
  at 96 completed, 4 running, 1 explicit failure; USalign task `1708992_8`
  remained active.
- At 2026-09-28T06:50:37+03:00, KUACC refinement shard `3140402` advanced to
  running task `_647`; tasks `_648–999` remained association-limit pending.
  GTalign remained at 96 completed, 4 running, 1 explicit failure; USalign
  task `1708992_8` remained active and downstream local stages remained held.
- At 2026-09-28T06:50:04+03:00, KUACC refinement shard `3140402` advanced to
  running task `_624`; tasks `_625–999` remained association-limit pending.
  GTalign remained at 96 completed, 4 running, 1 explicit failure; USalign
  task `1708992_8` remained active and downstream local stages remained held.
- At 2026-09-28T06:49:28+03:00, KUACC refinement shard `3140402` advanced to
  running task `_600`; tasks `_601–999` remained association-limit pending.
  GTalign remained at 96 completed, 4 running, 1 explicit failure; USalign
  task `1708992_8` remained active and downstream local stages remained held.
- At 2026-09-28T06:52:34+03:00, KUACC refinement shard `3140402` advanced to
  running task `_714`; tasks `_715–999` remained association-limit pending.
  GTalign remained at 96 completed, 4 running, 1 explicit failure; USalign
  task `1708992_8` remained active and downstream local stages remained held.
- At 2026-09-28T06:47:52+03:00, KUACC refinement shard `3140402` advanced to
  running task `_549`; tasks `_550–999` remained association-limit pending.
  GTalign remained at 96 completed, 4 running, 1 explicit failure; USalign
  task `1708992_8` remained active and downstream local stages remained held.

### Current live execution — 2026-09-28T08:08:17+03:00

- GTalign common refinement checkpoint audit reports 120 completed, 4 running,
  and 2 explicit input-normalization failures. All completed records contain
  `fiberdock`, `external_rosetta`, `dockq_fiberdock`, and `dockq_rosetta`.
- USalign production batches 1–4 have validated completion markers; batches
  5–8 remain active with growing pipeline logs. The compact/DockQ/aggregate
  dependencies remain held until the provider array is complete.
- The existing KUACC refinement lane owns the TMalign/MultiProt selected
  candidates; its current wave is advancing through shard 21. No duplicate
  refinement or cleanup is authorized while consumers remain active.

### Live poll — 2026-09-28T08:09:43+03:00

- GTalign now has 121 completed, 4 running, and 2 explicit failures; completed
  checkpoints still pass required-stage validation.
- USalign provider batches 1–4 are complete with exit markers; batches 5–8
  remain active and continue growing logs.
- KUACC shard 21 is active through task 73 while later work remains
  association-limited; no compact downstream marker is present yet.

### Live poll — 2026-09-28T08:11:01+03:00

- GTalign remains at 121 completed, 4 running, and 2 explicit failures with
  complete stage keys on every completed checkpoint.
- USalign batches 5–8 are still active and their logs grew; batches 1–4 remain
  the only provider batches with validated exit markers.
- KUACC shard 21 advanced through active task 116; no downstream compact
  marker is available and cleanup remains ineligible.

### Live poll — 2026-09-28T08:12:06+03:00

- GTalign remains at 121 completed, 4 running, and 2 explicit failures; active
  logs are changing and completed checkpoint contracts remain valid.
- USalign batches 1–4 remain validated complete, while batches 5–8 are active
  without exit markers.
- KUACC shard 21 progressed through task 142, with task 119 completing; no
  downstream compact marker or cleanup gate is available.

### Live poll — 2026-09-28T08:12:53+03:00

- GTalign remains at 121 completed, 4 running, and 2 explicit failures; active
  logs remain live and completed-stage validation is clean.
- USalign batches 1–4 remain complete and validated; batches 5–8 have no exit
  markers yet.
- KUACC shard 21 advanced through active task 178; downstream aggregation and
  cleanup remain ineligible.

### Contract audit and repair preparation — 2026-09-28T08:20:22+03:00

- GTalign refinement reports 124 completed, 4 running, and 3 explicit
  failures. The new `medium_1ijk_021` failure is a source/native side
  partition mismatch, not a total-chain mismatch.
- The selected manifest contains 268 rows with source-side 2+1 versus native
  1+2 chain partitions. A separate run-scoped repair wrapper preserves source
  partitions and restores the native DockQ mapping; its real-candidate smoke
  test passed and it rejects the known 1wq1 total-chain mismatch.
- The original array remains the sole active owner. Repair submission waits for
  the original failed-row set to become final, preventing overlap.

### Live poll — 2026-09-28T08:21:32+03:00

- GTalign reached 132 completed and 4 running checkpoints. Four failures are
  recorded: one repairable cross-partition candidate and three explicit
  `1wq1` total-chain mismatches.
- The original refinement array remains active; the corrected repair wrapper is
  validated but intentionally not submitted until ownership is released.

### Repair manifest prepared — 2026-09-28T08:24:41+03:00

- GTalign has 134 completed, 4 running, and 4 explicit failures.
- A hashed 268-row cross-partition repair manifest is retained at
  `gtalign_common_refinement_19855/repair/manifest.json`. It is disjoint by
  contract from valid completed rows, but submission waits for array 1709167
  to terminate so the original owner cannot overlap it.

### Live poll — 2026-09-28T08:25:31+03:00

- GTalign reached 135 completed and 4 running checkpoints; the four explicit
  failures remain classified and stage validation is clean for completed rows.
- The original array remains active. The repair manifest is ready but not
  submitted; KUACC shard 21 advanced through task 543 and downstream compact
  jobs remain held.

### Live poll — 2026-09-28T08:26:36+03:00

- GTalign reached 136 completed and 4 running checkpoints; the four explicit
  failures remain unchanged and completed-stage validation is clean.
- Array 1709167 remains active, while KUACC shard 21 advanced through task 581.
  The repair manifest remains ready but unsubmitted and cleanup remains
  ineligible.

### Live poll — 2026-09-28T08:29:25+03:00

- GTalign reached 137 completed and 4 running checkpoints; the four explicit
  failures remain classified and completed-stage validation is clean.
- Array 1709167 remains active, while KUACC shard 21 advanced through task 616.
  The repair manifest remains ready but unsubmitted and cleanup remains
  ineligible.

### Live poll — 2026-09-28T08:30:11+03:00

- GTalign reached 139 completed and 4 running checkpoints; the four explicit
  failures remain classified and completed-stage validation is clean.
- Array 1709167 remains active, while KUACC shard 21 advanced through task 688.
  The repair manifest remains ready but unsubmitted and cleanup remains
  ineligible.

### Live poll — 2026-09-28T08:32:19+03:00

- GTalign reached 141 completed and 4 running checkpoints; all completed
  checkpoints contain the required input-normalization, refinement, and DockQ
  stages, with four explicit failures retained.
- Array 1709167 remains active; the 268-row repair manifest remains ready but
  unsubmitted, and cleanup remains ineligible while downstream consumers are
  held.

### Live poll — 2026-09-28T08:34:11+03:00

- GTalign reached 145 completed and 3 running checkpoints; all completed
  checkpoints pass the required nested-stage audit.
- Failures are now classified as one repairable equal-total cross-partition
  `medium_1ijk_021` candidate and four explicit `medium_1wq1_045` total-chain
  mismatches. Array 1709167 remains active, so repair submission is still held.

### Live poll — 2026-09-28T08:35:13+03:00

- GTalign reached 146 completed and 4 running checkpoints; all completed
  checkpoints pass the required nested-stage audit.
- The five failure classifications remain unchanged, and array 1709167 remains
  active. Repair submission and cleanup remain held.

### Live poll — 2026-09-28T08:36:59+03:00

- GTalign remains at 146 completed and 4 running checkpoints, with five
  classified failures.
- Bounded repair array 1709749 (three shards) was submitted with dependency
  `afterany:1709167:1709178`; it is pending and cannot overlap the original
  refinement or first aggregate. Submission provenance is retained in the run
  package.

### Live poll — 2026-09-28T08:38:03+03:00

- GTalign reached 150 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 is still dependency-held by the original refinement and
  first aggregate; no repair output exists yet. USalign batches 5–8 remain
  active.

### Live poll — 2026-09-28T08:38:53+03:00

- GTalign reached 151 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T08:39:30+03:00

- GTalign reached 152 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T08:40:07+03:00

- GTalign reached 154 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T08:40:46+03:00

- GTalign reached 155 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T08:41:52+03:00

- GTalign reached 157 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T08:42:30+03:00

- GTalign reached 158 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T08:52:48+03:00

- GTalign reached 159 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T08:53:54+03:00

- GTalign reached 160 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T08:54:58+03:00

- GTalign reached 161 completed and 4 running checkpoints; all completed
  checkpoints pass the nested-stage audit and five explicit failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active.

### Live poll — 2026-09-28T09:00:05+03:00

- GTalign reached 166 completed and 4 running checkpoints; all 166 completed
  checkpoint records pass the corrected nested-stage audit and five explicit
  failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active; KUACC refinement also
  remains active.

### Presentation revision — 2026-10-01

- Revised `/scratch/rshadi25/GitHub/PRISM-prescript/docs/PRISM-benchmark.pptx`
  with denominator-aware benchmark groups, self-comparison controls, DockQ/
  Rosetta score interpretation and cross-checks, relaxed-threshold per-case
  counts, and the MultiProt candidate-count explanation.
- Added the 1s78A/1s78D PyMOL input, native, transformed, and refined figures;
  PNG/PSE artifacts are in `tmp/agent/20261001-prism-presentation/pymol/`.
- Final artifact has 40 slides with references last. Python-PPTX bounds/evidence
  QA, wording audit, Python syntax check, and ZIP integrity/uniqueness checks
  pass. No commit or production workload was requested or performed.

### Presentation typography revision — 2026-10-01

- Rebuilt `docs/PRISM-benchmark.pptx` from the preserved pre-revision deck using
  a standardized 24 pt slide-title band, 11.5 pt subtitle band, 8.5 pt source
  footer, compact 15–16.5 pt card text, and run-level font persistence.
- Replaced the generic title slide with `PRISM Docking Benchmark` and the
  subtitle `Pipeline comparison, scoring controls, and refinement-energy analysis`.
- Compactened the DockQ–energy reversal and low-energy/poor-DockQ tables so
  the values remain readable inside the 10 × 5.625 inch canvas.
- Normalized retained historical-slide titles to 24 pt and legacy body text to
  14.5 pt while preserving the original wording and figures.
- Final typography/package audit passes: 40 slides, saved run-level font sizes,
  all slide objects within bounds, expected numerical/text sections present,
  and all seven referenced PyMOL figures available. No production workload or
  source-code change was performed.

### Live poll — 2026-09-28T09:01:36+03:00

- GTalign reached 167 completed and 4 running checkpoints; all 167 completed
  checkpoint records pass the corrected nested-stage audit and five explicit
  failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active; KUACC refinement also
  remains active.

### Live poll — 2026-09-28T09:03:13+03:00

- GTalign reached 169 completed and 4 running checkpoints; all 169 completed
  checkpoint records pass the corrected nested-stage audit and five explicit
  failures remain.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active; KUACC refinement also
  remains active.

### Live poll — 2026-09-28T09:04:01+03:00

- GTalign reached 170 completed and 4 running checkpoints; all 170 completed
  checkpoint records pass the corrected nested-stage audit and five explicit
  failures remain.

### Comparable alignment contract — 2026-09-28

- Added an opt-in `alignment_gate_mode=common_match_coverage` contract to the
  active PRISM source. It applies the same minimum matched-residue count,
  size-adjusted interface coverage, inclusive boundary rule, orientation, and
  clash policy to TMalign, USalign, GTalign, and MultiProt.
- The default `native` mode is unchanged. Comparable mode deliberately does
  not gate on TM-score because MultiProt's stored RMSD-derived proxy is not a
  TMalign-compatible TM-score; provider scores remain available for analysis.
- CLI/environment controls are `--alignment-gate-mode` and
  `PRISM_ALIGNMENT_GATE_MODE`; resolved mode is retained in the audit
  threshold dictionary.
- Focused comparison/configuration tests pass (50), compileall passes, and
  `git diff --check` passes. No production rerun has been submitted for this
  contract yet.
- Repair array 1709749 remains dependency-held by 1709167 and 1709178, with no
  repair output yet. USalign batches 5–8 remain active; KUACC refinement also
  remains active.
