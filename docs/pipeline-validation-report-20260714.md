# PRISM pipeline validation report

## Scope

This report uses the repository benchmark files `T_Rigid.csv`, `T_medium.csv`,
and `T_difficult.csv`: 162 rigid, 60 medium, and 35 difficult rows. The
repository cohort is not the paper’s authoritative 88-case BM3 cohort.

All source identity is keyed by `dataset_row_id`; normalized PDB pairs are
audit fields only. Curated `r_u`/`l_u` files are pipeline inputs and curated
`r_b`/`l_b` files define native truth. Files in `benchmark/data/pdbs` are audit
comparators and are never substituted automatically.

## Current functional status

| Arm | Status | Evidence | Limitation |
|---|---|---|---|
| Root TM-align → current filters → Rosetta | `verified_end_to_end` for smoke | positive `1RGH_B + 1A19_B`, template `1b27AD`; existing 257-row run produced 143 Rosetta models | July full run is observational and unpaired |
| `working_version/TMalign` → historical filters → Rosetta | `verified_end_to_end` for one case | retained `rigid_fixed` model and Rosetta output | no retained native DockQ/iRMSD evaluation for this arm |
| `working_version/multiprot` → FiberDock | `execution_observed; scientifically unscoreable; primary arm blocked` | Explicit reduce.3 substitution plus launcher parent-directory fix produced final `.fiberdock.pdb` and `.intRes.txt` with controller exit 0 | both partner segments use chain `B` with a residue-number reset; standard two-chain scoring returns null DockQ; substitution is not historical-equivalent |
| GTalign → current filters → Rosetta | `exploratory` | parser and alignment integration exist; GTalign 0.19.0 is installed | no retained matched end-to-end speed/quality estimate |

The root and working-version TM-align binaries are byte-identical:
`a42447cdef7ff84f1c7710e57241ebc0644c9d5cb418c81d63411fdcbc82a364`.
Differences between those arms therefore belong to filters, source handling,
configuration, and downstream processing, not the binary itself.

## Cross-pipeline difference ledger

| Dimension | Current TM-align → Rosetta | Historical MultiProt → FiberDock | Evidence/status |
|---|---|---|---|
| Inputs and chains | Curated current `r_u/l_u` role files; current chain-qualified selectors | Derived MultiProt staging; historical selectors and compatibility adapter | The source gate requires curated roles; July outputs used a different observational staging path. |
| Template set | Current calculated-template inventory and current filter path | Historical MultiProt candidate generation and historical filters | Raw candidate streams were not matched; counts are composite evidence only. |
| Tool versions | TM-align 20220412; Rosetta module `rosetta/2022.42` | MultiProt 1.6; FiberDock payload hash recorded; bundled helper runtime incomplete | `runtime_manifest.json`; `references/nprot.2011.367.md`; FiberDock evidence above. |
| Parameters and acceptance | Current TM-score/coverage/clash filters, then Rosetta acceptance | MultiProt matching/filter thresholds, then FiberDock energy/refinement checks | The positive fixture changed from zero to one candidate only under the recorded 40% diagnostic intervention. |
| Intermediate outputs | Alignment JSON, transformations, Rosetta models | MultiProt `.res` files, transformed candidates, FiberDock energy/intermediate files | Per-task evidence is retained in the legacy smoke/probe directories. |
| Prediction counts | 143 Rosetta model rows in the July observational run | 2 FiberDock model rows, for one pair, in the July observational run | No shared scoreable pair; these are not paired quality estimates. |
| DockQ and iRMSD | Existing reports contain scores for some current models | Retained FiberDock output is not scoreable as a two-chain complex; the historical Rosetta model and chain-fixed paths named by the old CSV are unavailable here | Paired DockQ/iRMSD difference: not estimable; no-model metrics must remain null until chain-preserving artifacts are recovered or regenerated. |
| Interface size | Not comparable from the unpaired July outputs | Not comparable from the unpaired July outputs | Must be recomputed from the same frozen native mapping in the confirmatory pilot. |
| Runtime and throughput | Batch-level logs exist; no per-pair matched runtime estimate was frozen | Batch-level logs exist; helper failure changes terminal work | Runtime comparison remains unresolved until isolated one-task-per-pair records are replayed. |

The existing July aggregate means and raw prediction counts therefore document
an implementation difference, not a causal accuracy difference. The first
valid metric comparison requires the same curated inputs, template panel,
top-20 budget, native-independent ranking, and evaluator for both arms.

## FiberDock evidence

The checked-out executable is 64-bit and has hash
`5dc7408b0f370ff175abd924f9732c38ce4a017b0e905d4a5a8e45b1d3fe2d13`.
The direct probe used energy-only settings (`onlyEnergyCalculation=1`, no
backbone refinement, no refined-complex output), exited 0, and wrote a
nonempty `resFile.ref`. This confirms component loading and energy-only
behavior only.

The integrated positive diagnostic used `1RGH_B + 1A19_B / 1b27AD` with a
recorded 40% threshold intervention. It generated one candidate and reached
FiberDock intermediates, but returned
`pipeline_blocked_full_refinement_capability` with no final refined PDB.
`nma` ran successfully; the configured dynamic `reduce.2.21.030604` failed
because the host lacks its 32-bit loader/runtime. `FiberDock.32` and
`reduce.3` were not used and are not silently substituted.

The remaining validation experiment is an unchanged-helper microprobe for
`reduce.2`, followed by a complete bundled-example chain and one-controller
replay. Full FiberDock status remains blocked until nonempty hydrogenated
inputs, meaningful energies, final PDB output, and interface output are all
present with exit code 0.

### Reduce.3 compatibility experiment (2026-07-14)

The unchanged helper path was not available on the host: the bundled
`reduce.2.21.030604` is a dynamically linked 32-bit ELF executable whose
required loader and `libstdc++.so.5` were not present. The derived staging
script now accepts `--fiberdock-reduce-helper` for an explicitly labelled
exploratory substitution; it records the original and effective hashes and
keeps the environment capability flag fail-closed.

The experiment used the same positive diagnostic fixture (`1RGH_B + 1A19_B`,
template `1b27AD`, diagnostic MultiProt threshold 40%) and an isolated
workspace whose `external_tools` points to the derived environment. The first
attempt was invalid as a substitution test because the workspace still
resolved its copied helper; its stderr is retained under
`tmp/agent/20260714-fiberdock-reduce3-smoke/task-2/retry-1355540-r0/`.
The corrected attempt was KUTEM job `1355544`, which reached non-empty `.HB`,
NMA, `.fib`, and `.ref` artifacts but stopped at the legacy missing-parent
directory call. The launcher-side `mkdir -p fiberdock_output/smoke` fix was
then tested by job `1355545`: controller exit code `0`, final
`1b27AD_pdb1_0_pdb2_0.fiberdock.pdb` (131462 bytes), final
`.intRes.txt` (1115 bytes), and non-empty hydrogenated/NMA/intermediate files.
The effective helper hash is
`6d066f88bff740627d7c1d2fb0200d326978fe70a0c041a2528c8853e682b1ce`, while
the original reduce.2 hash is retained in the environment manifest as
`1f12c9d3931d95dc549f0bb7e140c9462764f13b575241d5f24c6d5b3e1e115f`.

This confirms that the derived adapter can execute the complete legacy
controller path for one diagnostic pair. It does not confirm that reduce.3 is
algorithmically or numerically equivalent to reduce.2, and it does not open
the primary FiberDock arm. The official Reduce repository describes `reduce2`
as the maintained successor to the original Reduce but does not document the
historical PRISM helper path as a reduce.3 drop-in replacement. PRISM’s
protocol documents FiberDock as the flexible-refinement and energy-ranking
stage, and the BADock repository documents FiberDock as a wrapper-driven
refinement stage. These sources support the conservative separation between
execution success and historical-equivalence claims:
<https://github.com/rlabduke/reduce>,
<https://pmc.ncbi.nlm.nih.gov/articles/PMC7384353/>, and
<https://github.com/badocksbi/BADock>.

The final PDB is nevertheless invalid for the frozen evaluator: all 1962 raw
ATOM records use chain `B`, and the residue numbering resets from 96 to 1 for
the second partner. Biopython therefore reads one chain with 96 residues and
974 atoms rather than two partners. The standard scorer returns
`dockq=null`, `dockq_global=null`, and an invalid `BB:AB` mapping. The old CSV
row records DockQ `0.8541875620` and iRMSD `0.612`, but both its historical
model path and its `rosetta_output_1_chainfixed` path are absent from the
current workspace; those values cannot validate the retained FiberDock PDB
or be independently rescored here, so they are excluded from the new
comparison.

### Exploratory chain-preserving repair probe

The explicit diagnostic repair `benchmark/scripts/split_legacy_partner_pdb.py`
split the retained FiberDock file at its single residue-number reset and
rewrote the two segments as chains `A` and `B`. The derived record is retained
under `tmp/agent/20260714-fiberdock-reduce3-smoke-parentfix/chain-repair-probe/`
with source hash `cabce7b1a9b6bd670bafee56c8d0049f056fc8198a246f3913d66f77ffd2a590`.
Biopython then reports `(A: 96 residues, 974 atoms)` and `(B: 89 residues,
988 atoms)`, and the frozen scorer accepts the mapping. The exploratory
derived file scores DockQ `0.7447319`, DockQ iRMSD `1.2151`, grouped iRMSD
`1.204`, interface Fnat `0.8`, and two clashes.

This is evidence that a chain-preserving output/postprocessing path can make
the file scoreable; it is not evidence that the split recovers the historical
Rosetta model, the intended FiberDock partner order, or reduce.2 equivalence.
The result is diagnostic only and is excluded from benchmark estimates. The
next repair experiment must provide chain-distinct inputs to FiberDock and
compare native output with this postprocessed view before any batch run.

Confirmed: the launcher parent-directory defect was fixed; reduce.3 works in
the tested helper invocation; one isolated controller run emitted final
FiberDock artifacts with exit code 0. The output-integrity gate also confirms
that those artifacts are not scientifically scoreable.

Likely: the remaining primary blockers are output-chain preservation and
historical helper equivalence, because the substituted run completed the
expected stages but failed the two-chain evaluator contract.

Unresolved: whether the chain collision is introduced by staging, FiberDock,
or the legacy output writer; whether the unavailable chain-fixed artifact can
be recovered; energy/ranking equivalence; and behavior across the benchmark
cohort. These require a chain-preserving replay and matched reduce.2/reduce.3
fixtures or a permitted native reduce.2 runtime.

## GTalign comparison status

GTalign 0.19.0 has CPU and GPU payloads. Existing large GTalign output trees
are stale/unprovenanced relative to the current 946 fully resolvable template
assets and cannot support a quality claim. The new harness records both
query- and reference-normalized TM-scores and unique run directories. Raw
alignment-output hashes are not yet attached to every retained pair record,
so the alignment artifacts do not support a stronger provenance claim.

The confirmatory comparison will use a deterministic 12-row pilot, followed
by the 240 strict-valid rows after all source/evaluator gates pass. The
practical noninferiority margins are a maximum loss of 0.02 best-GlobalDockQ
at 20, 5 percentage points model-success rate, and 0.5 Å grouped iRMSD.

## Confirmed, likely, unresolved

## Final confirmation evidence

- Fresh environment validation on `ai11` imported Python 3.11.13, Biopython
  1.84, NumPy 1.26.4, and pandas 2.3.3. `gtalign_cpu -h` reported GTalign
  0.19.00. The runtime manifest and all frozen executable hashes matched.
- The final focused verification suite passed `45` tests. The retained
  current TM-align smoke completed all stages and wrote alignment,
  transformation, and Rosetta output directories; it reported `Passed pairs
  0` for that plumbing fixture. The corresponding GTalign run also exited 0
  and wrote isolated records containing both normalized TM-scores and a raw
  output hash.
- Retained KUTEM evidence is pair/task isolated. Compatibility tool probes
  `1355097` passed four tools; the final historical probe array `1355128`
  passed NACCESS, POPS, MultiProt, and FiberDock. The earlier historical v5
  NACCESS task failed with `libgfortran.so.3` missing and is not treated as
  the final repaired result.
- The historical positive row for `1RGH_B + 1A19_B` with template `1b27AD`
  records DockQ `0.8541875620` and iRMSD `0.612`, but its referenced
  `rosetta_output_1_chainfixed` PDB is unavailable. Re-scoring the retained
  raw Rosetta and FiberDock PDBs fails the two-chain contract and returns
  null DockQ; the historical values are therefore not accepted as validation
  of the retained files.
- Slurm controller access was unavailable from the current `ai11` session;
  no new job was submitted during this confirmation. Existing task-level
  `exit.json`, logs, hashes, and requested-resource records remain the
  authoritative KUTEM evidence.

Confirmed:

- The benchmark contains 257 repository rows, not the paper’s 88-case list.
- The 17 source-gate chain-contract disagreements remain unresolved.
- Root and working TM-align executables are identical.
- Existing July method-level results have no shared scoreable pair.
- FiberDock energy-only execution works; a reduce.3-substituted controller run
  completed, but historical full-refinement capability and scoreable output
  remain unverified.
- The retained FiberDock final PDB fails the raw two-chain output contract;
  final artifact existence is not a valid scoreability result. The historical
  Rosetta artifacts required to test the same contract are unavailable.
- An explicit boundary-split repair makes the FiberDock file scoreable
  (DockQ `0.7447319`) but remains exploratory because partner provenance and
  tool equivalence are not established.
- GTalign has not yet been shown faster without quality loss.

Likely explanations requiring matched experiments:

- Candidate multiplicity, filter thresholds, template inventories, and refiner
  acceptance rules explain part of the current-versus-historical count gap.
- GTalign batching may improve alignment throughput, but this is not yet an
  end-to-end result.

Unresolved:

- Authoritative orientations for 17 rows.
- Availability of a compatible exact `reduce.2` runtime.
- Origin and repair of the duplicated-chain FiberDock/Rosetta output, or
  recovery of the historical chain-fixed artifact.
- Whether chain-distinct FiberDock inputs preserve the same partner boundary
  without postprocessing.
- GTalign speed/quality outcome under common inputs and Rosetta refinement.
- Whether any FiberDock replacement is scientifically equivalent.

## Final benchmark batch replay and comparison (2026-07-14)

The requested final batch was run as an observational replay of the retained
benchmark outputs, not as a new confirmatory model-generation run. The input
manifest is the existing 257-row repository benchmark manifest; the replay
contains 145 retained model rows (143 current TM-align/Rosetta rows and 2
legacy MultiProt/FiberDock rows). Missing source cases and absent historical
chain-fixed artifacts were not silently regenerated or substituted.

Each array index processed one shard in its own task directory. The strict
replay used KUTEM array job `1355946`; the alignment-enabled secondary audit
used job `1355947`. Both used the requested 1-node, 1-task, 2-CPU, 2-GB,
5-minute profile. The scoring environment was explicitly separated from the
pipeline environment: `gtalign_env` ran the wrapper and
`/scratch/tmp/prism-dockq-env/bin/python` supplied DockQ.

### Results against the previous report

| Arm | Previous report | Strict no-align replay | Alignment-enabled audit |
|---|---:|---:|---:|
| MultiProt/FiberDock model rows | 2 | 2; DockQ n=0; iRMSD n=0; errors=2 | DockQ n=2, mean 0.726814, best 0.741588; iRMSD n=2, mean 1.2585, best (min) 1.254; errors=0 |
| TM-align/Rosetta model rows | 143 | 143; DockQ n=0; iRMSD n=0; errors=143 | DockQ n=119, mean 0.115281, median 0.018902, best 0.644710; iRMSD n=127, mean 19.7861, median 18.946, best (min) 1.521; errors=29 |

The previous current-arm report listed DockQ n=103, mean 0.0436796,
median 0.012, best 0.785 and iRMSD n=128, mean 19.8055, median 18.9975,
with 101.335 reported as “best.” The hardened collector recomputes best
iRMSD as the minimum: 1.521 for the previous current rows and 1.254 for
the previous legacy rows. The previous legacy report listed DockQ n=2, mean
0.727, best 0.742 and iRMSD n=2, mean 1.2585. The alignment-enabled
audit reproduces the legacy values to rounding and exactly reproduces the
legacy iRMSD summary, but it is not a like-for-like evaluator comparison for
the current arm because the previous current report used a different
no-align/mapping path. The current DockQ increase is therefore an evaluator
configuration difference, not evidence of improved pipeline quality.

The strict replay is the valid confirmatory evaluator contract. It returned
null structural metrics because all 145 retained models lacked the complete,
explicit residue correspondence required for `--no-align`; this is a
fail-closed mapping result, not a claim that every model has zero quality.
The alignment-enabled run is retained as a diagnostic audit only. Its score
intersection with the previous aligned set is 50 pairs; it has no current-only
scoreable pair relative to that prior set and one prior-only row (`rigid_0122`).

Prediction counts were unchanged because this batch rescored retained models:
143 current rows (137 ready) and 2 legacy rows (2 ready). Strict scoring
produced no valid interface-size or structural-metric denominator; the aligned
audit was used only for DockQ/iRMSD reconciliation and does not establish a
new interface-size estimate. In the hardened rerun, strict shard runtimes were
18--22 seconds and aligned shard runtimes were 44--72 seconds, with all tasks completed under the
5-minute/2-GB contract. These are scoring-wrapper runtimes, not alignment,
refinement, or end-to-end pipeline throughput, so they cannot support a
GTalign speed claim.

Replay artifacts:

- Strict: `tmp/agent/20260714-observational-score-replay-strict-v4/collected/COMPARISON.md`
  and `replay_summary.json` (job `1355946`).
- Alignment-enabled audit:
  `tmp/agent/20260714-observational-score-replay-aligned-v4/collected/COMPARISON.md`
  and `replay_summary.json` (job `1355947`).
- Per-shard inputs, command files, parameters, logs, outputs, and exit records
  are under the corresponding `tasks/task-*` directories.

### Final status

Confirmed: the array batch completed all ten shards in both evaluator modes;
the wrapper used separate task directories; the legacy aligned scores are
reproducible; strict no-align scoring correctly prevents unsupported metric
promotion; and the previous current-arm mean cannot be compared directly to
the new aligned mean.

Likely: the remaining current-arm score differences are caused by evaluator
mapping/alignment behavior and historical output provenance, not by an
isolated change in TM-align or Rosetta. This remains a likely explanation
because no matched one-factor pipeline-generation experiment has been run.

Unresolved: whether the 143 current models can be regenerated with complete
chain/residue correspondence; whether the missing chain-fixed historical
artifacts can be recovered; and whether a full, source-gated benchmark run
with a functional historical FiberDock refinement can produce a causal
pipeline comparison. No accuracy, interface-size, or end-to-end throughput
claim is accepted from this observational replay.

## PyRosetta and GTalign isolated comparison (2026-07-14)

### PyRosetta arm

The PyRosetta arm is implemented separately from the stable external-Rosetta
pipeline. `src/pyrosetta_refinement.py` lazy-loads PyRosetta, records input and
output hashes, rejects stale outputs, and publishes a refined pose only after
successful completion. `benchmark/scripts/run_pyrosetta_refinement.py` supports
one-pose and manifest modes; `benchmark/jobs/pyrosetta_refinement_array.sbatch`
provides isolated KUTEM tasks.

The `gtalign_env` probe reported `ModuleNotFoundError: No module named
'pyrosetta'`. A one-pose smoke using the positive `1b27AD`/`1rghB`/`1a19B`
fixture returned code 3 with status `unavailable`, preserved the input SHA256,
and did not create a refined PDB. An authorized wheel was not available locally:
the official West mirror exposed a 1.659 GB cp311 wheel but the download was
cancelled before completion at approximately 203 kB/s; the East mirror failed
TLS certificate validation and a bounded retry timed out. No PyRosetta
environment was created or modified. See [PyRosetta downloads](https://www.pyrosetta.org/downloads)
and [PyRosetta licensing](https://www.pyrosetta.org/home/licensing-pyrosetta).

Confirmed: the separate arm is fail-closed and ready for a licensed,
hash-pinned wheel, but no PyRosetta refinement result exists and no PyRosetta
quality comparison is valid.

Likely: the installation blocker is package staging/network access, not an
adapter or input-structure failure, because the import gate fails before
refinement and no output is published.

Unresolved: whether an authorized wheel imports and initializes on this cluster,
and whether PyRosetta produces a chain-preserving pose with comparable DockQ,
iRMSD, interface size, runtime, and prediction counts.

### GTalign versus TM-align

Hardened KUTEM job `1356100` used the CPU GTalign binary because the GPU binary
cannot run on a CPU-only KUTEM node. Both tools processed the same real PRISM
inputs (`1fgnH` against template `1kcaCH`, chains C and H), producing two paired
records each and no missing pair. TM-align took 0.0272 s; GTalign CPU took
1.8055 s. Mean absolute TM-score difference was 0.01811 (90th percentile
0.03076), and mean absolute match-count difference was 0.5. This job stopped at
the alignment boundary; it did not produce transformed poses, refinements,
DockQ, iRMSD, interface-size, or end-to-end prediction-count measurements.

| Pilot | Pair workload | TM-align | GTalign CPU | Interpretation |
|---|---:|---:|---:|---|
| 10x10 | 100 pairs, 50 residues | 2.349 s | 8.115 s | TM-align faster on this small pilot |
| 50x50 | 2500 pairs, 50 residues | 47.019 s / 53.17 pairs/s | 46.543 s / 53.71 pairs/s | Essentially parity; GTalign ~1.01x faster |

These pilots measure alignment-stage wall time and output production only. They
do not support an accuracy claim. In the 50x50 pilot, GTalign wrote 50 query
output files, each containing 50 ranked reference records, for 2,500 parsed
records total; therefore 53.71 pairs/s is an effective all-pairs-searched rate,
not an output-file rate. The real smoke confirms that GTalign can be
used as an operational alignment backend on KUTEM, but does not show that it
improves the PRISM pipeline or preserves downstream quality.

Confirmed: GTalign CPU is functional for the tested synthetic batches and the
real two-chain alignment smoke; the GPU binary failed on KUTEM with
`cudaErrorInsufficientDriver`, so it was not used for the CPU comparison.

Likely: GTalign batching may help for larger workloads, but the available
evidence supports only alignment-throughput parity at 2500 synthetic pairs.

Unresolved: downstream transformation/refinement compatibility, multichain pose
quality, DockQ/iRMSD/interface-size parity, failure rate on the 257-row cohort,
and whether a V100 GPU run changes the speed conclusion. No claim that GTalign
improves overall accuracy is accepted.

## PyRosetta installation update (2026-07-15)

PyRosetta is now installed and importable in `/home/rshadi25/.conda/envs/gtalign_env`.
The verified distribution is `2026.3+releasequarterly.5e498f1409`, using Python
3.11.13. `pyrosetta.init('-mute all')` completed successfully, and the module
path is:

`/home/rshadi25/.conda/envs/gtalign_env/lib/python3.11/site-packages/pyrosetta/__init__.py`

The first refinement smoke exposed two installed-version API differences in the
adapter: `set_docking_local_refine` requires a boolean argument, and
`DockingProtocol` uses `set_highres_scorefxn` rather than `set_scorefxn`. The
adapter now supports these APIs while retaining compatibility with the test
double and older setter form.

The corrected one-pose smoke completed successfully using the positive
`1b27AD`/`1rghB`/`1a19B` assembled pose:

- status: `success`
- runner return code: `0`
- PyRosetta total score: `290.7085762721651`
- output SHA256: `043def6cfc9b7bbf585ecd0c84ae451da796803f9a7f40ae48f4948ba060e449`
- output: `tmp/agent/20260715-pyrosetta-installed/one-pose-v3/refined.pdb`

This verifies installation, initialization, protocol construction, refinement,
and output publication for one pose. It does not yet validate benchmark-scale
quality, chain preservation across all rows, DockQ/iRMSD parity, or runtime
against external Rosetta.

Focused PyRosetta tests pass: `6 passed`. The PyRosetta arm is therefore
`installed; one-pose verified; benchmark comparison pending`.

The default `-mute all` refinement is stochastic: independent default smokes
returned total scores `290.7085762721651` and `288.687558768033`. Two seeded
smokes using `-mute all -constant_seed -jran 12345` both returned
`294.41123996335625`. Their coordinates and scores agree, but their PDB byte
hashes differ because PyRosetta embeds the run-specific temporary output path
in the pose-energy-table comments. Benchmark runs must therefore freeze the
random-seed options and either normalize that comment before byte-level
comparison or compare structural/score content separately.
The rerun copied the immutable smoke seed into the unique workspace
`tmp/agent/prism-pipeline-smoke-1356100` and wrote outputs under
`tmp/agent/20260714-pyrosetta-gtalign-comparison/prism-smoke-1356100`.
