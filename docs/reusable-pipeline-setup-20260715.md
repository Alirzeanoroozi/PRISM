# Reusable PRISM pipeline setup and status

Last verified: 2026-07-15.  This guide records executable paths and evidence
for reuse; it does **not** reopen the source or evaluator gates for a
confirmatory benchmark claim.

## Shared setup

Run from the repository root:

```bash
cd /scratch/rshadi25/GitHub/PRISM-prescript
source /opt/ohpc/pub/compiler/conda3/latest/etc/profile.d/conda.sh
conda activate /home/rshadi25/.conda/envs/gtalign_env
export PRISM_ROOT=$PWD
export PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python
export DOCKQ_PYTHON=$PRISM_ROOT/benchmark/prism_processed/env/prism_score_env/bin/python
```

`gtalign_env` is the primary Python 3.11 runtime.  Its observed versions are
Python 3.11.13, NumPy 1.26.4, pandas 2.3.3, Biopython 1.84, and PyRosetta
`2026.3+releasequarterly.5e498f1409`.  DockQ is available through the
repository-local scoring prefix at
`benchmark/prism_processed/env/prism_score_env/bin/python` (Python 3.9.23,
NumPy 1.26.4, DockQ 2.1.3). Invoke it as `"$DOCKQ_PYTHON" -m DockQ`; the
prefix's standalone `bin/DockQ` launcher has a stale historical shebang.
The `environment.yaml` recipe still declares Python 3.11.13, so the exact
release scoring identity must be reconciled before a new benchmark batch.
`environment.yaml` pins the non-licensed base dependencies; PyRosetta remains
an authorized separately installed package, described by
`environment-pyrosetta.yaml`.

Important roots:

| Purpose | Path |
| --- | --- |
| Current driver and source | `prism.py`, `src/` |
| Current TM-align | `external_tools/TMalign` |
| GTalign CPU | `/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu` |
| GTalign GPU (not validated) | `/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign` |
| Current templates | `templates/interfaces/`, `templates/interfaces_lists/` |
| Current processed structures | `processed/pdbs/` |
| BM5/5.5 extension tables | `benchmark/data/T_Rigid.csv`, `T_medium.csv`, `T_difficult.csv` |
| Cached audit PDBs | `benchmark/data/pdbs/` |
| Curated archive/source gates | `benchmark/`, `tmp/agent/20260713-investigation-implementation/source-gate-aggregate-final/` |
| Legacy working copy | `working_version/multiprot/` |
| FiberDock payload | `working_version/multiprot/external_tools/fiberdock/` |
| Durable validation report | `docs/pipeline-validation-report-20260714.md` |
| Reusable output root | `tmp/agent/<YYYYMMDD>-<purpose>/` |

Never replace curated role files with the full-PDB cache.  The repository
tables contain 257 BM5/5.5-extension rows (162 rigid, 60 medium, 35 difficult),
not the paper's authoritative 88-case BM3 cohort.

## Pipeline status and commands

| Pipeline / scope | Status | Verified boundary | Not established |
| --- | --- | --- | --- |
| Current TM-align → external Rosetta | **Partially working** | Isolated plumbing/smoke and previous positive diagnostic | Source-gated 257-row causal benchmark |
| PyRosetta refinement arm | **Partially working** | Import/init and one canonical-pose refinement | Matched benchmark quality and chain-preservation coverage |
| GTalign CPU alignment arm | **Partially working** | Real paired alignment smoke and synthetic throughput pilots | PRISM downstream quality/replacement claim |
| DockQ replay/evaluator | **Ready and verified** for contract-preserving scoring replay | DockQ JSON normalization, raw-chain gate, provenance collection | A fair current-versus-legacy method estimate |
| Legacy MultiProt standalone | **Ready and verified** for isolated MultiProt probes | Python 2.7.15 compatibility staging and KUTEM tool probes | Historical full pipeline equivalence |
| Legacy MultiProt → FiberDock | **Not working** as a primary scientific arm | Energy-only/experimental controller execution | Native full refinement and scoreable two-chain output |
| Exact-paper BM3 88-case reproduction | **Not yet tested** | None; source is unavailable | Authoritative cohort reconstruction |

### Current TM-align → external Rosetta

Inputs are five-character chain-qualified target IDs, a staged template list,
the `processed/pdbs/` structures, and current template assets.  Expected
outputs are task-local `processed/alignment/`, `processed/transformation/`,
`processed/rosetta_refinement/`, and `run.log` directories.

Local plumbing check:

```bash
PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  bash benchmark/scripts/run_prism_pipeline_smoke.sh
```

For Rosetta refinement in Slurm, load the validated module in the job script:

```bash
source /opt/ohpc/admin/lmod/lmod/init/bash
module load rosetta/2022.42
cd <isolated-workspace>
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py 2>&1 | tee run.log
```

Production defaults are TM-score 0.5, 15 matched residues, 50% matched
residues, 20% difference allowance, 3.0 A clash distance, 5 maximum clashes,
and scaffold compatibility 5.0.  Environment threshold overrides are
diagnostic only.  A completed zero-pair smoke verifies plumbing, not docking
quality.

### Optional PyRosetta refinement

The adapter is `src/pyrosetta_refinement.py`; the runner is
`benchmark/scripts/run_pyrosetta_refinement.py`.  PyRosetta needs no GPU for
the verified CPU smoke.  Freeze the random seed for every benchmark task:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  benchmark/scripts/probe_pyrosetta_environment.py \
  --output tmp/agent/<run-id>/pyrosetta-probe.json

/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  benchmark/scripts/run_pyrosetta_refinement.py \
  --input-pdb benchmark/jobs/joblist_014/flexibleRefinement/1b27AD_1rghB_0_1a19B_0.rosetta.pdb \
  --output-pdb tmp/agent/<run-id>/pyrosetta/refined.pdb \
  --partners A_B \
  --init-options '-mute all -constant_seed -jran 12345' \
  --json tmp/agent/<run-id>/pyrosetta/result.json \
  --require-available
```

Expected outputs: `refined.pdb`, `refined.pdb.pyrosetta.json`, and result JSON
with `status: success`.  The one-pose smoke succeeded.  Two seeded runs gave
the same score (`294.41123996335625`); byte hashes can differ only because
PyRosetta embeds the temporary output path in PDB energy-table comments.

For up to ten independent manifest rows, use the isolated KUTEM array.  The
manifest must have `input_pdb` and optional `partners` columns; each array task
writes its own `task-N/` directory with command, parameters, result, logs, and
`exit.json`.

```bash
RUN_ROOT=$PRISM_ROOT/tmp/agent/<run-id>/pyrosetta-array \
TASK_MANIFEST=$PRISM_ROOT/<manifest.csv> \
PYROSETTA_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
PYROSETTA_INIT_OPTIONS='-mute all -constant_seed -jran 12345' \
sbatch benchmark/jobs/pyrosetta_refinement_array.sbatch
```

### GTalign CPU comparison arm

Use the CPU executable on KUTEM; the GPU binary has not been validated and
fails in CPU-only allocations.  The real PRISM alignment smoke produces a
task-local comparison summary, raw outputs, parameters, and `exit.json`:

```bash
RUN_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
sbatch benchmark/jobs/prism_gtalign_smoke_compare.sbatch
```

The verified command is `gtalign_cpu -h`, reporting version 0.19.00.  In the
real `1fgnH`/`1kcaCH` KUTEM smoke (job 1356100), both paired chain records were
produced; TM-align took 0.0310 s and GTalign CPU 1.8042 s.  Synthetic 50x50
alignment was throughput-parity only (53.17 vs 53.71 searched pairs/s).  No
DockQ, iRMSD, interface-size, failure-rate, or full-pipeline accuracy claim is
supported yet.

For a synthetic, isolated runtime pilot:

```bash
RUN_ROOT=$PRISM_ROOT/tmp/agent/<run-id>/alignment-pilot \
RUN_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
QUERIES=50 REFS=50 RESIDUES=50 TMALIGN_LIMIT_PAIRS=2500 \
sbatch benchmark/jobs/alignment_backend_comparison.sbatch
```

### DockQ and iRMSD scoring/replay

Use the scorer only with declared model/native receptor and ligand chains.
Invalid raw PDB partner structure is rejected before structural scores are
calculated; null metrics are intentional, not zeroes.

```bash
$DOCKQ_PYTHON benchmark/scripts/score_single_prism_pair.py \
  <model.pdb-or-directory> <native.pdb> \
  --model-receptor A --model-ligand B \
  --native-receptor <native-receptor-chains> \
  --native-ligand <native-ligand-chains> \
  --score-python $DOCKQ_PYTHON \
  --dockq-json-dir tmp/agent/<run-id>/dockq-json \
  --out-csv tmp/agent/<run-id>/scores.csv
```

The retained strict and alignment-enabled replay summaries are under
`tmp/agent/20260714-observational-score-replay-{strict,aligned}-v4/collected/`.
They are validity audits, not paired method estimates: 143 current and two
legacy model rows had no shared scoreable pair, and evaluator regimes differ.

### Legacy MultiProt and FiberDock

The derived, non-destructive MultiProt setup is under
`tmp/agent/20260713-investigation-implementation/multiprot-environment-v4/`.
When reusing it, activate only the derived copy:

```bash
source tmp/agent/20260713-investigation-implementation/multiprot-environment-v4/activate.sh
"$PRISM_MULTIPROT_PYTHON" --version
"$PRISM_MULTIPROT_ENV/bin/multiprot" --version
```

For a fresh reproducible stage rather than reusing an old derived result:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  benchmark/scripts/install_multiprot_environment.py \
  --source-root working_version/multiprot/external_tools/multiprot \
  --output-root tmp/agent/<run-id>/multiprot-environment \
  --python-executable /home/rshadi25/.conda/envs/tmalignRosetta/bin/python2.7
```

FiberDock's checked-in payload includes `FiberDock`, `FiberDock.32`, `nma`,
`reduce.2.21.030604`, and `reduce.3.23.130521`.  Native full refinement is
blocked: the required helpers are 32-bit and the compatible loader plus
`libstdc++.so.5` is absent.  An explicit reduce.3 substitution can reach an
exploratory final PDB but is not equivalent to historical reduce.2 and its
retained output collapses both partners to chain B with a residue reset.
Therefore do not use FiberDock PDBs for primary DockQ/iRMSD/interface claims.

## Batch contract, validation, and cleanup

All KUTEM arrays use the approved 2 CPU, 2 GB, 5-minute, 1–10 task profile.
Set a new `RUN_ROOT` for each submission; never share workspace/output paths
between tasks.  Each task must retain input hashes, parameters, command,
stdout/stderr, outputs, runtime, and `exit.json`.  Aggregate only after all
child tasks terminate; Slurm completion is not scientific success.

Focused validation after the PyRosetta API/default-array fixes:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_pyrosetta_refinement.py \
  tests/test_probe_pyrosetta_environment.py \
  tests/test_runtime_and_gtalign_contracts.py \
  tests/test_observational_score_replay.py
```

The last focused run passed 14 tests.  Preserve all source, benchmark inputs,
curated archives, validated reports, environments, and provenance logs.  Do
not delete `tmp/agent/` evidence without a reviewed cleanup manifest; the
existing `tmp/agent/20260714-pipeline-validation/cleanup_manifest.tsv` records
a prior no-delete decision.

## Confirmed limitations and next experiment

Confirmed: the exact BM3 88-case source list is unavailable; 17 BM5/5.5 rows
have unresolved curated chain-contract/orientation disagreement; and the
FiberDock primary arm is blocked by runtime and chain-integrity evidence.

Likely explanation for the old/current report difference: unequal staging,
filters, output integrity, and evaluator regimes dominate the observed
aggregate gap.  That is not a causal aligner or refiner conclusion.

Recommended next step: obtain an authoritative decision for the 17 source-gate
rows, then run one matched, source-gated shard with the same canonical poses
through external Rosetta and seeded PyRosetta, followed by GTalign versus
TM-align candidate/filter replay and frozen DockQ evaluation.
