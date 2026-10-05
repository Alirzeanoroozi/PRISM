# Continuation Handoff: TM-align Biological Ranking

Updated: 2026-07-12
Repository: `/scratch/rshadi25/GitHub/PRISM-prescript`

## Objective

Improve the biological usefulness of the current TM-align + Rosetta pipeline
and implement validated learning-based candidate ranking. Keep the legacy
MultiProt + FiberDock pipeline as a separate comparison baseline. Do not claim
ML improvement without independent native complexes and grouped held-out
evaluation.

## Pipeline Versions

### Current pipeline

- Entry point: `prism.py`
- Alignment: `src/alignment.py` with `external_tools/TMalign`
- Surface extraction: `src/surface_extract.py`, default NACCESS
- Optional surface backend: FreeSASA via `--surface_backend freesasa`
- Transformation/filtering: `src/transformation.py`
- Refinement: `src/rosetta_refinement.py`
- Smoke runner: `benchmark/scripts/run_prism_pipeline_smoke.sh`
- Default alignment backend: TMalign
- Experimental alternative: GTalign; SoftAlign remains experimental and is
  excluded from the current-vs-legacy comparison unless explicitly requested.

### Legacy pipeline

- Repository: `/scratch/rshadi25/GitHub/prism-oldversion`
- Alignment: MultiProt
- Refinement: FiberDock
- Surface selection is configured in legacy `prism.ini`.
- Do not mix legacy outputs into current-pipeline ML training data.

## Environments

- Current pipeline/local verification:
  `/home/rshadi25/.conda/envs/gtalign_env/bin/python`
  - Python 3.x, pandas/numpy/Bio available
  - scikit-learn and PyTorch were verified available at different checks
- DockQ scoring environment:
  `/scratch/tmp/prism-dockq-env/bin/python`
  - Python 3.10
  - DockQ 2.1.3 installed and verified
- Current test environment:
  `/scratch/tmp/prism-current-test-py311/bin/python`
  - Python 3.11
  - used for focused pytest checks
  - PyTorch/DockQ may be absent here; use the dedicated environments above
- FreeSASA is an explicit optional backend; its interpreter is passed with
  `--freesasa_python` or `PRISM_FREESASA_PYTHON`.

## Stable working pipeline recipe

The canonical stable-run notes are in `docs/STABLE_PIPELINE.md`. In summary,
use `benchmark/scripts/run_prism_pipeline_smoke.sh` for an isolated local
smoke check, and use a Slurm launcher that first executes:

```bash
source /opt/ohpc/admin/lmod/lmod/init/bash
module load rosetta/2022.42
```

The validated batch resource recipe is `cosbi`, 2 CPUs, 4G, and 30 minutes,
with `/home/rshadi25/.conda/envs/gtalign_env/bin/python`. Keep the production
thresholds at their source defaults (`TM=0.5`, minimum matches `15`, match
percentage `50`, difference allowance `20`, clash distance `3`, maximum
clashes `5`, scaffold threshold `5.0`). The relaxed settings used in rescue
experiments are diagnostic-only and must be recorded in the candidate audit.

## Important Directories

- Source: `src/`
- Tests: `tests/`
- Current templates/interfaces: `templates/interfaces/` and
  `templates/interfaces_lists/`
- Current downloaded PDBs: `processed/pdbs/` (four-letter filenames)
- Current surface outputs: `processed/surface_extraction/`
- Current alignments: `processed/alignment/`
- Current transformations: `processed/transformation/`
- Current Rosetta outputs: `processed/rosetta_refinement/structures/`
- Current Rosetta scores: `processed/rosetta_refinement/energies/`
- Derived experiments: `tmp/agent/`
- Benchmark native complexes: `benchmark/prism_processed/results/native_bound_complexes_t_rigid/`
- ML design: `docs/superpowers/specs/2026-07-11-biological-tmalign-ranking-design.md`
- ML workflow: `docs/ML_TRAINING.md`
- Pilot report: `docs/biological-ranking-pilot-20260711.md`

## Implemented Learning Workflow

- Candidate audit: `src/candidate_audit.py`
- Deterministic biological baseline: `src/candidate_ranker.py`
- DockQ/native-like labels and grouped split validation:
  `src/ranking_data.py`
- Ranking metrics: `src/ranking_metrics.py`
- Optional residue-contact PyTorch scaffold: `src/residue_contact_model.py`
- Candidate table builder: `benchmark/scripts/build_candidate_table.py`
- DockQ/iRMSD label attachment: `benchmark/scripts/attach_native_labels.py`
- Baseline CSV ranking: `benchmark/scripts/rank_candidate_table.py`
- Optional Stage 1 trainer: `benchmark/scripts/train_reranker.py`
- Low-resource Slurm trainer: `benchmark/jobs/train_reranker.sbatch`
- Live transformation auditing is opt-in with
  `PRISM_CANDIDATE_AUDIT_PATH`.

Training rules:

- Use only rows with native labels; retained failures/unlabeled rows are not
  negative training examples.
- Require `native_complex_id` and at least two independent native complexes.
- Prefer sequence-cluster/grouped splits to prevent decoy leakage.
- Keep the learned model disabled by default until held-out DockQ/iRMSD results
  beat the deterministic baseline.

## Validated Evidence

### Real `1ahw` current-pipeline pilot

Model chains `A,H` were mapped to native chains `B,C`.

| Candidate | DockQ | iRMSD |
|---|---:|---:|
| `1h5bAB_o1` | 0.031 | 16.937 |
| `3lqmAB_o1` | 0.004 | 27.240 |
| `3lqmAB_o2` | 0.006 | 31.949 |
| `2f0xEH_o2` | 0.004 | 33.521 |

All four are below the native-like DockQ threshold `0.23`. The deterministic
baseline selected `1h5bAB_o1`, which was best by both metrics, but this is only
one complex and is not a general accuracy claim.

Pilot artifacts:

- `tmp/agent/20260711-biological-ranking/current-1ahw-candidates.csv`
- `tmp/agent/20260711-biological-ranking/current-1ahw-ranked.csv`
- `tmp/agent/20260711-biological-ranking/current-1ahw-labeled.csv`
- `tmp/agent/20260711-biological-ranking/current-1ahw-dockq-original.csv`
- `tmp/agent/20260711-biological-ranking/current-1ahw-irmsd.csv`

### Current `1fqjBE` diagnostic

Valid target IDs are `1TNDC` and `1FQIA`; underscore forms such as `1TND_C`
and `1FQI_A` are rejected by the current downloader.

Default run `1343749` completed with exit `0` but `Passed pairs 0` because
TM-align scores below `0.5` fail the default alignment gate.

Relaxed alignment run `1343806` completed with transformed outputs but no
accepted pair. Clash counts were `15` (`o1`) and `32` (`o2`), while the default
maximum is `5`. Diagnostic threshold overrides are available through:

```text
PRISM_TM_SCORE_THRESHOLD
PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE
PRISM_DIFF_PERCENTAGE
PRISM_CLASHING_DISTANCE
PRISM_MAX_CLASHING_COUNT
```

Defaults must not be changed for production runs.

## Slurm/Connector Notes

- Submit and monitor from `login01` whenever possible.
- Heavy work must run on compute nodes through Slurm.
- Codex connectors/browser authentication are more reliable outside compute
  nodes; a local Codex client plus login-node SSH orchestration is preferred.
- `ai01` intermittently failed to contact the Slurm controller with:
  `Unable to contact slurm controller (connect failure)`.
- Jobs observed:
  - `1343749`: AI default current test, completed successfully, zero pairs
  - `1343750`: KUTEM default current test, pending by priority at last check
  - `1343806`: AI relaxed-alignment test, completed, clash-rejected candidates
  - `1343808`: KUTEM relaxed-alignment test, pending by priority at last check
  - `1343726`: failed because obsolete smoke CLI arguments were passed
  - `1343728`: failed because `--generate_templates false` evaluated as true and
    target IDs contained underscores

## Next Agenda

1. From `login01`, check `squeue`/`sacct` and collect jobs `1343750` and
   `1343808` logs and artifacts.
2. Run the opt-in `PRISM_MAX_CLASHING_COUNT=100` diagnostic only as a rescue
   experiment; inspect whether Rosetta produces a model and score it with
   `/scratch/tmp/prism-dockq-env/bin/python`.
3. Stage at least one additional independent current-compatible native complex.
4. Build a combined candidate table with `native_complex_id`, DockQ, and iRMSD.
5. Train only after at least two native complexes exist; report grouped held-out
   metrics and compare against the deterministic baseline.
6. Keep `1WDW_BD:A` out of the current run until multi-chain query support is
   explicitly implemented and tested; `1V8Z_AB` is incompatible with the
   current single-chain target contract.
