# Stable Current PRISM Pipeline

This is the reproducible current TM-align + Rosetta recipe for
`PRISM-prescript`. Keep diagnostic threshold overrides separate from this
default path.

## Local smoke check

Run from the repository root:

```bash
PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  bash benchmark/scripts/run_prism_pipeline_smoke.sh
```

The smoke runner stages an isolated workspace under `tmp/agent/`, uses the
default `prism.py` behavior (`--generate_templates` is omitted), and records
`run.log`, alignment, transformation, and Rosetta directories. A successful
zero-pair smoke run is still a valid pipeline health check; it does not imply
that the benchmark pair is biologically positive.

## Slurm Rosetta run

Rosetta must be loaded explicitly in every batch launcher:

```bash
source /opt/ohpc/admin/lmod/lmod/init/bash
module load rosetta/2022.42
python_bin=/home/rshadi25/.conda/envs/gtalign_env/bin/python
cd <isolated-workspace>
"$python_bin" prism.py 2>&1 | tee run.log
```

The validated resource recipe is `cosbi`, `2` CPUs, `4G`, and `00:30:00`;
submit and monitor from `login01` where possible. Record the Slurm job ID,
partition, node, launcher, input pair, template list, and work directory.

The stable production defaults are retained in source:

- TMalign TM-score threshold: `0.5`
- minimum matched residues: `15`
- minimum match percentage: `50.0`
- percentage difference allowance: `20.0`
- clash distance: `3.0`
- maximum clashes: `5`
- surface scaffold compatibility threshold: `5.0`

Use five-character single-chain IDs such as `1TNDC` and `1FQIA`; transformed
inputs are materialized as four-letter PDB files under `processed/pdbs/`.

## Optional deterministic ranking

Ranking is disabled by default and does not change the stable pipeline. To
retain only the deterministic top candidates per input pair before refinement,
add `--rank --top-k <positive-integer>` to a normal command. A ranked run
automatically writes a fresh audit file under `processed/candidate_audit/`; use
`--candidate-audit-path <path>` only when an explicit append-only audit
destination is required. The current baseline ranks TM-score, coverage, and a
known-clash penalty. It is a resource-reduction experiment, not a validated
learned-ranker replacement.

## Validation and scoring

Focused current-pipeline checks:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_transformation_audit.py \
  tests/test_transformation_clashes.py \
  tests/test_transformation_thresholds.py \
  tests/test_train_reranker.py \
  tests/test_ranking_data.py \
  tests/test_candidate_ranker.py
```

Score completed structures with:

```bash
/scratch/tmp/prism-dockq-env/bin/python \
  benchmark/scripts/score_single_prism_pair.py <structures-dir> <native.pdb> \
  --model-receptor A --model-ligand B \
  --native-receptor <native-chain-1> --native-ligand <native-chain-2> \
  --score-python /scratch/tmp/prism-dockq-env/bin/python \
  --out-csv <scores.csv>
```

## Diagnostic rescue only

The following overrides were used to investigate weak alignments and are not
stable production settings:

```bash
export PRISM_MINIMUM_RESIDUE_MATCH_COUNT=8
export PRISM_TM_SCORE_THRESHOLD=0.35
export PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE=35.0
export PRISM_DIFF_PERCENTAGE=35.0
export PRISM_MAX_CLASHING_COUNT=100
export PRISM_CANDIDATE_AUDIT_PATH="$work/candidate_audit.jsonl"
```

The reproducible rescue launcher is retained at
`tmp/agent/20260712-legacy-template-diagnostic/run_1bpb_3k77_legacytemplate_cosbi.sbatch`.
It loaded `rosetta/2022.42` and produced DockQ `0.770` / iRMSD `0.922` for
`1BPB/3K77`; that result validates the diagnostic path, not a change to the
stable thresholds.
