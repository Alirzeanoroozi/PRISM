# PRISM Benchmark Analysis Workspace

This folder is the clean, reproducible entry point for benchmark analysis from PRISM outputs.

It intentionally exposes only the benchmark-relevant assets:

- `results/`
  - benchmark outputs, reports, validation outputs, native bound complexes, and summary CSVs
- `env/`
  - environment package specification for the benchmark stage

This workspace now contains only the benchmark environment, benchmark outputs, and documentation. The benchmark CSVs remain in `benchmark/data`, and the PRISM raw archive should be downloaded from Google Drive when reproducing the analysis.

## Required Inputs

- Benchmark CSVs from `benchmark/data`:
  - `benchmark/data/T_Rigid.csv`
  - `benchmark/data/T_medium.csv`
  - `benchmark/data/T_difficult.csv`
- PRISM raw archive from Google Drive:
  - `https://drive.google.com/file/d/1znDfQGJ1SpKnhAX-zzAUC5M8KVg71qgm/view?usp=sharing`
- Extracted PRISM raw content:
  - `benchmark/prism_processed/prism_raw`
- Chain-fixed PRISM models produced during the run:
  - `benchmark/prism_processed/prism_raw_chainfixed`

## Required Scripts

These scripts live in `benchmark/scripts`, not inside `benchmark/prism_processed`:

- `benchmark/scripts/fix_model_chain_names.py`
- `benchmark/scripts/analyze_prism_rigid_results.py`
- `benchmark/scripts/analyze_prism_all_benchmarks.py`
- `benchmark/scripts/dockq.py`
- `benchmark/scripts/irmsd.py`
- `benchmark/scripts/generate_prism_benchmark_report.py`
- `benchmark/scripts/validate_prism_pipeline.py`
- `benchmark/scripts/submit_prism_analysis_all.sbatch`

## Environment

The benchmark stage needs a Python environment with:

- `numpy`
- `biopython`
- `matplotlib`
- `reportlab`
- `dockq`

Package versions used for the current benchmark run are frozen in:

- `env/requirements_prism_benchmark.txt`
- `env/prism_score_env` is the copied benchmark environment used for the existing results

If creating a fresh virtual environment:

```bash
python3 -m venv benchmark/prism_processed/env/prism_score_env_fresh
source benchmark/prism_processed/env/prism_score_env_fresh/bin/activate
pip install --upgrade pip
pip install -r benchmark/prism_processed/env/requirements_prism_benchmark.txt
```

Check the critical binary:

```bash
benchmark/prism_processed/env/prism_score_env/bin/DockQ --help
```

If you create a fresh environment, replace the `--score-python` path in the commands below with:

```bash
benchmark/prism_processed/env/prism_score_env_fresh/bin/python
```

## Reproducible Flow

Workflow summary:

1. Download and extract the PRISM raw archive.
2. Fix wrong or repeated model chain IDs and record the corrected chain mapping.
3. Match PRISM models to benchmark rows using only `PDB ID 1` and `PDB ID 2`.
4. Score matched models against native bound complexes with `DockQ` and `iRMSD`.
5. Aggregate per-model results into pair-level summaries for rigid, medium, and difficult sets.
6. Generate plots and the final benchmark report.
7. Validate the produced outputs for consistency and sampled score reproducibility.

### One-command flow

```bash
mkdir -p benchmark/prism_processed/prism_raw
unzip /path/to/prism_raw.zip -d benchmark/prism_processed/prism_raw

python3 -u benchmark/scripts/fix_model_chain_names.py \
  --rosetta-root benchmark/prism_processed/prism_raw/rosetta_output_1 \
  --fixed-root benchmark/prism_processed/prism_raw_chainfixed \
  --map-csv benchmark/prism_processed/results/chain_fix_map.csv

python3 -u benchmark/scripts/analyze_prism_all_benchmarks.py \
  --rosetta-root benchmark/prism_processed/prism_raw/rosetta_output_1 \
  --data-dir benchmark/data \
  --results-root benchmark/prism_processed/results \
  --score-python benchmark/prism_processed/env/prism_score_env/bin/python \
  --metrics dockq,irmsd \
  --sets rigid,medium,difficult \
  --model-fix-map-csv benchmark/prism_processed/results/chain_fix_map.csv \
  --score-timeout-sec 5 \
  --verbose-every 50 \
  --download-native-missing

benchmark/prism_processed/env/prism_score_env/bin/python -u \
  benchmark/scripts/generate_prism_benchmark_report.py \
  --results-root benchmark/prism_processed/results \
  --out-dir benchmark/prism_processed/results/benchmark_reports \
  --data-dir benchmark/data

python3 -u benchmark/scripts/validate_prism_pipeline.py \
  --results-root benchmark/prism_processed/results \
  --score-python benchmark/prism_processed/env/prism_score_env/bin/python \
  --sample-size-per-set 5 \
  --timeout-sec 10 \
  --out-dir benchmark/prism_processed/results/validation
```

If your starting point is the Google Drive raw archive:

```bash
wget --no-check-certificate -O benchmark/prism_processed/prism_raw.zip \
  "https://drive.google.com/uc?export=download&id=1znDfQGJ1SpKnhAX-zzAUC5M8KVg71qgm"

mkdir -p benchmark/prism_processed/prism_raw
unzip benchmark/prism_processed/prism_raw.zip -d benchmark/prism_processed/prism_raw
```

Supported archive types:

- `.zip`
- `.tar.gz`
- `.tgz`
- `.tar`

The archive should be extracted into `benchmark/prism_processed/prism_raw`. If the extracted content contains a `rosetta_output_1/` folder, use that folder as `--rosetta-root`. Otherwise use `benchmark/prism_processed/prism_raw` itself.

### Slurm submission

If you want to run the full benchmark on HPC, use:

```bash
sbatch benchmark/scripts/submit_prism_analysis_all.sbatch
```

Optional overrides:

```bash
sbatch --export=ALL,SETS=rigid,medium,difficult,SCORE_TIMEOUT_SEC=5 \
  benchmark/scripts/submit_prism_analysis_all.sbatch
```

The batch script uses the same paths documented in this README:

- input models: `benchmark/prism_processed/prism_raw/rosetta_output_1`
- chain-fixed models: `benchmark/prism_processed/prism_raw_chainfixed`
- benchmark CSVs: `benchmark/data`
- results: `benchmark/prism_processed/results`
- scoring environment: `benchmark/prism_processed/env/prism_score_env`

### Step-by-step

1. Download the PRISM raw archive from the Google Drive link above and extract it into `benchmark/prism_processed/prism_raw`.

2. Fix chain naming issues in PRISM model files and write the correction map:

```bash
python3 -u benchmark/scripts/fix_model_chain_names.py \
  --rosetta-root benchmark/prism_processed/prism_raw/rosetta_output_1 \
  --fixed-root benchmark/prism_processed/prism_raw_chainfixed \
  --map-csv benchmark/prism_processed/results/chain_fix_map.csv
```

3. Run benchmark analysis for all benchmark sets:

```bash
python3 -u benchmark/scripts/analyze_prism_all_benchmarks.py \
  --rosetta-root benchmark/prism_processed/prism_raw/rosetta_output_1 \
  --data-dir benchmark/data \
  --results-root benchmark/prism_processed/results \
  --score-python benchmark/prism_processed/env/prism_score_env/bin/python \
  --metrics dockq,irmsd \
  --sets rigid,medium,difficult \
  --model-fix-map-csv benchmark/prism_processed/results/chain_fix_map.csv \
  --score-timeout-sec 5 \
  --verbose-every 50 \
  --download-native-missing
```

4. Generate the report and plots:

```bash
benchmark/prism_processed/env/prism_score_env/bin/python -u \
  benchmark/scripts/generate_prism_benchmark_report.py \
  --results-root benchmark/prism_processed/results \
  --out-dir benchmark/prism_processed/results/benchmark_reports \
  --data-dir benchmark/data
```

5. Validate the produced outputs:

```bash
python3 -u benchmark/scripts/validate_prism_pipeline.py \
  --results-root benchmark/prism_processed/results \
  --score-python benchmark/prism_processed/env/prism_score_env/bin/python \
  --sample-size-per-set 5 \
  --timeout-sec 10 \
  --out-dir benchmark/prism_processed/results/validation
```

## Main Outputs

Per-set benchmark outputs:

- `results/prism_rigid_analysis_all_jobs/per_prediction.csv`
  - one row per analyzed model-complex mapping in the rigid set
- `results/prism_rigid_analysis_all_jobs/pair_summary.csv`
  - aggregated rigid-set summary per `(PDB ID 1, PDB ID 2, complex)` group
- `results/prism_medium_analysis_all_jobs/per_prediction.csv`
  - one row per analyzed model-complex mapping in the medium set
- `results/prism_medium_analysis_all_jobs/pair_summary.csv`
  - aggregated medium-set summary per `(PDB ID 1, PDB ID 2, complex)` group
- `results/prism_difficult_analysis_all_jobs/per_prediction.csv`
  - one row per analyzed model-complex mapping in the difficult set
- `results/prism_difficult_analysis_all_jobs/pair_summary.csv`
  - aggregated difficult-set summary per `(PDB ID 1, PDB ID 2, complex)` group

Shared outputs:

- `results/chain_fix_map.csv`
  - map of original model files to corrected chain assignments and fixed model paths
- `results/benchmark_reports/cross_benchmark_report.md`
  - human-readable benchmark report in Markdown
- `results/benchmark_reports/cross_benchmark_report.pdf`
  - PDF version of the benchmark report
- `results/benchmark_reports/error_files.csv`
  - list of files that produced DockQ or iRMSD errors
- `results/predicted_models_summary.csv`
  - summary of how many unique predicted models belong to each benchmark set
- `results/repeated_pairs_different_templates.csv`
  - summary of PDB ID pairs that have multiple predicted models from different templates
- `results/validation/pipeline_validation.md`
  - readable validation report for the generated outputs
- `results/validation/pipeline_validation.json`
  - machine-readable validation results

## Notes

- Matching is based on `PDB ID 1` and `PDB ID 2` from model filenames, not the first filename token.
- Chain fixing is required because some PRISM model PDBs encode multiple segments under a wrong or repeated chain ID.
- Native bound complexes do not need to be uploaded in advance if you run analysis with `--download-native-missing`.
- `benchmark/prism_processed` is intentionally results/environment/documentation focused. The executable benchmark scripts remain in `benchmark/scripts`, and benchmark CSVs remain in `benchmark/data`.
