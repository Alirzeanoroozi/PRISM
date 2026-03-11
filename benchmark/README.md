# Benchmark Analysis

Run the benchmark pipeline with the Slurm submission script:

```bash
sbatch benchmark/scripts/submit_prism_analysis_all.sbatch
```

This job downloads the PRISM raw archive, extracts it into `benchmark/prism_processed/prism_raw`, fixes problematic model chain names, runs benchmark scoring for the rigid, medium, and difficult sets, downloads missing native bound complexes when needed, and writes the processed benchmark outputs under `benchmark/prism_processed/results`.

## Pipeline

1. Download PRISM raw archive from Google Drive.
2. Extract the archive and detect `rosetta_output_1`.
3. Rewrite model PDB chain IDs when PRISM outputs have incorrect repeated chain naming.
4. Match prediction files to benchmark rows using only `PDB ID 1` and `PDB ID 2`.
5. Download native bound complexes if they are missing locally.
6. Run `DockQ` and `iRMSD` scoring for matched model/complex comparisons.
7. Aggregate per-prediction scores into pair-level summaries.

## Main Outputs

Per benchmark set:

- `benchmark/prism_processed/results/prism_rigid_analysis_all_jobs/per_prediction.csv`
  - one row per analyzed model-to-complex mapping in the rigid set
- `benchmark/prism_processed/results/prism_rigid_analysis_all_jobs/pair_summary.csv`
  - aggregated rigid-set summary per matched `(PDB ID 1, PDB ID 2, complex)`
- `benchmark/prism_processed/results/prism_medium_analysis_all_jobs/per_prediction.csv`
  - one row per analyzed model-to-complex mapping in the medium set
- `benchmark/prism_processed/results/prism_medium_analysis_all_jobs/pair_summary.csv`
  - aggregated medium-set summary per matched `(PDB ID 1, PDB ID 2, complex)`
- `benchmark/prism_processed/results/prism_difficult_analysis_all_jobs/per_prediction.csv`
  - one row per analyzed model-to-complex mapping in the difficult set
- `benchmark/prism_processed/results/prism_difficult_analysis_all_jobs/pair_summary.csv`
  - aggregated difficult-set summary per matched `(PDB ID 1, PDB ID 2, complex)`

Shared outputs:

- `benchmark/prism_processed/results/chain_fix_map.csv`
  - map from original model files to corrected chain assignments and fixed-model paths
- `benchmark/prism_processed/results/native_bound_complexes_t_rigid/`
- `benchmark/prism_processed/results/native_bound_complexes_t_medium/`
- `benchmark/prism_processed/results/native_bound_complexes_t_difficult/`
  - downloaded native bound complexes used during scoring
- `benchmark/prism_processed/results/predicted_models_summary.csv`
  - summary of how many unique prediction files belong to each benchmark set
- `benchmark/prism_processed/results/repeated_pairs_different_templates.csv`
  - repeated `(PDB ID 1, PDB ID 2)` groups that have multiple templates
- `benchmark/prism_processed/results/benchmark_reports/cross_benchmark_report.md`
  - markdown benchmark report
- `benchmark/prism_processed/results/benchmark_reports/cross_benchmark_report.pdf`
  - PDF benchmark report
- `benchmark/prism_processed/results/benchmark_reports/error_files.csv`
  - model/complex rows that failed scoring, with reasons
- `benchmark/prism_processed/results/validation/pipeline_validation.md`
- `benchmark/prism_processed/results/validation/pipeline_validation.json`
  - validation outputs for aggregation and sampled score recomputation

## Notes

- Benchmark CSV inputs remain in `benchmark/data`.
- Benchmark scripts remain in `benchmark/scripts`.
- The detailed reproducibility guide is in `benchmark/prism_processed/README.md`.
