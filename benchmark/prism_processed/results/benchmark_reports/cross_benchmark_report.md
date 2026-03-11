# PRISM Benchmark Analysis Report (Rigid / Medium / Difficult)

## Scope

This report summarizes PRISM prediction analysis against `T_Rigid.csv`, `T_medium.csv`, and `T_difficult.csv` using `DockQ` and `iRMSD`.

## Method Logic (Why Scores Are Not Read From PRISM Output Files)

- The target metrics in this analysis are **benchmark comparison metrics** (`DockQ`, `iRMSD`) that require a comparison between:
- a PRISM predicted complex model (`.pdb`) and
- the benchmark native bound complex PDB with a specific chain mapping from the benchmark CSV.
- PRISM output files (e.g. Rosetta score files / intermediate files) do not directly provide these benchmark-comparison metrics with the required benchmark-native chain mapping.
- Therefore, the pipeline recomputes `DockQ` and `iRMSD` from the model PDB + downloaded native bound PDB + CSV-defined chain roles.
- The chain fields reported in `per_prediction.csv` are taken from the model filename (or inferred fallback if omitted), as requested.

## Output Locations

- `benchmark/prism_processed_results/prism_rigid_analysis_all_jobs`
- `benchmark/prism_processed_results/prism_medium_analysis_all_jobs`
- `benchmark/prism_processed_results/prism_difficult_analysis_all_jobs`
- Plots: `benchmark/prism_processed_results/benchmark_reports/plots`

## Summary Table

| Set | Complexes in CSV | Total Models | Matched Models | Unmatched Models | DockQ Scored Rows | iRMSD Scored Rows | DockQ Errors | iRMSD Errors |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Rigid | 162 | 572 | 371 | 201 | 378 | 374 | 21 | 25 |
| Medium | 60 | 572 | 134 | 438 | 120 | 118 | 14 | 16 |
| Difficult | 35 | 572 | 67 | 505 | 67 | 67 | 0 | 0 |

Abbreviation guide:
- `DockQ` = DockQ score.
- `iRMSD` = interface RMSD.
- `Scored Rows` = matched model-complex rows where numeric score was produced.

## Predicted Model Membership by Set

| Set | Models Belonging to Set | Percent of All Predicted Models |
|---|---:|---:|
| Rigid | 371 | 64.86% |
| Medium | 134 | 23.43% |
| Difficult | 67 | 11.71% |

Abbreviation guide:
- `Set` = benchmark subset (`Rigid`, `Medium`, `Difficult`).
- `Models Belonging to Set` = unique model files whose PDB ID pair matches at least one row in that set.

## Complex-Level Coverage

| Set | Complexes With Found Models | Complexes With >=1 DockQ | Complexes With >=1 iRMSD | Matched Models Missing DockQ | Matched Models Missing iRMSD |
|---|---:|---:|---:|---:|---:|
| Rigid | 110 | 109 | 107 | 21 | 25 |
| Medium | 34 | 34 | 33 | 14 | 16 |
| Difficult | 17 | 17 | 17 | 0 | 0 |

Abbreviation guide:
- `Cx` = complexes.
- `DQ` = DockQ.
- `IR` = iRMSD.
- `Match no DQ/IR` = matched rows lacking numeric DockQ/iRMSD output.

## Code Logic

1. Match model files to benchmark rows using only `(PDB ID 1, PDB ID 2)` partner names from model filename, ignoring the first template/target token.
2. Assign model partner chains to `PDB ID 1` and `PDB ID 2` roles, including swapped partner-name order handling.
3. Expand each complex row to all native receptor-ligand chain pairs (e.g., `1AHW_AB:C` -> `AC` and `BC`).
4. For each native pair, run DockQ and iRMSD against the downloaded native bound complex PDB.
5. Run reciprocal iRMSD order as well, keep per-model best (`min`) across directions and native pairs.
6. Aggregate by `(PDB ID 1, PDB ID 2, complex)` group: mean/variance/best across predictions/templates.

## Key Insights

- Rigid: 371/572 unique models matched the benchmark set (64.9%); DockQ scored 94.7% of matched rows and iRMSD scored 93.7%.
- Medium: 134/572 unique models matched the benchmark set (23.4%); DockQ scored 89.6% of matched rows and iRMSD scored 88.1%.
- Difficult: 67/572 unique models matched the benchmark set (11.7%); DockQ scored 100.0% of matched rows and iRMSD scored 100.0%.
- Best DockQ across all analyzed rows was 0.859 in Rigid (3e1zAB_2nnrA_0_3bpfA_0.rosetta_0001_0001.pdb; pair 2nnr|3bpf, model chains B/A, native pair A:B).
- Best (lowest) iRMSD across all analyzed rows was 0.000 A in Rigid (1xu1BS_1u5yB_0_1xutA_0.rosetta_0001_0001.pdb; pair 1u5y|1xut).

## Test Case (1AHW)

- Benchmark row: `Complex=1AHW_AB:C`, `PDB ID 1=1FGN_LH`, `PDB ID 2=1TFH_A`.
- Model A: `1ahwBC_1fgnH_0_1tfhA_0.rosetta_0001_0001.pdb`
- Model chains used: receptor=`H`, ligand=`A`.
- Native pairs evaluated: `AC,BC`.
- DockQ per pair: `AC=0.036409;BC=0.514546`.
- iRMSD forward per pair: `AC=11.022;BC=1.706`.
- iRMSD reverse per pair: `AC=11.022;BC=1.706`.
- Best selected: DockQ=`0.5145464053742148` (Medium), iRMSD=`1.706`.

Core logic used in this case:

```python
# complex row: 1AHW_AB:C -> native pairs: AC, BC
native_pairs = [('A','C'), ('B','C')]
# model 1ahwBC_1fgnH_0_1tfhA_0 -> model chains for PDB ID1/2 roles
model_r, model_l = 'H', 'A'
dockq_vals = []
irmsd_vals = []
for ref_r, ref_l in native_pairs:
    dockq = DockQ(model, native_1ahw, mapping=f'{model_r}{model_l}:{ref_r}{ref_l}')
    ir_fwd = iRMSD(model, model_r, model_l, native_1ahw, ref_r, ref_l)
    ir_rev = iRMSD(model, model_l, model_r, native_1ahw, ref_l, ref_r)
    dockq_vals.append(dockq)
    irmsd_vals.extend([ir_fwd, ir_rev])
best_dockq = max(dockq_vals)
best_irmsd = min(irmsd_vals)
```
- Model B (alt partner naming): `1ahwBC_1jptH_0_1tfhB_0.rosetta_0001_0001.pdb` -> model chains `H/B`, native pairs `HT,LT`.

## Plots

Generated plot files:
- `benchmark/prism_processed_results/benchmark_reports/plots/rigid_status_bar.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/rigid_dockq_histogram.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/rigid_irmsd_histogram.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/rigid_dockq_vs_irmsd.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/medium_status_bar.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/medium_dockq_histogram.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/medium_irmsd_histogram.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/medium_dockq_vs_irmsd.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/difficult_status_bar.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/difficult_dockq_histogram.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/difficult_irmsd_histogram.png`
- `benchmark/prism_processed_results/benchmark_reports/plots/difficult_dockq_vs_irmsd.png`

Plot confirmation notes:
- Per-set status bars use `Matched`, `DockQ Scored`, `iRMSD Scored`, `DockQ Errors`, `iRMSD Errors` from each set summary.
- Histograms and scatter use only rows with numeric metric values in `per_prediction.csv`.

## Rigid Set Details

- Pair-summary rows: `110`
- DockQ scored rows: `378` (mean = `0.214`, max = `0.859`)
- iRMSD scored rows: `374` (mean = `11.358`, min = `0.000`)
- DockQ CAPRI classes: `Acceptable=43, High=11, Incorrect=250, Medium=74`
- DockQ failure reasons: `No interface found in native for selected mapping=1, Other DockQ error=20`
- iRMSD failure reasons: `Other iRMSD error=20, iRMSD script runtime failure=5`

Top 5 pairs by best DockQ:
- `2nnr|3bpf` (`2OUL_A:B`): best DockQ = `0.859`, n = `1`
- `1a19|1rgh` (`1AY7_A:B`): best DockQ = `0.854`, n = `1`
- `3pc6|3pc7` (`3PC8_A:C`): best DockQ = `0.843`, n = `6`
- `1ccp|1ycc` (`2PCC_A:B`): best DockQ = `0.840`, n = `3`
- `1c3d|2gom` (`3D5S_A:C`): best DockQ = `0.834`, n = `2`

Top 5 pairs by best (lowest) iRMSD:
- `1u5y|1xut` (`1XU1_ABD:T`): best iRMSD = `0.000`, n = `2`
- `3pc6|3pc7` (`3PC8_A:C`): best iRMSD = `0.600`, n = `6`
- `1a19|1rgh` (`1AY7_A:B`): best iRMSD = `0.612`, n = `1`
- `1c3d|2gom` (`3D5S_A:C`): best iRMSD = `0.642`, n = `2`
- `3f5v|3rvt` (`3RVW_CD:A`): best iRMSD = `0.642`, n = `1`

## Medium Set Details

- Pair-summary rows: `34`
- DockQ scored rows: `120` (mean = `0.101`, max = `0.713`)
- iRMSD scored rows: `118` (mean = `13.507`, min = `0.000`)
- DockQ CAPRI classes: `Acceptable=16, Incorrect=96, Medium=8`
- DockQ failure reasons: `Other DockQ error=14`
- iRMSD failure reasons: `Other iRMSD error=14, iRMSD script runtime failure=2`

Top 5 pairs by best DockQ:
- `2g75|2ghv` (`2DD8_HL:S`): best DockQ = `0.713`, n = `2`
- `1iam|1mq9` (`1MQ8_A:B`): best DockQ = `0.587`, n = `1`
- `1l2z|1qgv` (`1SYX_A:B`): best DockQ = `0.584`, n = `1`
- `1hpt|2cga` (`1CGI_E:I`): best DockQ = `0.539`, n = `1`
- `1n0v|1xk9` (`1ZM4_A:B`): best DockQ = `0.531`, n = `2`

Top 5 pairs by best (lowest) iRMSD:
- `1h0c|2c0m` (`3R9A_AC:B`): best iRMSD = `0.000`, n = `2`
- `2g75|2ghv` (`2DD8_HL:S`): best iRMSD = `0.573`, n = `2`
- `1auq|1fvu` (`1IJK_A:BC`): best iRMSD = `1.712`, n = `2`
- `1mjn|3hi5` (`3HI6_XY:B`): best iRMSD = `1.867`, n = `1`
- `1jb1|2hpr` (`1KKL_ABC:H`): best iRMSD = `1.933`, n = `1`

## Difficult Set Details

- Pair-summary rows: `17`
- DockQ scored rows: `67` (mean = `0.193`, max = `0.701`)
- iRMSD scored rows: `67` (mean = `13.004`, min = `1.088`)
- DockQ CAPRI classes: `Acceptable=22, Incorrect=41, Medium=4`

Top 5 pairs by best DockQ:
- `1jpe|1jzo` (`1JZD_AB:C`): best DockQ = `0.701`, n = `8`
- `1a6z|1cx8` (`1DE4_AB:CF`): best DockQ = `0.555`, n = `2`
- `1qfk|1tfh` (`1FAK_HL:T`): best DockQ = `0.476`, n = `3`
- `3h13|4jj7` (`3H11_BC:A`): best DockQ = `0.466`, n = `20`
- `1qup|2jcw` (`1JK9_B:A`): best DockQ = `0.458`, n = `2`

Top 5 pairs by best (lowest) iRMSD:
- `1jpe|1jzo` (`1JZD_AB:C`): best iRMSD = `1.088`, n = `8`
- `1fnl|3ave` (`1E4K_AB:C`): best iRMSD = `1.428`, n = `4`
- `1a6z|1cx8` (`1DE4_AB:CF`): best iRMSD = `2.513`, n = `2`
- `1qfk|1tfh` (`1FAK_HL:T`): best iRMSD = `2.569`, n = `3`
- `1buy|1ern` (`1EER_A:BC`): best iRMSD = `2.841`, n = `3`

## Failure-Reason Notes

- `DockQ` failures are mostly due to: `No interface found in native for selected mapping`.
- This means the requested chain mapping exists syntactically, but DockQ cannot detect an interface in the native structure for that chain pair under its criteria.
- `iRMSD` failures are fewer and usually stem from chain/alignment/runtime issues for specific rows.
- Full error lists are exported to: `benchmark/prism_processed_results/benchmark_reports/error_files.csv`.
- Example error: set=`rigid`, metric=`dockq`, file=`2wp3OT_1tfhB_0_1jptL_0.rosetta_0001_0001.pdb`, complex=`1JPS_HL:T`, reason=`model_chain_missing_in_pdb:B`.

