# DiffMaSIF Replacement Experiment

This folder is an isolated experiment for testing whether a `DiffMaSIF`-style
surface representation can replace the residue-level relative accessibility that
PRISM currently gets from `Naccess`.

## Why this is only an experiment

PRISM consumes residue-level relative ASA in two places:

- `src/hotspot.py` for HotPoint-style buried residue filtering
- `src/surface_extract.py` for target surface residue selection

`Naccess` produces residue-level solvent accessibility directly. `DiffMaSIF`
instead works on a molecular surface point cloud with learned geometric and
chemical features. That means it is not a drop-in replacement for `.rsa`
output. The experiment here maps surface points back to residues and produces a
PRISM-compatible proxy score.

## Files

- `export_residue_accessibility.py`
  Converts a CSV surface point cloud into residue-level accessibility proxy
  scores by assigning each point to the nearest residue atom.
- `compare_with_naccess.py`
  Runs the existing PRISM `Naccess` path and compares it against the point-cloud
  proxy on the same target or template structure.
- `generate_geometric_surface_points.py`
  Generates a coarse geometry-only solvent-surface point cloud from a PDB file.
  This is not DiffMaSIF, but it lets you test whether surface-point-to-residue
  mapping is strong enough to approximate PRISM's current `Naccess` signal.
- `evaluate_replacement_candidate.py`
  Batch harness that generates geometric surface points, compares them against
  `Naccess`, and writes a summary CSV.

## Expected surface point input

The adapter expects a CSV file with at least these columns:

- `x`
- `y`
- `z`

Optional columns:

- `score`
- `feature`

If `score` is present, the comparison output also reports score-weighted
residue statistics.

## Example workflow

1. Generate a surface point cloud externally from a MaSIF / DiffMaSIF-style
   preprocessing pipeline.
2. Save the point cloud as CSV.
3. Compare it against the current PRISM `Naccess` output:

```bash
python3 tests/diffmasif_replacement/compare_with_naccess.py \
  --structure-kind target \
  --structure-id 1fgnA \
  --surface-points /path/to/diffmasif_surface_points.csv \
  --assignment-cutoff 4.0 \
  --output tests/diffmasif_replacement/output/1fgnA_comparison.json
```

Template example:

```bash
python3 tests/diffmasif_replacement/compare_with_naccess.py \
  --structure-kind template \
  --structure-id 1a28AB \
  --surface-points /path/to/diffmasif_surface_points.csv \
  --assignment-cutoff 4.0 \
  --output tests/diffmasif_replacement/output/1a28AB_comparison.json
```

Geometry-only local smoke test:

```bash
python3 tests/diffmasif_replacement/evaluate_replacement_candidate.py \
  --structure-kind target \
  --structure-id 1fgnA
```

Real MaSIF surface test from a generated `.ply` mesh:

```bash
python3 tests/diffmasif_runtime/export_ply_vertices_to_csv.py \
  --ply /scratch/rshadi25/GitHub/masif/data/masif_site/data_preparation/01-benchmark_surfaces/1FGN_H.ply \
  --output tests/diffmasif_runtime/output/1FGN_H.surface_points.csv

python3 tests/diffmasif_replacement/compare_with_naccess.py \
  --structure-kind target \
  --structure-id 1FGNH \
  --surface-points tests/diffmasif_runtime/output/1FGN_H.surface_points.csv \
  --assignment-cutoff 4.0 \
  --output tests/diffmasif_replacement/output/1FGN_H.masif_surface.comparison.json
```

## Output meaning

The generated JSON includes:

- `naccess_relative_asa`
- `diffmasif_proxy_relative_asa`
- `per_residue`
- `summary`

The proxy score is not a physical ASA value. It is a normalized residue surface
coverage derived from assigned surface points:

- more assigned surface points -> more exposed residue
- fewer assigned surface points -> more buried residue

This is enough to test whether a surface-model signal can substitute for the
`Naccess` thresholds used by PRISM before doing a deeper pipeline rewrite.

## How to interpret the local smoke test

The geometry-only generator is weaker than DiffMaSIF because it has no learned
features. If even this coarse point-cloud representation cannot recover the
`Naccess` burial and surface thresholds with decent agreement, a DiffMaSIF-based
replacement is unlikely to be plug-compatible without re-tuning PRISM. If the
agreement is decent, the next step is to feed real DiffMaSIF surface points into
the same comparison script and see whether they do better.

Use `--assignment-cutoff 4.0` as the starting point for this repo. Smaller
values under-assign surface points to residues for the generated shell points.

## Current results in this repo

For target `1FGNH`, the geometry-only proxy remains the stronger local result:

- geometry-only proxy: buried agreement `0.8318`, surface agreement `0.8738`, MAE `8.94`
- real MaSIF surface vertices: buried agreement `0.8131`, surface agreement `0.8224`, MAE `11.83`

That means a MaSIF-generated surface can reproduce PRISM's `Naccess`-driven
residue decisions to a useful degree, but this single run is not strong enough
to claim a drop-in replacement yet. The next step is to test more targets and,
if needed, retune the residue thresholds away from the current `Naccess`-based
`15` and `20` cutoffs.
