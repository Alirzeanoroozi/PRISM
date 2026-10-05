# DiffMaSIF Runtime Setup

This folder is for getting a real `MaSIF` / `DiffMaSIF` runtime working next to
PRISM, without wiring it into the main pipeline yet.

## Goal

Produce a real surface-point output from a local `MaSIF`-family checkout, then
feed that output into:

- `tests/diffmasif_replacement/compare_with_naccess.py`

That lets you evaluate `DiffMaSIF` as an alternative to `Naccess` using the
same residue-level comparison harness already added to this repo.

## What you need

The official `LPDI-EPFL/masif` setup expects at least:

- `Python 3.6`
- `reduce`
- `MSMS 2.6.1`
- `Biopython`
- `PyMesh`
- `PDB2PQR`
- `multivalue`
- `APBS`
- `open3D`
- `TensorFlow 1.9`

Common environment variables expected by the upstream code:

- `APBS_BIN`
- `MULTIVALUE_BIN`
- `PDB2PQR_BIN`
- `REDUCE_HET_DICT`
- `PYMESH_PATH`
- `MSMS_BIN`

Note:
- The current cloned `MaSIF` code uses `source/triangulation/xyzrn.py` during
  `computeMSMS`, so `PDB2XYZRN` is not required for the checked-in Python path.

## Suggested local layout

Example:

```text
/scratch/rshadi25/GitHub/
  PRISM/
  masif/
```

Then run the wrapper from `PRISM` and point it at the external checkout.

## What this folder provides

- `check_runtime.py`
  Validates the runtime before you spend time debugging the upstream pipeline.
- `run_diffmasif_test.sh`
  Thin wrapper that checks the environment, points at a PRISM structure, and
  prints the next command you should run in the external `masif` checkout.
- `bootstrap_conda_env.sh`
  Creates a starting conda environment for the legacy Python stack.
- `masif_env_template.sh`
  Shell template for the required external-binary environment variables.
- `masif_env_masif_py36.sh`
  Workspace-specific env file with the paths that are already known after local
  installation into `masif_py36`.
- `export_ply_vertices_to_csv.py`
  Converts an upstream MaSIF `.ply` surface into a CSV with `x,y,z` columns that
  PRISM's replacement harness can consume.
- `batch_masif_compare.py`
  Runs MaSIF preprocessing for multiple chain-qualified target ids, converts the
  generated `.ply` surfaces to CSV, compares them against `Naccess`, and writes
  an aggregate summary.

## Basic validation

```bash
python3 tests/diffmasif_runtime/check_runtime.py \
  --masif-root /scratch/rshadi25/GitHub/masif \
  --pdb processed/pdbs/1fgn.pdb
```

## Wrapper example

```bash
bash tests/diffmasif_runtime/run_diffmasif_test.sh \
  /scratch/rshadi25/GitHub/masif \
  processed/pdbs/1fgn.pdb \
  L
```

## Conda bootstrap

This is only the Python side. It does not install `MSMS`, `APBS`, `PDB2PQR`,
`multivalue`, or `reduce`.

```bash
bash tests/diffmasif_runtime/bootstrap_conda_env.sh masif_py36
```

Then inspect and fill:

```bash
tests/diffmasif_runtime/masif_env_template.sh
```

## What the wrapper does

1. Confirms the external checkout path exists.
2. Confirms the PRISM PDB exists.
3. Checks upstream runtime dependencies and env vars.
4. Writes a staging directory under:

```text
tests/diffmasif_runtime/output/
```

5. Prints the next command to run in the external `masif` checkout.

## Why it stops short of full execution

The upstream `MaSIF` repository has multiple workflows and legacy dependency
constraints. This repo does not yet include a pinned local clone or a guaranteed
working environment for those tools. The wrapper therefore focuses on
repeatable validation and staging, which is the reliable next step.

## After you have real surface points

If you export a CSV with `x,y,z` columns, run:

```bash
python3 tests/diffmasif_replacement/compare_with_naccess.py \
  --structure-kind target \
  --structure-id 1FGNH \
  --surface-points /path/to/diffmasif_points.csv \
  --assignment-cutoff 4.0 \
  --output tests/diffmasif_replacement/output/1FGNH.diffmasif.json
```

If MaSIF gives you a `.ply` surface first, convert it with:

```bash
python3 tests/diffmasif_runtime/export_ply_vertices_to_csv.py \
  --ply /path/to/4ZQK_A.ply \
  --output /path/to/4ZQK_A.surface_points.csv
```

## Batch benchmark

After activating `masif_py36` and sourcing `masif_env_masif_py36.sh`, you can
run a small multi-target benchmark like this:

```bash
python3 tests/diffmasif_runtime/batch_masif_compare.py \
  --masif-root /scratch/rshadi25/GitHub/masif \
  --structure-id 1FGNH \
  --structure-id 1TFHA \
  --structure-id 1TFHB \
  --structure-id 2AI9A \
  --structure-id 2AI9B
```

This writes:

- `tests/diffmasif_replacement/output/masif_batch_summary.json`
- `tests/diffmasif_replacement/output/masif_batch_summary.csv`
