# PRISM-prescript CLI Reference

This reference describes the parser and runtime wiring on local branch
`feature/prism-cli-parity` at compatibility commit `20ebf4c3cc3`. It separates
three different claims:

- **Parser-supported:** `prism.py:build_parser()` accepts the syntax.
- **Runtime-wired:** `prism.py:main()` passes the value to its backend.
- **Scientifically validated:** a retained run establishes the behavior under
  a declared environment and input contract.

Parser compatibility does not by itself establish end-to-end or scientific
validation.

## Working directory and interpreter

Run commands from:

```bash
cd /scratch/rshadi25/GitHub/PRISM-prescript
```

Use the validated interpreter:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py [OPTIONS]
```

The examples below shorten that path to `python` for readability. For actual
runs, retain the explicit interpreter unless the active environment has been
verified.

Inspect the live parser without executing the pipeline:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --help
```

## Basic options

### Select input CSV

Accepted spellings:

```bash
python prism.py --inputs_csv inputs.csv
python prism.py --inputs-csv inputs.csv
```

The CSV must contain `Receptor` and `Ligand` columns.

Status: **parser-supported, partially runtime-wired**. The selected path is
passed to `src.pdb_download.py:pdb_downloader()`. However,
`src/transformation.py:transformer()` still reads its module-level
`PRISM_INPUTS_CSV`/`inputs.csv` path. Until propagation is unified, a custom
CLI path can make download and transformation read different CSVs.

Safe current workaround for an isolated run:

```bash
PRISM_INPUTS_CSV=/absolute/path/pairs.csv \
python prism.py --inputs_csv /absolute/path/pairs.csv [OTHER OPTIONS]
```

Use the same absolute path in both controls and record it in the run manifest.
Do not modify the repository's shared `inputs.csv` merely to route a test.

### Generate templates

Bare flag:

```bash
python prism.py --generate_templates
python prism.py --generate-templates
```

Explicit boolean, retained for compatibility:

```bash
python prism.py --generate_templates true
python prism.py --generate_templates false
```

The default is `false`.

### Limit templates

```bash
python prism.py --template-limit 100
python prism.py --template_limit 100
```

If omitted, prescript uses the complete `templates/calculated_templates.txt`
manifest. A supplied limit must be positive.

## Surface extraction

### NACCESS

NACCESS remains the prescript default:

```bash
python prism.py --surface-backend naccess
python prism.py --surface_backend naccess
```

### FreeSASA

```bash
python prism.py \
  --surface-backend freesasa \
  --freesasa-python /path/to/python-with-freesasa
```

Both spellings work:

```text
--surface-backend / --surface_backend
--freesasa-python / --freesasa_python
```

NACCESS is the stable reference backend. Treat FreeSASA as an explicit
alternative and compare surface/candidate outputs before interpreting yield
differences.

## Structural alignment

### TMalign

TMalign is the default:

```bash
python prism.py --aligner tmalign
```

The executable remains configurable through `PRISM_TMALIGN`; there is no
TMalign-path CLI option in the current parser.

### GTalign

```bash
python prism.py \
  --aligner gtalign \
  --gtalign-path /home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu
```

Optional controls:

```bash
--gtalign-dev-min-length 3
--gtalign-pre-score 0.0
--gtalign-speed 0
--gtalign-refinement 3
```

Underscore aliases also work:

```text
--gtalign_dev_min_length
--gtalign_pre_score
--gtalign_speed
--gtalign_refinement
```

The path and controls are passed to `src.alignment_gtalign.py:align_gtalign()`.
GTalign CPU and GPU outputs are not currently authorized as interchangeable.

### MultiProt

```bash
python prism.py \
  --aligner multiprot \
  --multiprot-path external_tools/multiprot.Linux \
  --multiprot-workers 8
```

- `--multiprot-path` selects the executable at call time.
- `--multiprot-workers` controls the bounded alignment worker count.

Both values are runtime-wired to
`src.alignment_multiprot.py:align_multiprot()`. MultiProt's legacy binary and
runtime constraints, fragment-length bias, and separate score/gate contract
still apply.

## Refinement

Prescript refines by default with external Rosetta. CLI spelling parity did
not change this repository-specific default.

### External Rosetta

```bash
python prism.py \
  --refine \
  --refiner external_rosetta
```

`--refine` is optional because refinement defaults to enabled. Rosetta module,
database, and executable configuration must already be valid in the batch
runtime.

### Skip refinement

```bash
python prism.py --no-refine
```

This stops after transformation or optional ranking. It is useful for
alignment/transformation diagnostics and for preparing candidates for a
separate controlled refinement run.

`--no-refine --compare` is generally not a meaningful fresh-run combination:
the comparison stage expects refined output conventions. Use comparison only
when the intended scoreable models are explicitly staged and identified.

### PyRosetta

```bash
python prism.py \
  --refine \
  --refiner pyrosetta \
  --pyrosetta-output-dir processed/pyrosetta_refinement \
  --pyrosetta-init-options='-mute all -constant_seed -jran 12345'
```

Use `=` when initialization options begin with `-`; otherwise `argparse` may
treat their tokens as PRISM options. The output root and initialization string
are runtime-wired to `src.pyrosetta_refinement.py:refine_pairs()`.

### FiberDock

```bash
python prism.py \
  --refine \
  --refiner fiberdock \
  --fiberdock-dir external_tools/fiberdock
```

The directory is applied before calling
`src.fiberdock_refinement.py:refine_pairs()`. FiberDock requires its binary,
NMA, Reduce/hydrogen-generation helpers, and parameter scripts under the
supplied runtime contract. Historical equivalence remains unresolved.

## Candidate ranking

Ranking is disabled by default. Unlike the new bare boolean flags described
above, `--rank` currently requires an explicit boolean value.

### Baseline ranking

```bash
python prism.py \
  --rank true \
  --rank-method baseline \
  --top-k 5
```

Optional minimum score:

```bash
--rank-min-score 0.25
```

The deterministic baseline keeps the best candidates per receptor–ligand
group. It is an opt-in resource-selection experiment, not established quality
improvement.

### PRODIGY ranking

```bash
python prism.py \
  --rank true \
  --rank-method prodigy \
  --top-k 5 \
  --prodigy-executable /path/to/prodigy
```

Optional controls:

```bash
--prodigy-output-dir processed/ranking/prodigy
--prodigy-distance-cutoff 5.5
--prodigy-acc-threshold 0.05
--prodigy-temperature 25.0
--prodigy-timeout 120
```

Select an explicit candidate-audit path when required:

```bash
--candidate-audit-path processed/candidate_audit/run.jsonl
```

Without it, ranked runs receive a generated run-scoped audit path. Disable
ranking explicitly with:

```bash
--rank false
```

## DockQ comparison

Comparison is opt-in. The bare flag and explicit boolean forms all parse:

```bash
python prism.py --compare
python prism.py --compare true
python prism.py --compare false
```

Parallel comparison:

```bash
--compare-jobs 8
```

Skip DockQ's internal structural alignment with either form:

```bash
--dockq-no-align
--dockq-no-align true
```

Use `--dockq-no-align` only when model and native structures are already in
the intended common frame.

## Combined parser-compatible example

GTalign, FreeSASA, PRODIGY top-3 ranking, PyRosetta, and DockQ:

```bash
PRISM_INPUTS_CSV=/absolute/path/inputs.csv \
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py \
  --inputs_csv /absolute/path/inputs.csv \
  --surface-backend freesasa \
  --freesasa-python /home/rshadi25/.conda/envs/gtalign_env/bin/python \
  --aligner gtalign \
  --gtalign-path /home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu \
  --template-limit 100 \
  --rank true \
  --rank-method prodigy \
  --top-k 3 \
  --prodigy-executable /path/to/prodigy \
  --refine \
  --refiner pyrosetta \
  --pyrosetta-init-options='-mute all -constant_seed -jran 12345' \
  --compare \
  --compare-jobs 4
```

This command is **parser-validated**, not a universal ready-to-run recipe. It
still requires staged inputs/templates, a working FreeSASA interpreter,
GTalign, PRODIGY, PyRosetta, native structures/mappings for comparison, and an
isolated output/evidence contract.

Submit computationally significant pipeline runs through Slurm. Do not launch
alignment, refinement, scoring, or broad template loops directly on a login
node.

## Version identification

Read-only local inspection:

```bash
git branch --show-current
git show --stat 20ebf4c3cc3
```

Expected active branch for this reference:

```text
feature/prism-cli-parity
```

The local feature branch currently has no configured upstream. Local Git
metadata alone therefore does not prove whether equivalent commits exist on a
remote. This reference intentionally contains no history-reversal commands.

## Parser verification

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_prism_cli_parity.py
```

Also inspect current help after any CLI update:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --help
```
