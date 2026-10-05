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

For the complete maintained callable surface, see the generated
[PRISM function inventory](PRISM_FUNCTION_INVENTORY.md). Regenerate it after
adding or moving current `src/` functions:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  tools/build_pipeline_function_inventory.py
```

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

Status: **parser-supported and runtime-wired**. The selected path is passed to
both `src.pdb_download.py:pdb_downloader()` and
`src/transformation.py:transformer()`. The environment variable
`PRISM_INPUTS_CSV` remains a default fallback for callers that invoke the
transformer directly.

Use a separate CSV for an isolated comparison without modifying the shared
`inputs.csv`:

```bash
python prism.py --inputs-csv /absolute/path/pairs.csv [OTHER OPTIONS]
```

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

### Select templates explicitly

Use either a plain-text manifest or explicit six-character template IDs:

```bash
python prism.py --template-list /absolute/path/template_subset.txt --no-refine
python prism.py --templates 1zx4AB 1ngmEF --no-refine
python prism.py --templates 1zx4AB,1ngmEF --template-limit 1 --no-refine
```

The two selection options are mutually exclusive. Explicit selection replaces
the generated/default panel, and `--template-limit` is applied after
selection. The current template contract is a four-character PDB ID followed
by two chain IDs, for example `1zx4AB`.

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

## Orientation comparison

The transformation stage supports native and fixed-orientation comparison
runs through `--orientation`:

```bash
# Native MultiProt-like behavior: evaluate both implicit assignments.
python prism.py --orientation native --no-refine

# Fixed-assignment comparison runs.
python prism.py --orientation o1 --no-refine
python prism.py --orientation o2 --no-refine
```

Omitting `--orientation` deterministically selects `native`, which evaluates
both assignments (`o1` and `o2`),
retaining their labels in transformation filenames and candidate audits.
`o1` and `o2` constrain only the transformation-stage branch; all other
inputs, thresholds, aligner settings, and refinement settings should be held
constant when comparing runs. This is an orientation restriction, not an
instruction to ignore template-chain assignment, which would invalidate
chain-specific contacts and transforms.

Preserve each run in a separate isolated output/work directory when comparing
generated PDBs: the stable transformation filenames include `o1`/`o2`, but a
fixed `o1` run can otherwise overwrite the native run's `o1` artifacts.

Status: **parser-supported and runtime-wired**. The comparison contract is
covered by focused tests; scientific equivalence to every legacy run remains
an empirical validation question.

For the recommended diagnostic order, use
`notebooks/pipeline_orientation_comparison.ipynb`. It freezes the production
baseline, records asset/source hashes and orientation partner availability,
builds alignment and independent gate ledgers, replays the clash grid, checks
refinement evidence, compares ranking methods only on an explicitly supplied
same-set candidate panel, and keeps US-align as a preflight-only separate arm.
When native-like recovery is supplied, the notebook also requires an
independent native-label source path and SHA-256; PRODIGY affinity is never a
native-like label. The manifest includes the template interface-list JSON used
for chain-specific coverage; incomplete alignment return/status or raw-output
hash evidence keeps dependent gate outcomes unknown.

## Transformation thresholds

Transformation thresholds are environment-backed for compatibility and can be
overridden per run from the CLI. Omitted options retain the current defaults:

```text
minimum-residue-match-count       15
minimum-residue-match-percentage  50.0
minimum-hotspot-match-number      1
diff-percentage                   20.0
template-residue-count            50
contact-count-threshold           5
clashing-distance                 3.0 Å
max-clashing-count                5
scaffold-threshold                5.0 Å
tm-score-threshold                0.5
MultiProt match count / coverage   10 / 30.0%
alignment-gate-mode               native
```

All controls are available as hyphenated CLI options, for example:

```bash
python prism.py \
  --minimum-residue-match-count 12 \
  --minimum-residue-match-percentage 35 \
  --diff-percentage 25 \
  --clashing-distance 2.5 \
  --max-clashing-count 8 \
  --scaffold-threshold 5.0 \
  --tm-score-threshold 0.4 \
  --no-refine
```

The MultiProt-specific controls are
`--multiprot-minimum-residue-match-count` and
`--multiprot-minimum-residue-match-percentage`. These are separate from the
TMalign/GTalign TM-score gate because MultiProt uses native correspondence
count and interface coverage. Threshold changes are diagnostic until assessed
against a fixed benchmark with raw, qualified, contact, clash, and downstream
outcomes recorded.

### Comparable cross-aligner mode

For a controlled TMalign/USalign/GTalign/MultiProt threshold comparison, use:

```bash
python prism.py \
  --alignment-gate-mode common_match_coverage \
  --minimum-residue-match-count 15 \
  --minimum-residue-match-percentage 50 \
  --diff-percentage 20 \
  --orientation native \
  --no-refine
```

`common_match_coverage` applies the same matched-residue count and size-adjusted
interface-coverage gates to every aligner, using an inclusive boundary. It
does not apply `tm-score-threshold`, because MultiProt's RMSD-derived score is
not a TMalign-compatible TM-score. Provider scores remain in the alignment
records for post-hoc analysis. The default `native` mode is unchanged.

Candidate-audit rows also retain the resolved threshold dictionary and use
distinct terminal statuses for missing alignment, protocol rejection,
alignment-threshold rejection, transformation failure, and clash rejection.
The scaffold threshold is applied during surface extraction; its stable value
is 5.0 Å and diagnostic overrides should be recorded with the run manifest.

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
improvement. For any native-like recovery claim, provide labels from an
independent native reference and retain its source path plus SHA-256 in the
comparison ledger.

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
