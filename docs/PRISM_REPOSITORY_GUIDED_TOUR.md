# PRISM-prescript Repository Guided Tour

## Purpose and audience

This is a presenter-ready onboarding tour for researchers and junior developers.
It explains what PRISM-prescript does, how a command moves through the code,
which external tools perform the scientific work, and how to distinguish a
created file from a validated result.

The central question for the audience is:

> If a docking result appears under `processed/`, what evidence tells us which
> inputs, tools, filters, and execution attempt produced it?

The stable reference route is **NACCESS + TMalign + external Rosetta**. Other
backends and ranking modes are explicit alternatives and must not be presented
as equivalent or better without matched validation.

## Choose a tour length

### 15-minute orientation

1. Show the repository landmarks and `python prism.py --help`.
2. Trace `prism.py:main()` through the seven pipeline stages.
3. Open one alignment JSON, one transformed PDB, and one stage-status JSONL
   from a retained run.
4. Explain stable defaults versus opt-in capabilities.
5. End with the Phase 1 provenance consumer-gate problem.

### 45–60-minute guided session

| Time | Segment | Outcome |
| --- | --- | --- |
| 0–5 min | Scientific goal and repository map | Audience can state what PRISM predicts. |
| 5–12 min | CLI and configuration | Audience can identify defaults and opt-ins. |
| 12–30 min | Trace one receptor–ligand pair | Audience can connect stages to code and files. |
| 30–40 min | External tools and HPC execution | Audience knows what is safe to inspect locally and what needs Slurm. |
| 40–50 min | Validation, audits, and provenance | Audience can explain why file counts are insufficient. |
| 50–60 min | Current work and learning check | Audience can separate validated, experimental, and unfinished behavior. |

## Presenter preflight

Run lightweight inspection commands from the repository root. Do not launch
alignment, refinement, scoring, downloads, or batch loops on a login node.

```bash
pwd
git status --short
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --help
```

Expected observations:

- The CLI imports successfully under `gtalign_env`.
- The default aligner is `tmalign`, default surface backend is `naccess`, and
  default refiner is `external_rosetta`.
- Ranking and comparison are opt-in.
- The worktree may be dirty. Do not imply that retained results came from the
  current bytes unless their source manifest or hashes establish that link.

Prepare retained evidence before the session. Prefer a completed isolated run
under `tmp/agent/` with its launcher, `run.log`, stage-status JSONL, candidate
audit when applicable, and output files. Avoid depending on network downloads
or a live Rosetta run during the presentation.

## Repository landmarks

| Path | What to say | What to show |
| --- | --- | --- |
| `prism.py` | User-facing CLI and stage orchestration. | `main()`, parser definitions, `run_stage()`. |
| `src/` | Pipeline adapters, filters, transformations, refinement, evaluation, ranking, and provenance. | Follow the modules in the stage map below. |
| `templates/` / `new_template/template/` | Interface structures and derived template assets. | `pdbs`, `interfaces`, `interfaces_lists`, `contacts`, `hotspots`, `rsas`. |
| `processed/` | Mutable runtime outputs, not proof of completion by itself. | Alignment JSON, transformed models, refinement folders. |
| `benchmark/scripts/` | Reproducible launchers, scorers, audits, collectors, and diagnostics. | Stable smoke launcher and selected validators. |
| `tests/` | Unit, contract, and regression evidence. | Tests that exercise CLI, transformation, ranking, and provenance. |
| `.planning/` | Active reliability/provenance scope and decisions. | `PROJECT.md`, `STATE.md`, roadmap, Phase 1 artifacts. |
| `.agents/.../project-memory/references/` | Durable validated context and unresolved questions. | `summary.md`, `decisions.md`, `open_questions.md`. |
| `working_version/` | Retained legacy compatibility/reference tree. | Mention only; do not claim current equivalence. |

Ask the audience: **Which of these paths contains executable behavior, and
which contains evidence or project intent?** The useful distinction is source
code versus mutable outputs versus durable decisions.

## The pipeline in one view

```text
inputs.csv
   |
   v
PDB download/materialization --> optional template generation
   |
   v
surface extraction (NACCESS default; FreeSASA opt-in)
   |
   v
structural alignment (TMalign default; GTalign/MultiProt opt-in)
   |
   v
alignment/coverage/protocol gates + transform + clash filtering
   |
   v
optional candidate ranking (baseline or PRODIGY)
   |
   v
refinement (external Rosetta default; PyRosetta/FiberDock opt-in)
   |
   v
optional DockQ/iRMSD comparison and evidence audit
```

`prism.py:main()` is the best place to present this sequence. `run_stage()`
wraps selected stages and can append started/completed/failed records when
`PRISM_STAGE_STATUS_PATH` is configured. Surface extraction and optional
template generation are not yet uniformly wrapped by that stage ledger; this
is one reason the current reliability work is not finished.

## Command-to-code map

| Command or control | Functionality | Implementation section | Tool/backend | Principal output | Status and presentation note |
| --- | --- | --- | --- | --- | --- |
| `python prism.py --help` | Lists the current public CLI. | `prism.py`, parser block after `main()` | Python/argparse | Terminal help | Safe live demo. |
| `python prism.py` | Runs the stable default orchestration. | `prism.py:main()` | NACCESS, TMalign, external Rosetta | `processed/*` | Heavy/network-capable; use Slurm or retained evidence. |
| `PRISM_INPUTS_CSV=<csv>` | Selects the receptor–ligand table without editing `inputs.csv`. | `src/pdb_download.py:pdb_downloader()` and `src/transformation.py:transformer()` | pandas | Chain-qualified PDBs and pair loop | CSV requires `Receptor,Ligand`. |
| `--generate_templates true` | Analyzes PDB templates and regenerates interface assets. | `src/analyse_pdbs.py:run_analysis()`; `src/template_generate.py:template_generator()` | NACCESS and template utilities | Template assets and `calculated_templates.txt` | Optional; not part of the bounded default smoke route. |
| `--template-limit N` | Uses only the first N calculated templates. | `prism.py:main()` template-loading block | Python | Bounded downstream workload | Intended for smoke tests; N must be positive. |
| `--surface_backend naccess` | Computes relative accessibility and CA surface scaffolds. | `src/surface_extract.py:extract_surfaces()`; `src/naccess_utils.py` | `external_tools/naccess/naccess` | `processed/surface_extraction/*.asa.pdb`; `failures.tsv` | Stable default. |
| `--surface_backend freesasa --freesasa_python <python>` | Uses FreeSASA through a separate interpreter. | `src/naccess_utils.py`; `src/freesasa_runner.py` | FreeSASA | Same surface-stage contract | Explicit alternative; compare outputs before interpreting yield differences. |
| `--aligner tmalign` | Aligns each query surface to template-interface chains. | `src/alignment.py:align()` and `parse_tmalign()` | `external_tools/TMalign` or `PRISM_TMALIGN` | `processed/alignment/*.json` | Stable default. |
| `--aligner gtalign --gtalign_path <binary>` | Runs GTalign and parses hits into alignment records. | `src/alignment_gtalign.py:align_gtalign()` | GTalign CPU/GPU | `processed/alignment_gtalign/<run-id>/*.json` plus compatibility symlink | Opt-in; CPU/GPU divergence remains unresolved. Symlink statefulness is an open contract issue. |
| `--aligner multiprot` | Runs MultiProt, extracts matches, calculates transforms and true TM-score. | `src/alignment_multiprot.py:align_multiprot()` | `external_tools/multiprot.Linux` | `processed/alignment/*.json` | Diagnostic/opt-in; short-fragment bias and legacy runtime equivalence remain open. |
| transformation controls such as `PRISM_TM_SCORE_THRESHOLD` | Applies alignment-specific gates, transforms both partners, and rejects clashes. | `src/transformation.py:transformer()`, `alignment_score_passes()`, `alignment_passes_thresholds()`, `process_pair_for_template()` | NumPy/Biopython and template assets | `processed/transformation/*.pdb`; optional candidate-audit JSONL | Defaults are scientific baseline controls; diagnostic overrides are not production defaults. |
| `PRISM_FILTER_MODE=published_protocol` | Enables protocol asset/hotspot filtering rather than geometry-only experimental mode. | `src/transformation.py`; `src/template_filtering.py` | contacts/hotspots/RSA assets | Candidate decisions in audit | Asset parity remains unresolved; label the mode clearly. |
| `--rank true --top-k N` | Ranks accepted candidates before refinement. | `src/candidate_selector.py:select_top_candidates()`; `src/candidate_ranker.py` | Deterministic baseline | Run-scoped `processed/candidate_audit/*.jsonl`; reduced candidate list | Opt-in resource-reduction experiment, not proven quality improvement. |
| `--rank true --rank-method prodigy --prodigy-executable <path>` | Scores transformed complexes by predicted affinity and keeps top candidates. | `src/prodigy_ranker.py`; called through `candidate_selector.py` | PRODIGY executable | `processed/ranking/prodigy/` inputs, stdout/stderr, JSON records | Opt-in and paired-smoke validated mechanically; independent biological evaluation remains open. |
| `--refiner external_rosetta` | Pre-packs, docks, score-filters, and assembles candidate models. | `src/rosetta_refinement.py:refiner()` | Rosetta 2022.42 executables/database | `processed/rosetta_refinement/` | Stable default, but per-candidate return/score-gate observability is incomplete. |
| `--refiner pyrosetta` | Refines candidates through the PyRosetta API. | `src/pyrosetta_refinement.py:refine_pairs()` | PyRosetta | `processed/pyrosetta_refinement/` | Explicit alternative; never an implicit fallback. |
| `--refiner fiberdock` | Adds hydrogens, generates modes/parameters, runs FiberDock, and records energies. | `src/fiberdock_refinement.py:refine_pairs()` | Reduce, NMA, FiberDock, helper scripts | `processed/fiberdock_refinement/structures` and `energies` | Opt-in; historical equivalence and parser integration remain unresolved. |
| `--compare true [--compare-jobs N]` | Scores refined outputs against native structures. | `src/compare.py:compare_pairs_from_outputs()`; `src/eval/` | DockQ and iRMSD code | `processed/summary.csv`; trimmed natives | Opt-in; mapping and hash identity matter more than filename presence. |
| `PRISM_STAGE_STATUS_PATH=<jsonl>` | Records stage start/completion/failure events. | `prism.py:record_stage_event()` and `run_stage()` | JSONL writer | User-selected JSONL | Useful evidence, but current lifecycle vocabulary/coverage is not yet the full Phase 2 contract. |
| `python -m src.provenance.run_evidence --help` | Shows the Phase 1 contract/attempt/artifact-ledger CLI. | `src/provenance/run_evidence.py:main()` | Standard-library provenance module | Contract JSON, manifest JSON, ledger TSV, closeout/validation records | Ongoing Phase 1 work; do not present as accepted integration yet. |

### Boolean CLI pitfall

The boolean options use a custom `parse_bool()`. For clarity in a presentation,
write explicit values such as `--rank true` and `--compare true`. The bare
`--rank` spelling shown in some older documentation is not accepted by the
current parser because the option expects a value.

## Trace one pair through the source

Use one retained receptor/ligand/template case. At every stop, ask the audience
to predict the next filename or record before revealing it.

### 1. Input and chain materialization

Open `src/pdb_download.py:pdb_downloader()`.

- Reads `PRISM_INPUTS_CSV` or `inputs.csv`.
- Normalizes PDB-plus-chain selectors such as `3i6eEF`.
- Downloads a four-character PDB only when it is not staged already.
- Writes chain-qualified structures under `processed/pdbs/`.

Question: **Why is the chain suffix scientifically important even though the
download URL uses only the four-character PDB ID?**

Expected answer: a benchmark row identifies molecular roles and selected
chains, while the archive entry may contain additional chains that must not be
silently included.

### 2. Surface scaffold

Open `src/surface_extract.py:extract_surfaces()` and `extract_surface()`.

- The default NACCESS route calculates relative solvent accessibility.
- Residues above the RSA threshold seed a nearby CA scaffold.
- The stage preserves an explicit empty PDB plus `failures.tsv` information
  for some failures instead of silently omitting the target.

Show, do not generate live:

```text
processed/pdbs/<target>.pdb
processed/surface_extraction/<target>.asa.pdb
processed/surface_extraction/failures.tsv
```

Question: **Is an `END`-only surface PDB success, failure, or a data outcome?**

Expected discussion: the file preserves a downstream contract, but its stage
status and failure record determine how it should be interpreted.

### 3. Structural alignment

Open the backend selected by `prism.py`:

- `src/alignment.py:align()` for TMalign.
- `src/alignment_gtalign.py:align_gtalign()` for GTalign.
- `src/alignment_multiprot.py:align_multiprot()` for MultiProt.

In an alignment JSON, locate the aligner name, TM-score contract, match count,
match dictionary, and transformation values. Do not compare a field named
`tm_score` across backends without checking its contract and normalization.

Question: **If two adapters both emit JSON, does that prove they implement the
same schema?**

Expected answer: no. A formal shared alignment-result contract is still an
open issue.

### 4. Gates, transformation, and clashes

Open `src/transformation.py` in this order:

1. Module-level default thresholds.
2. `alignment_score_passes()` for aligner-specific scoring.
3. `alignment_passes_thresholds()` for count, coverage, and hotspot gates.
4. `process_pair_for_template()` for orientation, transformation, clash
   accounting, output creation, and candidate audit.

Explain that transformed structures exist only after both partner alignments
and downstream geometry checks. Diagnostic environment overrides can change
candidate yield, so they belong in a recorded run contract.

Question: **What would be lost if we reported only the number of transformed
PDB files?**

Expected answer: the rejected denominator and the reason each candidate failed.

### 5. Optional ranking

Open `src/candidate_selector.py:select_top_candidates()`.

The baseline uses retained alignment/audit evidence. PRODIGY instead invokes
an external scorer and retains its command, input hash, output, return code,
and failure status. If a PRODIGY group is only partially scoreable, the group
is preserved rather than silently selecting from incomplete evidence.

Question: **Does reducing two candidates to one prove that ranking improved
docking quality?**

Expected answer: no. It proves candidate forwarding/load reduction; quality
requires frozen native labels and independent evaluation.

### 6. Refinement

Follow the selected branch from `prism.py` to one of:

- `src/rosetta_refinement.py:refiner()`
- `src/pyrosetta_refinement.py:refine_pairs()`
- `src/fiberdock_refinement.py:refine_pairs()`

External Rosetta is the stable reference but requires the Rosetta module and
configured executables/database. The code currently does not retain enough
per-candidate subprocess and energy-gate detail to explain every missing
canonical model; this is an active observability gap.

### 7. Optional evaluation

Open `src/compare.py:compare_pairs_from_outputs()` and then `src/eval/dockq.py`.

Stress four identifiers:

- dataset row identity,
- model and native hashes,
- receptor/ligand chain mapping,
- requested interface scope.

A score with the wrong chain mapping can be numerically valid yet
scientifically attached to the wrong comparison.

## Demonstration commands

### Safe, lightweight live commands

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py --help

rg -n "^def (main|run_stage|record_stage_event)" prism.py

rg -n "^def (pdb_downloader|extract_surfaces|align|align_gtalign|align_multiprot|transformer|refiner|refine_pairs|compare_pairs_from_outputs)" \
  prism.py src

/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  -m src.provenance.run_evidence --help
```

These import or inspect code; they do not intentionally launch scientific
compute. The provenance help command currently emits a harmless `runpy`
warning because the package imports the module before `-m` executes it.

### Present from retained evidence

```bash
sed -n '1,160p' <run-root>/run.log
sed -n '1,80p' <run-root>/<stage-status>.jsonl
sed -n '1,3p' <run-root>/<candidate-audit>.jsonl
```

Replace placeholders with a reviewed isolated run. Avoid dumping large JSONL
files; show one record and ask the audience to interpret its identity, status,
and reason fields.

### Submit through Slurm; do not run as a live login-node demo

The canonical smoke launcher is:

```bash
PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  bash benchmark/scripts/run_prism_pipeline_smoke.sh
```

It creates an isolated workspace under `tmp/agent/`, stages one input pair and
one template, and records `run.log`. It can still download, align, and invoke
Rosetta, so submit it in an appropriate Slurm allocation when running the full
path. A zero-candidate smoke can validate mechanical pipeline health; it does
not establish biological quality.

For many independent tasks under the `ai` QoS limit, the project-standard
strategy is one Slurm allocation with bounded internal workers rather than one
submitted job per candidate. Record job ID, partition/QoS/account, resources,
working directory, source identity, inputs, tools, and task-level outcomes.

## Verification and evidence order

When reviewing a run, inspect evidence in this order:

1. Launcher/command and `run.log`.
2. Declared inputs, templates, source revision/diff, environment, and tool
   paths/versions.
3. Terminal stage-status records.
4. Candidate-audit JSONL when transformation auditing or ranking is enabled.
5. Raw alignment output and parsed alignment JSON.
6. Transformed/refined PDB integrity and exact artifact identity.
7. Raw evaluator JSON, chain mapping, score table, and denominator audit.

Never infer completion from directory or file counts. Cancelled and failed
jobs can leave plausible partial artifacts.

Lightweight focused checks suitable for demonstrating the test boundary:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_prism_cli.py \
  tests/test_transformation_thresholds.py \
  tests/test_candidate_selector.py \
  tests/test_run_evidence.py
```

Run tests only where lightweight compute is permitted. For a short talk, show
a previously recorded result rather than spending the session waiting.

## Maturity map: what is present versus what is proven

### Stable reference behavior

- Python 3 stage orchestration in `prism.py`.
- NACCESS surface extraction.
- TMalign structural alignment.
- External Rosetta refinement when its runtime is configured.
- Stable thresholds remain unchanged for the reference path.

### Mechanically validated but opt-in

- Deterministic top-K candidate ranking and run-scoped audit creation.
- PRODIGY candidate scoring/ranking on a retained paired smoke.
- DockQ comparison adapter and extensive benchmark scoring/audit utilities.
- GTalign, MultiProt, PyRosetta, FiberDock, and FreeSASA adapters, each with
  separate runtime and scientific caveats.

### Ongoing or unresolved

- MultiProt gate source/test drift: current `alignment_score_passes()` requires
  true TM-score `>= 0.3`, at least 10 matches, and at least 30% match coverage,
  while `tests/test_transformation_thresholds.py` still contains an older
  expectation that a MultiProt record with `tm_score = 0.0` passes based on
  match count. The focused tour verification currently reports 30 passes and
  this one failure; resolve the intended contract before presenting the suite
  as green.
- Phase 1 run identity and artifact-ledger integration. The core module and
  focused tests exist, but fail-closed CLI behavior and full `prism.py`
  integration are not accepted yet.
- Complete stage/candidate lifecycle observability, especially per-candidate
  external-Rosetta return codes and score-gate reasons.
- A formal shared alignment adapter schema and removal of GTalign symlink
  statefulness.
- GTalign CPU/GPU equivalence.
- MultiProt short-fragment biological interpretation and historical
  MultiProt/FiberDock parity.
- Independent evidence that candidate ranking improves biological quality.
- Authoritative resolution of 17 audit-only benchmark rows.

## Documentation caveats to call out

- `docs/STABLE_PIPELINE.md` currently includes a scoring command using
  `/scratch/tmp/prism-dockq-env/bin/python`. Project memory records that this
  environment no longer has DockQ. The verified repository-local scoring entry
  is `benchmark/prism_processed/env/prism_score_env/bin/python`; verify
  `import DockQ` before submission.
- `docs/multiprot_tmscore_analysis_report.md` preserves an earlier proxy-score
  diagnosis. Current source computes true TM-score and uses a MultiProt-specific
  gate. Present the report as historical investigation, not current behavior.
- `ReadMe.md` is useful background, but operational commands must be checked
  against `prism.py --help`, current source, and validated project memory.

## Presenter recovery notes

| Symptom | Ask first | Safe response |
| --- | --- | --- |
| Import fails during `prism.py --help` | Which interpreter is active? | Use the validated `gtalign_env` interpreter; do not install into shared environments during the talk. |
| No candidates reach refinement | Which stage first records rejection? | Inspect stage status, alignment JSON, and candidate audit before changing thresholds. |
| Rosetta executable is missing | Was `rosetta/2022.42` loaded and were `PRISM_ROSETTA_*` paths recorded? | Stop and fix the batch runtime; do not treat this as a docking-quality result. |
| MultiProt exits unexpectedly | What is the execution context and seccomp state? | Use approved Slurm/unrestricted context; do not substitute legacy helper versions and call them equivalent. |
| DockQ fails or scores few rows | Are model/native hashes and chain mappings complete? | Inspect raw evaluator errors and denominator audit; do not silently drop rows. |
| Files exist after a failed job | Is there a normal return and terminal record? | Treat them as partial artifacts until the evidence contract says otherwise. |

## Socratic learning checkpoints

Pause after each checkpoint and ask one participant to explain the answer in
their own words.

1. **CLI:** Which defaults define the stable scientific reference path?
2. **Identity:** Why is `1abcA` not interchangeable with `1abc`?
3. **Alignment:** What makes two equally named `tm_score` fields potentially
   incomparable?
4. **Filtering:** Where can a candidate disappear, and what record should
   explain that disappearance?
5. **Ranking:** What evidence distinguishes lower refinement load from better
   docking quality?
6. **Evaluation:** Why are row ID, hashes, mapping, and interface scope all
   needed beside a DockQ number?
7. **HPC:** What information must accompany a Slurm job ID for reproducibility?
8. **Provenance:** What should happen if a downstream consumer sees a changed
   artifact or mismatched row identity?

## Closing learning recap

- **Concept mastered:** PRISM is a file-oriented scientific pipeline whose
  orchestration, external-tool adapters, and evidence contracts are separate
  responsibilities.
- **Mistake to avoid:** treating file presence, candidate count, or one score
  as proof of successful and scientifically comparable execution.
- **Practice exercise:** choose one retained candidate and draw its lineage
  from `inputs.csv` through alignment JSON, transformation audit, refined PDB,
  and evaluator record. Mark every point where identity or status could be
  lost.
- **Success test:** a participant can locate the responsible function for a
  stage, name its external tool and output, and state what evidence is needed
  before trusting the result.
