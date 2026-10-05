# PRISM current status — 2026-09-13

## Bottom line

The maintained canonical pipeline is operationally implemented and its focused
software tests/smokes pass, but no clean same-input/same-template/same-evaluator
comparison has scientifically validated USalign, PRODIGY, or a current exact
19K panel. Existing high DockQ and runtime results are historical observational
evidence, not a promotion-grade causal matrix.

Canonical source: `/scratch/rshadi25/GitHub/PRISM-prescript` at `1026a7fc0609e19f19c48a68510e7ca39d67573d`, branch `feature/prism-cli-parity`, with `380` dirty entries. It was not edited.

## Evidence state

| State | Meaning in this package |
|---|---|
| IMPLEMENTED | Code/path exists in canonical or an isolated candidate. |
| TESTED | Focused tests or deterministic probes passed. |
| VALIDATED | Runtime and quantitative output evidence exists within the stated scope. |
| REVIEWED | Source/evidence review completed; this does not upgrade scientific validation. |
| NOT_COMPARABLE / NOT_RUN | Inputs, evaluator, provenance, or execution are insufficient for causal comparison. |

Matrix row counts: {"COMPLETED": 21, "FAILED": 1, "NOT_COMPARABLE": 8, "NOT_RUN": 29}.

## Implemented and tested

* TMalign is the maintained default; GTalign CPU/GPU are opt-in; NACCESS and
  FreeSASA are surface choices; external Rosetta is the stable refiner.
* USalign and PRODIGY changes exist only in isolated worktrees. USalign's
  transform-producing invocation is `-outfmt -1 -m -`; its bounded real probe
  produced 214 matches, TM-score 0.98359, and replay RMSD 0.7500643.
* PRODIGY 2.4.0 bounded success selected one candidate (affinity -65.827) and
  its no-contact failure preserved the full candidate group. This is state
  handling, not docking-quality validation.
* The isolated shared-contract candidate passed 28 focused contract/adapter/
  import tests and 31 tests in total. It fail-closes unavailable optional
  backends and preserves unknown/rejected alignment states; it is not wired
  into canonical provider writers and is not promoted.
* Focused canonical/isolated tests recorded previously were 22/22, 20/20,
  and 9/9 respectively. The full canonical suite recorded 344 passed and 6
  skipped in the latest retained run.

## Runtime and quality evidence

* Logical 946-template alignment panel: 1 query × 946 templates × 2 sides =
  1,892 pairs; TMalign 22.27–23.08 s and GTalign 72.44–95.96 s across two
  sites. All pairs were present; mean absolute TM-score difference was about
  0.030. GTalign device was not recorded, so this is not GPU evidence.
* Historical 19,855-template BM5.5 run: GTalign GPU wall time 906 s, 60/64
  transforms, 20/22 scored models, mean best DockQ 0.883/0.880; TMalign had
  1,474 reported transforms and much larger alignment output. GTalign CPU
  completed with no predictions (27 alignment JSONs). These are historical
  artifacts with incomplete parity and known mapping/CPU-GPU caveats.
* Current exact on-disk denominators are 19,948 checked entries (SHA-256
  `740f68deae299995a016bf2f83abeb036713189d7e9e78b119c8a556f85e10a0`) and 19,062 calculated entries (SHA-256
  `7277d71a864be57e89fd50a8aeae13db06abbe6dc3d82e233941fc83936cc0aa`). The 946 panel is logical first-946 selection from a
  historical 19,855-entry list (selected-panel SHA-256 `9ac528a5de6bb03ae71893aea6c5ab339f5a1f4038c35749ca165ce8d005918f`). The new
  current checked-prefix 946 runs use a distinct selected-panel SHA-256
  `b77212e0aee62357bf0b2885105e709479f4d7a06907ef6162f08a6ab3d75a20`.

## New exact-panel alignment-only runs

| Arm | Status | Selected templates | Staging seconds | Search seconds | Return code |
|---|---|---:|---:|---:|---:|
| exact-checked_prefix_946-gtalign-gpu-1659356 | COMPLETED | 946 | 3.069310188293457 | 4.813395977020264 | 0 |
| exact-checked_full_19948-gtalign-gpu-1659357 | FAILED | 19948 |  |  | 1 |
| exact-calculated_19062-gtalign-gpu-1659358 | COMPLETED | 19058 | 70.65782117843628 | 73.95956349372864 | 0 |
| exact-checked_prefix_946-tmalign-1659448 | COMPLETED | 946 | 0.829352617263794 | 11.476548194885254 | 0 |
| exact-calculated_19062-tmalign-1659449 | COMPLETED | 19058 | 12.218315124511719 | 144.9500024318695 | 0 |
| exact-checked_prefix_946-usalign-1659360 | COMPLETED | 946 | 0.7320449352264404 | 21.49501657485962 | 0 |
| exact-calculated_19062-usalign-1659361 | COMPLETED | 19058 | 12.126771211624146 | 439.97117948532104 | 0 |

These jobs use one query and precomputed surface PDBs. They validate search
execution, transform-producing output/return handling, and stage timing only;
they do not validate transformation filtering, refinement, ranking, DockQ, or
scientific superiority.

## Decision at this checkpoint

* Fastest observed completed arm: historical GTalign GPU + external Rosetta at
  906 s on the 19,855-entry BM5.5 run. It is not yet a scientifically
  defensible promoted default because current-panel parity and stage timing are
  incomplete.
* Best observed quality subset: historical GTalign GPU + external Rosetta,
  mean best DockQ 0.928 (12 scored models); GTalign GPU + PyRosetta was 0.925
  (13 scored models). The small scored subset and mapping gaps limit this claim.
* Recommended reproducible default today: retain the documented TMalign +
  NACCESS + external Rosetta recipe. Use GTalign GPU + external Rosetta as the
  provisional large-library speed candidate after exact-panel validation.
* GTalign GPU is provisionally the preferred large-library search engine from
  observed runtime/output evidence; USalign is not yet a stable alternative.
  PRODIGY has no demonstrated end-to-end compute saving or acceptable quality
  retention. A ~20K search is operationally demonstrated historically but not
  yet promotion-ready on the current exact panel.
* Highest-value optimization: instrument and reduce repeated alignment output,
  structure parsing, and staging while proving candidate-set agreement.

## Blockers and uncertainty

The maintained path lacks a uniform AlignmentResult/no-drop schema across all
providers; external Rosetta does not preserve complete per-candidate return
codes; transformation audit does not always contain actual clash counts;
MultiProt's RMSD-derived proxy is not a TM-score; exact current-panel GTalign
GPU, TMalign, and USalign evidence is alignment-only, while current exact
transformation/filter, ranking, refinement, and evaluator arms are absent; and
no PRODIGY end-to-end timing or quality retention evidence exists. Job
reconciliation is in `RECONCILIATION.md`. The isolated contract candidate
remains partial because it lacks event/parent correlation, does not wrap every
stage, and is not connected to provider writers. Existing canonical
`run_evidence`, lineage, evaluator, and completion contracts remain
authoritative; the candidate is an additive projection and must not replace
them.

## Authoritative artifacts

* `PIPELINE_MATRIX.csv/json`: all retained, planned, incomplete, and historical arms.
* `TEMPLATE_PANELS.json` and `NO_DROP_SUMMARY.json`: frozen panel definitions and retained no-drop accounting.
* `PERFORMANCE_COMPARISON.csv`: exact 946 alignment evidence, historical 19,855 results, and explicit current-panel gaps.
* `RANKING_COMPARISON.csv`: bounded baseline/PRODIGY evidence and required missing no-ranking/top-k arms.
* `notebooks/prism_pipeline_comparison.ipynb`: portable artifact-driven report view.
* `PIPELINE_COMPARISON.md`, `STABLE_PIPELINES.md`, `VALIDATION_REPORT.md`, `NEXT_ACTIONS.md`, and `CANDIDATE_CHANGES.md`.
