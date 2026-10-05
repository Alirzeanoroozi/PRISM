# PRISM boundary and contract hardening

This ExecPlan is a living record for the current cross-repository integration
cycle. It must be updated after each meaningful implementation or validation
milestone.

## Purpose / Big Picture

Make `PRISM-prescript` the single maintained pipeline authority while safely
transferring only reviewed features from the experimental `PRISM` tree. The
result must preserve stable defaults, make backend failures observable, make
artifact identity auditable, and prevent incomplete evidence from becoming a
scientific success claim.

## Progress

- [x] Freeze both dirty repositories in run-scoped snapshot directories.
- [x] Record the repository boundary and canonical defaults in ADR-0003.
- [x] Repair the experimental `PRISM` documentation and candidate-path issues.
- [x] Complete fail-closed provenance validation in `PRISM-prescript`.
- [~] Complete stage and candidate lifecycle contracts.
- [x] Classify source, evidence, raw, and temporary files without destructive cleanup.
- [x] Run focused and full software validation.
- [x] Run isolated stable-baseline and backend smoke validations.
- [~] Complete orientation, FiberDock, MultiProt, Rosetta, and GTalign gates.
- [ ] Resolve benchmark authority issues and perform matched comparisons.
- [ ] Update project memory and reach the final release gate.

## Surprises & Discoveries

- The first snapshot attempt was interrupted; a second snapshot with run ID
  `20260912-workflow-freeze-retry1` was completed and preserved.
- The `PRISM-prescript` untracked manifest traverses nested generated and
  environment trees and is 256 MB uncompressed. It is retained as a gzip
  artifact rather than expanded repeatedly.
- Slurm inventory shows an externally running interactive job `1658869` on
  `ai12`; this plan must not cancel, resize, or reuse it.
- The current `PRISM` worktree deletes a guidance document still linked by its
  README. This must be repaired before release review.
- Optional-backend job 1658926 succeeded for PyRosetta and FiberDock but
  failed closed at DockQ because the selected environment's DockQ extension
  was compiled against NumPy 1.x while the environment provides NumPy 2.4.6.
- CPU and GPU GTalign completed on the same one-pair/template smoke, but their
  raw hashes, match counts, and some TM-scores differ; the arms remain
  provenance-separated.
- The first manifest CLI implementation expanded copied-environment files into
  a 262 MB contract. It was corrected to use bounded status entries plus
  content hashes; the corrected contract smoke produced 16 KB.
- The first template-inventory smoke followed a repository symlinked asset tree
  and was stopped before producing output. A bounded one-file inventory then
  completed successfully; linked template trees must be supplied through an
  explicit manifest or bounded directory.
- The evaluator integrity tests exposed a real compatibility gap in the dirty
  benchmark scorer: explicit chain overrides were ignored, the repository root
  was inferred from the current directory, and JSON output paths were not
  isolated. These were repaired and covered by the full suite.
- The controlled orientation notebook job 1658931 completed the three
  no-refinement pipeline arms in an isolated comparison root. Native default
  recorded two alignment-threshold rejections; `o1` and `o2` recorded one
  each. All alignment/transformation stages completed and refinement was
  explicitly skipped. Notebook post-processing was canceled after it began
  recursively hashing the large template/reference inventory; the partial
  arm evidence is preserved and is not a biological conclusion.

## Decision Log

- Decision: use `PRISM-prescript` as maintained pipeline and evidence
  authority. Rationale: it contains the current provenance, benchmark, scoring,
  and validation contracts. Date: 2026-09-12.
- Decision: treat `PRISM` as a reviewed feature source, not a merge source.
  Rationale: its modular orchestration refactor was reverted and its current
  worktree lacks equivalent durable validation. Date: 2026-09-12.
- Decision: preserve all current dirty and generated material until it is
  classified. Rationale: benchmark evidence, raw structures, references, and
  temporary artifacts are currently intermingled. Date: 2026-09-12.

## Outcomes & Retrospective

Phase 0 and the repository-boundary decision are complete. The experimental
`PRISM` branch now has restored documentation, explicit transformation
invariants, orientation-safe output names, observable DockQ command metadata,
and 70 passing hermetic tests. The prescript provenance gate has controlled
duplicate-ledger failure, declared expected-inventory checks, nonzero consumer
exit codes, bounded dirty-tree/tool identity, and terminal skipped-stage
records. The dirty prescript worktree is classified in the freeze snapshot;
no generated or reference material was deleted.

The focused prescript contract suite passes 68 tests. The complete prescript
suite passes 344 tests with 6 skips and one SciPy/NumPy compatibility warning.
The runtime manifest validator and a bounded run-identity smoke pass. The
stable TMalign smoke (1658925) and GTalign CPU/GPU smokes (1658928, 1658929)
are plumbing evidence with zero transformed survivors, not biological results.
The optional backend run (1658926) produced successful PyRosetta/FiberDock
records; its first DockQ attempt and explicit replay (1658930) failed nonzero
because the selected `gtalign_env` DockQ extension is incompatible with NumPy
2.4.6. A repository-local compatible scoring environment was then replayed
successfully in job 1658942. Orientation,
MultiProt true-TM, complete external-Rosetta observability, benchmark
authority reconciliation, and release-gate cleanup remain open.

## Context and Orientation

The maintained entry point is `prism.py`. The stable pipeline is NACCESS,
TMalign, and external Rosetta. Optional alignment/refinement/ranking arms must
remain explicit. The legacy MultiProt/FiberDock tree is reference-only.

The current freeze artifacts are:

- `PRISM/tmp/agent/20260912-workflow-freeze-retry1/`
- `PRISM-prescript/tmp/agent/20260912-workflow-freeze-retry1/`

## Plan of Work

First repair the experimental branch and define its adapter contracts. Then
finish the prescript provenance and lifecycle contracts. Only after focused
software validation passes should isolated Slurm smoke runs be launched. Final
benchmark comparisons require frozen inputs, assets, mappings, evaluators, and
complete per-candidate evidence.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM`, document backend input/output/status
   contracts and run isolated backend smoke tests.
2. From `/scratch/rshadi25/GitHub/PRISM-prescript`, complete provenance
   expected-inventory integration and add remaining corruption/identity/hash
   tests.
3. From `/scratch/rshadi25/GitHub/PRISM-prescript`, complete stage and
   candidate terminal-state records for external commands.
4. Build a classification manifest for the dirty worktree before moving or
   deleting any generated or reference file.
5. Run focused software checks, then isolated Slurm stable-baseline and
   backend smoke runs with run-scoped logs and manifests.
6. Complete FiberDock, MultiProt, external-Rosetta, GTalign CPU/GPU, and
   native/o1/o2 validation gates.
7. Resolve benchmark source authority and run only matched causal comparisons.
8. Update memory/planning records and perform the release gate.

## Validation and Acceptance

Acceptance requires clean `git diff --check`, passing focused and full tests,
valid manifests and artifact ledgers, explicit terminal stage/candidate
statuses, no accidental deletions, and no unexplained benchmark rows. A
biological claim additionally requires frozen inputs/assets, hashed model and
native mappings, cross-interface DockQ, and complete failure accounting.

## Idempotence and Recovery

Never use `git clean`, `git reset --hard`, or broad deletion. Every compute
run uses a fresh run root. Failed, canceled, timed-out, and superseded runs
remain preserved. If a step is interrupted, resume from the latest snapshot,
manifest, checkpoint, or handoff rather than regenerating or overwriting
evidence.

## Artifacts and Notes

The freeze records include status, recent history, tracked diff, index diff,
changed paths, compressed untracked paths, benchmark roots, binary hashes, and
cluster state. The active Slurm job recorded during the freeze is `1658869`.

## Interfaces and Dependencies

The work depends on the verified `gtalign_env` Python interpreter, GTalign,
TMalign, NACCESS, external Rosetta, PyRosetta, FiberDock, DockQ, Slurm, and
the prescript provenance/benchmark scripts. Live partition, account, QOS,
memory, and GPU capacity must be checked immediately before any submission.

## 2026-09-12 validation evidence

- Focused prescript validation: 74 passed across CLI, candidate, ranking,
  transformation, provenance, validation-gate, lifecycle, smoke-command, and
  evaluator-integrity tests; compilation and shell syntax passed.
- PRISM validation: 70 passed with compilation and scoped git diff checks.
- Stable baseline: job 1658925, TMalign/NACCESS/external-Rosetta
  no-refinement smoke, isolated root
  tmp/agent/20260912-stable-baseline-smoke-v3/; exit 0, all stages
  terminal, zero accepted candidates.
- GTalign CPU: job 1658928, executable SHA256
  62e253d7d8f372f734be1a07ce044ffb4c8aae07c134ad91032b329bb08737,
  isolated root tmp/agent/20260912-gtalign-cpu-smoke-v2/; exit 0,
  four alignment records with return code/raw hashes, zero candidates.
- GTalign GPU: job 1658929, executable SHA256
  f85d990beaf61e99599c0b65839d00cd29d900202f60e7f671366d4dfc21f1a1,
  isolated root tmp/agent/20260912-gtalign-gpu-smoke/; exit 0, four
  alignment records, zero candidates. CPU/GPU outputs differ and are not
  combined.
- Optional backend: job 1658926, isolated root
  tmp/agent/20260912-optional-backend-smoke/; PyRosetta status success with
  output hash and score 2324.3037505674843; FiberDock parsed
  fiberdock_energies.ref, energy 52119.41, and produced a non-empty
  structure. DockQ failed with an explicit NumPy ABI incompatibility.
- DockQ replay: job 1658930, isolated root
  tmp/agent/20260912-dockq-replay/; machine-readable status records
  `status=failed`, `return_code=1`, and the full ABI/import failure. It is
  preserved as a blocked evaluator attempt, not promoted to a score.
- Repository-local DockQ replay: job 1658942, isolated root
  `tmp/agent/20260912-dockq-repo-env/`; the environment's Python module entry
  point completed both raw DockQ and `score_single_prism_pair.py` with return
  code 0. DockQ `2.1.3` under NumPy `1.26.4` produced
  `GlobalDockQ=0.2116967149685021` for model hash
  `60db5d97de0004184c872d251a5f2c1b6627f153b6fd6906ae8bf98fd365eb3b`
  against native hash
  `a4342d54f0661ddf37538d710bd009af6eb28ae937625959be24502f9b5e3d01`,
  mapping `OA:GF`. The raw JSON hash is
  `985492a3f117371df1c70a352a56d0017dfa0f768074cbafc72badf1c3e4925f`.
  This is one observational evaluator replay, not a quality or ranking claim.
- DockQ pipeline override: job 1658943, isolated root
  `tmp/agent/20260912-dockq-override-smoke/`; the real `gtalign_env` process
  called `src.eval.dockq` with `DOCKQ_PYTHON` set to the repository-local
  interpreter and completed with the same score and model/native hashes.
- DockQ runtime regression: `tests/test_dockq_runtime.py`, together with
  `tests/test_compare.py` and `tests/test_model_output_integrity.py`, passed
  before the raw-JSON CLI addition; the final focused suite includes the new
  CLI test. Both replay launchers pass `bash -n`.
- DockQ CLI JSON retention: job 1658944, isolated root
  `tmp/agent/20260912-dockq-cli-json-smoke/`; the adapter returned zero,
  emitted one raw `model_1glcFG.dockq.json`, and reported DockQ
  `0.2116967149685021`.
- Orientation notebook smoke: job 1658931, isolated root
  tmp/agent/20260912-orientation-study/; native-default, `o1`, and `o2`
  stages completed with `--no-refine`, but the job was canceled during
  post-processing because the provenance asset index was still traversing the
  large template/reference tree. The three arm stage/candidate ledgers remain
  available; notebook comparison tables were not completed.
- Run identity smoke: contract and manifest under
  tmp/agent/20260912-run-identity-smoke/; one template file was hashed,
  contract size was 16604 bytes, and the manifest recorded the attempt ID,
  output root, command, dirty-tree hashes, and explicit tool fingerprints.
- Full prescript software suite: `344 passed, 6 skipped` in 122.94 seconds;
  one warning reports SciPy's NumPy version constraint. The command was
  bounded to 300 seconds and completed with exit 0.
- Global prescript `git diff --check` remains blocked by six pre-existing
  modified benchmark CSV first-line whitespace records and one pre-existing
  trailing-space line in `src/rosetta_refinement.py`; evidence CSVs were not
  rewritten silently.

## Immediate next-action prompt

Use the following bounded prompt for the next execution cycle:

> Continue from `PRISM-prescript` as the maintained source of truth. Do not
> touch `PRISM` wholesale, the legacy tree, canonical benchmark outputs, raw
> references, environments, or Slurm job `1658869`. Keep all work in a fresh
> `tmp/agent/<run-id>/` root and preserve failed/canceled evidence.
>
> 1. Reconcile the DockQ runtime contract. Compare the declared
> `environment.yaml` identity with the verified
> `benchmark/prism_processed/env/prism_score_env/bin/python` identity. The
> latter is currently Python 3.9.23, DockQ 2.1.3, NumPy 1.26.4 and succeeds
> only through `python -m DockQ`; its standalone launcher has a stale shebang.
> Either build a fresh user-managed Python 3.11.13 environment from the
> recipe or formally register the existing prefix. Do not modify
> `gtalign_env`. Record interpreter path, module path, package versions,
> launcher/interpreter hashes, command, and environment status.
>
> 2. Re-run the one-pair evaluator control in a new root using model/native
> hashes and an explicit mapping. Run both raw DockQ JSON and
> `benchmark/scripts/score_single_prism_pair.py`; require nonzero failure on
> either missing JSON, nonzero return code, missing score, invalid raw chain
> contract, or mapping mismatch. Preserve `scored`, `failed`, and
> `not_scoreable` outcomes separately.
>
> 3. Run a small hash-joined cohort, not the full benchmark. Freeze the input
> rows, native source, receptor/ligand roles, model hashes, evaluator mode,
> and output root. Emit one row per candidate with dataset row ID, chain map,
> model/native hashes, DockQ JSON hash, GlobalDockQ, requested cross-interface
> DockQ, grouped iRMSD, command, return code, and terminal status. Audit
> duplicate identities, missing natives, swapped model/native paths, and
> silent row loss before aggregation.
>
> 4. Only after the scorer cohort passes should you resume the scientific
> arms: corrected FiberDock parser replay, calibrated MultiProt true-TM
> orientations, external-Rosetta per-candidate return/output observability,
> and separate GTalign CPU/GPU arms. Do not combine CPU/GPU results while
> hashes, arguments, or parser outputs differ.
>
> 5. Fix the orientation notebook's symlink-safe/bounded asset inventory,
> then rerun native-default, `o1`, and `o2` with identical inputs and
> `--no-refine`. Inspect input parity, alignment return/hash evidence, gate
> ledgers, clashes, and unknown statuses before refinement, native DockQ,
> threshold sweeps, or ranking. Treat zero-pair and survivor counts as
> plumbing observations only.
>
> 6. Update project memory after each milestone, run focused tests first,
> validate through Slurm, preserve exact job IDs and output roots, and do not
> make claims about ranking quality, speedup, causal backend superiority,
> CPU/GPU interchangeability, or complete benchmark validation until the
> final release gates pass.
