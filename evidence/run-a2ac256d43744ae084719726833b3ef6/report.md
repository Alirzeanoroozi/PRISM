# PRISM-prescript alignment-runtime repair — Run report (CORRECTED)

Run: `run-a2ac256d43744ae084719726833b3ef6`
Workstream: `prism-prescript-alignment-runtime-repair`
Date: 2026-09-08 (continuation-03 correction pass)
Execution workspace: `/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-a2ac256d43744ae084719726833b3ef6` (branch `valar/run-a2ac256d43744ae084719726833b3ef6`)
Copilot session: `78fae6af-9164-5f65-b7e5-bd9bc8792dfd`
Hosting Slurm job: 1656898 (RUNNING, ai26) — code evidence below was regenerated inside this job

## Correction notice (continuation-03)

The previous closeout incorrectly stated "USalign executable is not available in this environment". **That statement is wrong for the environment and correct only for the adapter's default contract path.** Corrected facts (OBSERVATION, verified 2026-09-08 inside job 1656898):
- `/home/rshadi25/.conda/envs/gtalign_env/bin/USalign` exists, is executable (136-byte bash wrapper: `exec /scratch/rshadi25/GitHub/Template-based-structure-aligners/old_pipeline/usalign/USalign "$@"`).
- The wrapped binary exists, is **ELF 64-bit x86-64**, 861,320 bytes, `file`-verified.
- `USalign -h` → "US-align (Version 20241108) ... Usage: USalign PDB1.pdb PDB2.pdb [Options]", exit 0 (help flag only; no structural alignment run on the login node).
- `check_usalign_runtime()` (the adapter's own function) → `can_run: true` for both the wrapper and the real binary path; `can_run: false` / `executable_not_found` for the default contract path `external_tools/USalign`, which genuinely does not exist in the worktree. So the adapter's fail-closed default is correct behavior; a USalign run requires the explicit contract: `--usalign-path /home/rshadi25/.conda/envs/gtalign_env/bin/USalign` (or `PRISM_USALIGN`).
- Seccomp mode inside job 1656898 (ai26): 0.

## Bounded objective (unchanged)

1. `prism.py` rejected `--aligner usalign`; USalign path disconnected (the protein-DNA runner already invokes `prism.py --aligner usalign`).
2. MultiProt blocked on the observed Seccomp/32-bit runtime combination.

Constraints honored: canonical checkout untouched; no USalign/PRODIGY download or build; no hosted models; no Seccomp bypass; stable TMalign/NACCESS/external-Rosetta defaults unchanged.

## OBSERVATION (directly observed, reproducible)

- Root cause of the CLI failure: `src/alignment_multiprot.py` (commit `56aae1fd860`) imports `_has_ca_atoms`/`_write_empty_alignment` from `src/alignment`; no file on any branch defines them (grep + `git log --all -S`). `import prism` crashed before argparse, which subsumed the `--aligner usalign` rejection.
- USalign availability: see Correction notice above. `external_tools/` in the worktree (committed on this branch: `git ls-tree HEAD -- external_tools/`) contains only TMalign, TMalign.cpp, naccess — no USalign, no MultiProt.
- MultiProt bounded search (2026-09-08): worktree `external_tools/` negative; `working_version/` is **not committed on this branch** (`git ls-tree HEAD` shows no such tree) and is absent from the worktree; `/scratch/tmp` contains no MultiProt binary (only a historical figure `suppfig-smoke/s3f_multiprot_complete_tm_score_vs_rmsd.png` and this run's own synthetic pytest ELF fixtures); canonical-checkout access is denied from this sandbox. **The bounded search is negative; this does NOT establish that MultiProt is globally absent** — project memory documents a prior unrestricted Seccomp:0 validation (job 1416942, 47/56 pairs) with the default `external_tools/multiprot.Linux`, so the binary existed in some context previously. Where it currently lives in the canonical dirty checkout: **UNKNOWN** (access denied).
- `check_multiprot_runtime()` results (job 1656898, seccomp mode 0): documented default `external_tools/multiprot.Linux` → `exists: false, can_run: false, reason: executable_not_found`; legacy path `working_version/Multiprot-new/prism-fiberdock-cli/` → same (path absent in worktree). (`evidence/multiprot-runtime-checks.json`)
- Real-tool smoke: **not run**. Exact reason: no safe existing small local PDB pair was available to this sandbox — the worktree contains no `.pdb` files at all (fresh worktree, no data dirs); canonical-checkout reads (templates_test/, new_template/, processed/) are permission-denied; `/scratch/tmp` PDB candidates are dangling symlinks/test fixtures (one points to a placeholder path). Per the bounded instruction, the `-h`/runtime evidence is retained instead.
- Slurm discrepancy: manifest history shows job 1656884 RUNNING 18:10:50Z → **FAILED 18:34:51Z**, then retry 1656898 RUNNING (this session). `logs/slurm.log` (14 lines) contains only framework preflight lines (RuntimeWarning + profile/preflight) — **no exit code or stderr for the failed job is captured in run logs**, and `scontrol show job 1656884` now returns "Invalid job id specified" (record purged), so the failed job's exit code is **UNKNOWN**.
- Graphify in-memory trace (cache_root=None, parallel=False; no write to project `graphify-out/`): 6 bounded files → **76 nodes / 137 edges**; adapter containment edges (`alignment_usalign.py -contains-> align_usalign()`, same for gtalign/multiprot/tmalign) and call edges (`_align_one() -> check_multiprot_runtime()/_parse_multiprot_solution()/_write_empty_alignment()`, `align() -> parse_tmalign()`) extracted. Static limit: `main()` dispatches adapters through `run_stage(..., lambda: align_*)`, so the lambda-indirected call edges do not appear in the graph — dispatch correctness is covered by the explicit-branch test instead. Known graphify limitation observed: one cross-module false-positive `calls` edge (name-based resolution). Artifact: `evidence/graphify-command-trace.json`. (Continuation-03 note: `cache_root=None` still writes graphify's default AST cache to `./graphify-out/` in the process CWD; the trace is re-run with CWD set to `validation/graphify-scratch/` under this run's evidence namespace so no project directory is written — result identical at 76/137. See decisions.jsonl entry 012.)

## INFERENCE

- The MultiProt adapter was committed broken; `import prism` has been dead on this branch since `56aae1fd860`.
- The failed Slurm job 1656884 did not alter code state: the worktree diff and the full 34/34 test result were regenerated inside the current job 1656898, and the retry chain preserved the same worktree/session. (The failure's cause itself remains UNKNOWN — no exit code captured.)

## HYPOTHESIS (open)

- `parse_usalign_output()` conforms to the **real** USalign 20241108 stdout layout. It is written against the documented layout and passes synthetic tests only; a real alignment was not run (no safe PDB pair), so parser conformance is **not validated**.
- Whether the legacy MultiProt binary still exists somewhere in the canonical dirty checkout or elsewhere on the cluster: UNKNOWN.

## IMPLEMENTED (isolated worktree only, unchanged by continuation-03)

- `src/alignment.py` (+42): `_has_ca_atoms()`, `_write_empty_alignment(status, extra)` — root-cause fix.
- `src/alignment_usalign.py` (new): adapter with explicit executable contract (`--usalign-path` > `PRISM_USALIGN` > `external_tools/USalign`; cwd → repo root → `shutil.which`); fail-closed `RuntimeError` with searched path + no-fallback statement; per-pair `status="failed"` records with diagnostics/reproduce command; `parse_usalign_output()`; `_build_match_dict()` in the PRISM `chain.AA.resnum` convention; `tm_score_contract="standard_length_normalized"`.
- `prism.py` (+21/-2): `usalign` in `--aligner` choices; explicit dispatch branch **before** the GTalign catch-all `else`; `--usalign-path`/`--usalign-workers`.
- `src/alignment_multiprot.py` (+135/-19): `_read_seccomp_mode()`, `_elf_class()`, `check_multiprot_runtime()`; `_align_one()` fail-closed records `skipped_runtime_incompatible` (ELF32 + seccomp ≥ 2, exit-159 signature, reproduce command) / `skipped_unavailable`. `PRISM_MULTIPROT_FORCE` not set and not required.
- `tests/conftest.py` (new; find_spec-guarded stubs for 4 uncommitted optional modules — no-op in canonical checkout), `tests/test_alignment_usalign.py` (14 tests), `tests/test_alignment_multiprot_runtime.py` (17 tests).

## TESTED (fresh evidence, job 1656898, 2026-09-08)

- `python3 -m pytest tests/ -q` → **34 passed in 1.78s** (`evidence/pytest-focused-2026-09-08-corrected.txt`; identical to the pre-correction run).
- `python3 -m py_compile` over all 9 touched/related files → OK (`evidence/py-compile-corrected.txt`).
- `check_usalign_runtime()` on wrapper + real binary + default paths → `evidence/usalign-runtime-validated.json`.
- `check_multiprot_runtime()` on documented default + legacy path → `evidence/multiprot-runtime-checks.json`.
- `USalign -h` → Version 20241108, exit 0.
- Graphify in-memory extract/build → `evidence/graphify-command-trace.json` (76 nodes / 137 edges).

## VALIDATED (what the evidence actually supports)

- Import crash fixed; CLI parses and dispatches `--aligner usalign` to the USalign adapter (explicit branch, not GTalign catch-all).
- USalign executable availability **via the explicit contract path** (`/home/rshadi25/.conda/envs/gtalign_env/bin/USalign`): validated by existence/exec-bit/`file`/`-h` and by the adapter's own `check_usalign_runtime()` (can_run=true). The default path correctly fails closed.
- Fail-closed behavior for both adapters (unit + stub-integration tests).
- MultiProt runtime incompatibility represented by explicit, reproducible, fail-closed records without `PRISM_MULTIPROT_FORCE`.
- Stable defaults unchanged in the diff: aligner `tmalign`, surface `naccess`, refiner `external_rosetta`.

## NOT VALIDATED (stated, not claimed)

- **No live USalign alignment was run** (no safe existing PDB pair in this sandbox) → parser conformance with real USalign stdout remains a hypothesis. Next action: with a small real PDB pair on the target node, run one `USalign PDB1 PDB2` via `--usalign-path .../bin/USalign` and confirm `parse_usalign_output()` recovers TM-score/R/t/pairs.
- **No live MultiProt execution**; its current location is UNKNOWN to this sandbox; the Seccomp/32-bit block is represented by diagnostics + synthetic tests, not a live reproduction.
- Slurm job 1656884 failure cause: UNKNOWN (record purged, no exit code in logs).

## REVIEWED

- Diff re-reviewed in continuation-02/03 for imports/paths/JSON schema/error handling; no new issues.
- All writes confined to the isolated worktree and this run directory; no merge/push/reset/delete; canonical checkout untouched.

## Evidence index

- `evidence/pytest-focused-2026-09-08.txt`, `evidence/pytest-focused-2026-09-08-corrected.txt`
- `evidence/py-compile.txt`, `evidence/py-compile-corrected.txt`
- `evidence/usalign-runtime-validated.json` (continuation-03)
- `evidence/multiprot-runtime-checks.json` (continuation-03)
- `evidence/graphify-command-trace.json` (continuation-03)
- `evidence/runtime-environment.json` (continuation-02; note: its "USalign absent" line refers to the default contract path only)
- `evidence/diff-stat.txt`, `evidence/git-status.txt`, `evidence/diff-tracked.txt`

## Next action (single, bounded)

On the target node with a small real PDB pair (e.g., from `templates_test/` in the canonical checkout):
`cd <worktree> && /home/rshadi25/.conda/envs/gtalign_env/bin/USalign pair1.pdb pair2.pdb > out.txt` then `python3 -c "import sys; sys.path.insert(0,'.'); from src.alignment_usalign import parse_usalign_output; print(parse_usalign_output(open('out.txt').read())[:2])"` — record conformance in run provenance before any production `--aligner usalign` run.

## ORCHESTRATOR RECONCILIATION — 2026-09-08T19:08:36Z

### OBSERVATION

Slurm job `1656898` is terminal `FAILED`, exit code `1:0`, after 27:28. The
framework manifest is reconciled to `current_state=FAILED` and
`validation_status=NOT_STARTED`; the reservation is released.

### EVIDENCE

- `evidence/copilot-objective-activation.json` is valid and records the exact
  persisted prompt SHA-256, session, run, objective artifact, and verified
  dispatch.
- `logs/copilot-output.jsonl` contains one prompt-matching dispatch,
  `session.task_complete`, and a `result` with process `exitCode=0`.
- `logs/slurm.log` reports the objective-evidence failure only because the
  captured `artifacts/qwen-worker.sh` uses `if ! python3 ...; then` while
  placing success actions in `then` and failure actions in `else`. A verifier
  exit 0 is therefore inverted into the shell failure path and exit 1.
- The isolated worktree status contains exactly the intended seven worker
  files; the canonical checkout remains at 373 pre-existing status entries.

### INFERENCE

The terminal Slurm failure is a framework-wrapper bookkeeping defect. It does
not establish that the Copilot provider failed or that the isolated code
changes are correct. The worker’s retained test evidence is still evidence of
`TESTED` behavior only.

### UNKNOWN

The isolated implementation has not been independently accepted or promoted.
Real USalign parser conformance, live MultiProt execution under the relevant
restricted runtime, and production-pipeline behavior remain unknown.

### STATUS BOUNDARY

- `IMPLEMENTED`: isolated worktree changes only.
- `TESTED`: retained 34-test and compile results.
- `VALIDATED`: explicit USalign availability, fail-closed paths, and static
  dispatch claims only, as previously recorded.
- `REVIEWED`: orchestrator reconciliation and canonical-state audit completed.
