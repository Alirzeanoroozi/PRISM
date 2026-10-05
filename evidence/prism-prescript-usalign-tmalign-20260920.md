# PRISM-prescript CPU-array evidence — 2026-09-20

## TMalign audit

- Remote run root: `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919/tmalign`
- Task ledgers: 257 total; 252 `completed`; 5 `running`; no ledger was classified as failed in the audit snapshot.
- A completed-ledger marker audit found 229 rows with the expected 79,420 alignment records, 79,420 successful records, zero failed records, and return code 0. The remaining completed ledgers require summary-level reconciliation before being treated as fully valid.
- Representative completed case `0000`: 79,420/79,420 successful alignments, zero failed records, 10 transformed model pairs, return code 0, empty validation errors, and a recorded summary hash.
- No refinement-named files were found under the TMalign run root. The submitted TMalign command uses `--no-refine`; these outputs are alignment/transformation evidence only, not refined-model outputs.

## KUACC USalign replacement

- `3121375` was confirmed `PENDING` with `QOSMaxCpuPerUserLimit`; `3121376` was confirmed `PENDING` on `afterany:3121375_*`.
- Only those two exact jobs were canceled after scheduler metadata was captured. KUACC TMalign `3121371` and MultiProt `3121373` were not modified.

## VALAR replacement

- Immutable panel: 19,855 templates; SHA-256 `4680d3eda8030861a40373cd193b0e8bef7c21a771a90553c4c49814e964b48d`.
- Manifest: `/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-valar-usalign-2cpu/batch_manifest.json`.
- Namespace: `/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-valar-usalign-2cpu/`.
- Configuration: VALAR `kutem`, 2 CPUs, 4 GiB, array throttle `%24`, 257 cases, 79,420 records per case, rank/PRODIGY/refinement disabled.
- First submission `1676217`/`1676218` failed before case execution because the framework `src/` path was absent from `PYTHONPATH`; all 257 case ledgers remained `not_started` and no attempt directories were created.
- Corrected wrappers export the framework `src/` path. The corrected array is `1676493`, with dependent aggregator `1676494`. Early smoke state: 24 running, 233 pending, 24 ledgers `running`, zero failed.

## Latest transformation count

- A ledger-selected attempt-root audit now finds 253 completed cases with
  34,230 transformed PDB files: 17,115 `*_L.pdb` files and 17,115
  `*_R.pdb` files, i.e. 17,115 finalized model pairs.
- Four cases are still running and currently contribute 198 additional
  in-progress PDB files. Therefore the live current-attempt filesystem total
  is 34,428 PDB files, while the reconciled completed total is 34,230.

## TMalign downstream smoke test — 2026-09-21

- Isolated source: case `rigid_1ahw_000`, template `2ec9TU`, orientation `o1`.
  Only the two transformed PDB parts, matching native complex, and case
  manifest were copied into
  `/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260921-tmalign-pyrosetta-dockq-smoke/`.
- PyRosetta refinement passed in `gtalign_env` in 27.1 seconds and produced a
  refined PDB with recorded metadata and output SHA-256
  `3dfefdd4281d6f3764d1409c894d7d47b231a5d5036d53580ffc97bb2dba5cad`.
- Direct DockQ/iRMSD on the original transformed pair passed: one scored row,
  GlobalDockQ `0.296464`, CAPRI `Acceptable`.
- DockQ/iRMSD on the PyRosetta-refined output passed: one scored row,
  GlobalDockQ `0.0395533`, CAPRI `Incorrect`.
- The first DockQ attempt exposed interpreter selection of NumPy 2.4.6 and
  was recorded as `score_failed`; rerunning with the repository's explicit
  NumPy-1-compatible `DOCKQ_BIN` launcher passed. This is a smoke-test
  observation only, not a conclusion about full-run refinement quality.

## USalign runtime/parser diagnosis — 2026-09-21

- Live run root:
  `/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-valar-usalign-2cpu/`.
- Array `1676493`: 24 tasks running on `rk01`, 233 pending with
  `JobArrayTaskLimit`; no started ledger is failed. Aggregator `1676494` is
  dependency-pending. No case summary or transformed model is present yet.
- Case `0000` was actively writing alignment records while inspected:
  `50,056/79,420` records in one snapshot, all sampled records with
  `status=success`, nonzero `match_count`, and `tm_score=0.0`. The sampled
  stage ledger has completed input and surface extraction and is still in
  alignment. Cases `0001`, `0002`, and `0023` showed approximately 50,549,
  51,808, and 74,500 records respectively.
- The runner maps `--aligner usalign` to the TMalign-compatible PRISM path,
  sets `PRISM_USALIGN_PATH` to
  `/home/rshadi25/.conda/envs/gtalign_env/bin/USalign`, and sets workers to
  `min(4, SLURM_CPUS_PER_TASK)`. With 2 CPUs, each case has two external
  processes. `run_cpu_pipeline.py` writes `run_summary.json` only after the
  whole case subprocess returns; transformation is later, so partial
  alignment JSON is not surfaced as final results.
- Slurm batch-step telemetry for actual task `1676495` showed
  `AveCPU=03:28:51`, `MaxRSS=124000K` at about 11:25 wall time, with 2 CPUs
  requested. This is evidence of low average CPU occupancy rather than a
  saturated CPU-bound task.
- Completed TMalign reference case `0235` used the same 79,420-record panel:
  3,338.0 s alignment, 3,754.0 s total, four workers, and 79,420 successful
  records. Direct same-input checks measured TMalign/USalign wall times of
  0.011/0.092 s on a small pair and 0.128/0.232 s on a larger pair.
- Exact output comparison: TMalign prints `TM-score= ... normalized by
  length of Chain_1/Chain_2`; USalign prints the corresponding
  `Structure_1/Structure_2` labels. `src/alignment.py` checks only for
  `Chain_1` and `Chain_2`, leaving USalign scores at zero while retaining
  mappings. This is a confirmed parser correctness defect. The same parser
  hard-codes the record provenance field to `TMalign`.
- Interpretation: scheduling is progressing; the delay is caused by the
  very large per-case process/record workload, slower USalign calls under the
  two-worker setting, shared-filesystem/process overhead, and final-result
  publication only at case completion. The zero-score parser defect is
  independent and makes any eventual USalign ranking/scientific result
  invalid until corrected.

## Continuation decision and manuscript comparability — 2026-09-21

- The current run should not be resumed or duplicated unchanged for
  scientific production. It is operationally progressing, but its score
  fields are invalid and no case has reached a validated final result.
- The installed launcher resolves to
  `/scratch/rshadi25/GitHub/Template-based-structure-aligners/old_pipeline/usalign/USalign`,
  reports Version 20241108, and the runner invokes pairwise default USalign
  with only `-m matrix`. It does not pass `-fast`, `-dir1`, `-dir2`, or
  `-mol prot`.
- Same-input benchmark: default versus `-fast` USalign was 0.240 versus
  0.132 wall seconds on a larger staged interface pair. The `-fast` output
  changed aligned length from 109 to 93 and changed the reported scores, so
  it cannot be adopted without a quality comparison.
- The present experiment is not evidence against the manuscript result that
  USalign can be faster than MultiProt. That result may depend on exact
  binary/version, `-fast` or directory-mode batching, input size and
  structure composition, template scope, worker count, and hardware. Here,
  79,420 pairwise launches per case and shared-filesystem temporary output
  create a different workload.
- Recommended controlled pilot before any full rerun: repair the USalign
  score-label/provenance parser; run default and `-fast` on identical sampled
  pairs used by MultiProt; record per-call wall time, CPU time, failure rate,
  match coverage, transformed candidates, and DockQ/iRMSD; then select the
  mode and resources from measured evidence. Current jobs were not modified.

## USalign ETA estimate — 2026-09-21

- Latest live sample: 24 active cases, approximately 11.97 h elapsed. Six
  high-progress cases had 0.1–0.9 h of alignment remaining; the other 18 had
  roughly 5.7–7.0 h remaining. The active wave should therefore finish in
  approximately 7–7.5 h including a small transformation/finalization tail,
  if rates remain stable.
- Aggregate progress was approximately 1,383,600/20,410,940 records at
  32.09 records/s. Extrapolating the remaining 19,027,340 records through
  the 24-task array throttle gives approximately 6.86 days for the complete
  array, before a modest post-alignment overhead. This estimate assumes
  later cases have comparable difficulty and is not a scheduler guarantee.

## Cancellation and cleanup — 2026-09-21

- Canceled only USalign-owned VALAR jobs `1676493`, task `1680421`, and
  aggregator `1676494`. No KUACC USalign job was visible, and no unrelated
  job was modified.
- Started removal of the generated run root
  `/scratch/rshadi25/GitHub/PRISM/tmp/agent/20260920-valar-usalign-2cpu/`.
  The removal process remains in uninterruptible shared-filesystem I/O, so
  the root is still present and cleanup is explicitly incomplete.
- Preserved durable wrappers and diagnostic documentation. No replacement
  workload was submitted.
- Post-cancellation contract checks passed: wrapper shell syntax, Python
  compilation, and USalign Version 20241108 help/version output.

## Manual USalign validity pilot — 2026-09-21

- Inputs were three single-chain interface pairs from
  `/scratch/rshadi25/GitHub/PRISM-prescript/templates_test/interfaces/`:
  `3hpgAF_A`/`3hpgAF_F`, `3izlAB_A`/`3izlAB_B`, and `1xhzCD_C`/`1xhzCD_D`.
  This matches the pairwise single-chain shape used by the alignment adapter.
- All TMalign, USalign default, and USalign `-fast` calls returned 0. Wall
  times in seconds were `(0.024, 0.051, 0.046)`, `(0.021, 0.047, 0.043)`,
  and `(0.017, 0.043, 0.042)` respectively.
- The same three pairs ran through MultiProt with return code 0, `2_sol.res`
  present, and Largest Solution values 37, 16, and 15; wall times were
  0.048, 0.040, and 0.037 seconds. This bounded pilot does not reproduce a
  large USalign speed advantage over MultiProt.
- The corrected parser regression test passed 8/8. USalign Structure-label
  scores are now recognized and the runner passes explicit `USalign`
  provenance. The compatibility `tm_score` remains the existing maximum of
  both normalized scores; both raw scores should be retained/compared if the
  project selects strict reference-normalized scoring.
- Manual output roots:
  `/scratch/tmp/usalign-manual-pairs.ndeW2Y/` and
  `/scratch/tmp/multiprot-manual.Z2CEOz/`.
