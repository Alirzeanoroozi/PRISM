# Restore Canonical BM5.5 Scoring for Current Ranking Outputs

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` current.

## Purpose / Big Picture

Produce reproducible BM5.5 quality labels for the completed current pipeline variants so the existing ranking option can be evaluated and trained. Preserve every stable pipeline and prior failed scoring artifact. Only derived staging manifests, canonical native assemblies, strict DockQ JSON, and new aggregate tables may be added.

## Progress

- [x] Audit the failed variant scorer and separate job completion from scientific score validity.
- [x] Identify the canonical strict bijective scorer already covered by focused tests.
- [x] Correct current-model discovery, durable row identity, and query-partner chain grouping in staging.
- [x] Assemble native complexes from curated BM5.5 archive role files with hashes and explicit failures.
- [x] Validate representative monomer/monomer and multichain/monomer cases through DockQ 2.1.3.
- [x] Submit the smallest complete scoring set through Slurm and audit its outputs.
- [ ] Scale canonical scoring to the full current variant (job 1393041 running; audit 1393064 dependent).
- [x] Update project memory and Graphify feedback with verified commands, artifacts, failures, and chronology.
- [x] Harden ranking-label attachment around `dataset_row_id`, model SHA256, source-gate eligibility, and requested cross-interface DockQ.
- [x] Reconstruct and audit all 235 batch-1 ranking feature rows from retained GTalign JSON without rerunning alignment.
- [ ] Run the full provenance-safe feature/label reconstruction and deterministic baseline evaluation (job 1393128, dependent on accepted score audit 1393064).
- [x] Reconstruct batch-1 ranking features, rank per benchmark row, and record per-row baseline/oracle metrics.

## Surprises & Discoveries

- Completed Slurm generation jobs do not imply valid benchmark scores; the ad hoc scorer produced failure rows rather than quality labels.
- `score_bm55_variant.py` inferred the native from a template token, discarded `dataset_row_id`, guessed chain order, and parsed human-readable DockQ output. It is retained only as failed-attempt evidence.
- Current external-Rosetta outputs end in `_rosetta_0001_0001.pdb`, while the staging adapter discovered only `_rosetta_0001.pdb`.
- Refinement sizes model partner groups from query receptor/ligand selectors. The staging adapter instead split observed chains using template partner counts, corrupting multichain mappings.
- The existing graph uses pre-path-qualified node IDs and traverses vendored test environments. Graph chronology is supporting evidence only; dated manifests, source locations, Slurm accounting, and current memory are authoritative.
- Batch 1 contains 636 Rosetta artifacts but only 235 distinct final poses after selecting the highest refinement suffix per pose.
- The first internally parallel Slurm run, job `1392965`, retained 235 explicit failures because the launcher resolved the DockQ virtual-environment Python symlink into its base interpreter.
- Preserving `/scratch/tmp/prism-dockq-env/bin/python` without `Path.resolve()` fixed the clean-batch environment. Job `1392974` scored all 235 final batch-1 poses in 2m54s with DockQ 2.1.3, 938 interface/global records, and zero model/native/raw-JSON hash mismatches.
- GlobalDockQ can be dominated by native receptor-internal interfaces in multichain cases. Ranking labels must use the requested receptor-ligand cross-interface metrics, reported separately from GlobalDockQ.
- The July generation runs did not retain candidate-audit JSONL, but they did retain authoritative GTalign JSON in each batch's `processed/alignment_gtalign/<run-id>/` directory. Learned-reranker features can be reconstructed without rerunning alignment if the adapter preserves batch/run, orientation, template-chain, model-hash, alignment-JSON hash, and raw-output-hash provenance.
- Accepted full preparation job 1393026 retained 6,539 model records: 5,828 valid strict-clean poses, 359 explicit strict-clean model-contract failures, 282 staged audit-only poses, and 70 audit-only model-contract failures. All 156 represented strict-clean natives assembled; eight native failures are audit-only.
- Every batch contributing staged models has exactly one retained GTalign run directory. The zero-JSON directories for batches 0005, 0006, 0007, and 0020 contribute no staged models and are not a reconstruction gap.
- Batch-1 feature reconstruction recovered and labeled 235/235 candidates. The ranking audit found zero model-hash, alignment-JSON-hash, raw-output-hash, aligner/status, feature, or label-contract failures; all labels use `dockq_cross_mean`.
- Batch-1 evaluation identified a ranking bug rather than a scoring bug: production audits never supplied template sizes, so coverage was absent and the baseline invented `mean_match_count / 50`. Passing real template sizes for future audits and using TM-score alone when historical coverage is unavailable changed batch-1 native-like top-1 from 1/10 to 3/10 (all three oracle-positive rows). This remains diagnostic pending the full independent-row evaluation.
- Offline ranking and reranker evaluation originally treated unrelated native
  complexes as one global list. Both now evaluate within `dataset_row_id` or
  native-complex groups. Batch 1 reconstructed and labeled all 235 candidates;
  3/10 rows had an oracle native-like candidate and baseline top-1 succeeded
  for 1/10.

## Decision Log

- Decision: use `dataset_row_id` as benchmark identity and preserve raw selectors and source row.
  Rationale: PDB/template names are not unique row identities and aliases exist.
  Date/Author: 2026-07-26 / Codex.
- Decision: use `score_bijective_benchmark_models.py` with complete within-partner bijections, raw DockQ JSON, hashes, and explicit failure rows.
  Rationale: it implements the verified multichain contract; the ad hoc variant scorer does not.
  Date/Author: 2026-07-26 / Codex.
- Decision: native structures must be derived from curated BM5.5 role files, not downloaded full-PDB substitutes.
  Rationale: archive role files are the verified benchmark truth and avoid extra-chain ambiguity.
  Date/Author: 2026-07-26 / Codex.
- Decision: stable pipeline outputs and failed scoring artifacts are read-only evidence.
  Rationale: repairs belong in isolated `tmp/agent/20260726-bm55-canonical-scoring/` outputs.
  Date/Author: 2026-07-26 / Codex.
- Decision: rank and evaluate independently within each durable benchmark row.
  Rationale: top-K is a per-query choice; a global list across unrelated
  complexes is biologically and operationally meaningless.
  Date/Author: 2026-07-26 / Codex.

## Validation and Acceptance

- Every staged row retains `dataset_row_id`, source row, raw selectors, model hash, and exact model partner groups.
- Model discovery accounts for both current external-Rosetta and PyRosetta layouts without duplicate staging.
- Observed model chains exactly equal expected query partner counts; extra or missing chains are explicit failures.
- Canonical native files are tied to curated archive members and hashes.
- DockQ runs with `/scratch/tmp/prism-dockq-env/bin/python`, records version 2.1.3, complete mapping, raw JSON, model/native hashes, GlobalDockQ, requested cross-interface scores, and grouped iRMSD.
- Representative focused tests and a Slurm scoring smoke pass before wider submission.

## Idempotence and Recovery

All generated artifacts use a fresh directory under `tmp/agent/20260726-bm55-canonical-scoring/`. Staging refuses a non-empty output root. Existing benchmark results and stable pipeline files are never overwritten or deleted.

## Outcomes & Retrospective

Pending corrected staging, canonical native assembly, and a valid strict-scoring smoke.

Batch-1 acceptance is complete. The canonical path staged 235 unique final poses across 10 strict-clean rows, assembled 10 row-specific natives from curated bound roles, and scored every pose successfully in job `1392974`. Median requested cross-interface DockQ is 0.005077 (range 0.002175-0.859206); median grouped iRMSD is 28.162 A (range 0.874-48.876). The next action is a full-variant preparation audit followed by full scoring only if identity, chain, and source-gate counts remain consistent.

Full preparation is now accepted. Full scorer 1393041 is running on 5,828 eligible strict-clean poses, with dependent audit 1393064 enforcing 6,539 total model records, no score failures, complete mapping/interface records, and hash integrity. The combined regression suite passes 44 tests with one expected skip.

The retained-alignment adapter is now fail-closed: it rejects ambiguous run paths, unsuccessful/non-GTalign payloads, missing raw-output hashes, and model-hash mismatches. Batch-1 outputs are `scored-batch-0001-v2/ranking-{candidates,labeled}.csv` plus `ranking-audit.json`. Full reconstruction/evaluation job 1393128 is chained with `afterok:1393064`, so an unaccepted score table cannot enter ranking analysis.
