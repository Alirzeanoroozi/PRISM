# PRISM aligner comparison evidence package

This run package is the compact evidence surface for the staged comparison of
TMalign, MultiProt, USalign, and GTalign GPU. It is intentionally split into
panel lanes:

- historical/working lane: 19,855 templates and 257 BM5.5 cases;
- exact-panel lane: 19,948 checked, 19,062 calculated, and 19,058 materialized
  templates, with four recorded materialization exclusions.

The historical TMalign/MultiProt ledgers and GTalign transformed scores are
reusable but are not silently pooled with exact-panel alignment-only evidence.
USalign production is resumable through Slurm jobs recorded in
`usalign_production_19855/production_manifest.json`; scheduler state is not
accepted as scientific completion. The stage ledger is the authoritative
reconciliation index.

Every downstream result must retain explicit status (`scored`,
`valid_unscored`, `scored_cross_only`, `not_scoreable`, or `score_failed`),
GlobalDockQ separately from requested cross-interface scores, both USalign TM
normalizations, hashes, commands, and panel provenance. Disposable raw and
structure artifacts are eligible for deletion only after their compact result
and validation manifest are durable; unresolved score failures remain
auditable.

The final package will add case-level/candidate-level quality tables,
candidate overlap, ranking/top-k, paired refinement deltas, timing/resource
tables, and a scientific interpretation after the active USalign and common
refinement stages validate. Until then, no overall best-aligner claim is
authorized.

Current downstream execution:

- GTalign common refinement uses 14,470 materialized models in array `1709167`;
  read-only corrected aggregation is dependency job `1709178`.
- USalign historical-panel production is array `1708992`, followed by compact
  replay `1709007` and transformed DockQ `1709046`.
- `eda_partial_historical_19855/` is a compact, explicitly incomplete reducer
  output for TMalign/MultiProt/GTalign; it is not the final matched analysis.
