# PRISM-prescript DockQ contract-correction report

## Verdict

`VALIDATED` for the parser/scorer software contract; `NOT PROMOTION-READY`
for new scientific benchmark conclusions until corrected BM55 scoring and EDA
are generated.

## Findings

- Complete mappings use raw DockQ 2.1.3 `GlobalDockQ`.
- `best_dockq` remains a named diagnostic sum and is never used as an
  unqualified global score.
- A parseable JSON document without complete `GlobalDockQ` is
  `valid_unscored`, with no interface-score fallback.
- Pairwise recovery after a complete mapping failure uses
  `score_scope=requested_cross_interfaces_only` and
  `score_status=scored_cross_only`; `GlobalDockQ` stays unavailable.
- Explicit no-align checks validate chain bijection, residue numbering,
  insertion codes, and residue identity before invoking DockQ.
- Raw JSON path/hash, exact argv, mapping, no-align mode, and CPU count are
  retained.

## Validation

- PRISM-prescript focused suite: 57 passed.
- Changed Python files compiled and benchmark Slurm scripts passed `bash -n`.
- Scoped whitespace check passed.
- Pinned DockQ help exposed `--json`, `--mapping`, `--no_align`, and
  `--n_cpu`.
- Independent final review: `PASS WITH CAVEATS`; no defect was found in the
  checked scope/status claims.

## Boundaries

The retained raw benchmark outputs and historical failed-attempt scripts were
not regenerated or overwritten. The canonical staged scorer is the required
route for the next corrected BM55 run; EDA must filter by status and scope.
