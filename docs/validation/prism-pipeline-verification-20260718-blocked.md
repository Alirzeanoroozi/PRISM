# PRISM pipeline verification status (blocked, 2026-07-18)

This is a validation-status record, not a benchmark-quality report. No PDF files
were read. Confirmatory claims are intentionally withheld.

| Gate | Status | Evidence / reason |
|---|---|---|
| 0 | complete | Baseline manifest: `tmp/agent/20260718-prism-verification/baseline/` |
| 1–2 | complete | Cancellation-aware stages and refined-pose integrity; AI `1368213`, COSBI `1368214` |
| 3 | complete (contract) | Canonical scorer retry `1368545`: 30 model rows, 22 scored, 8 explicit DockQ failures, 58 interface rows; all failure rows retain native/model hashes |
| 4 | complete (provenance) | AI `1368528` / COSBI `1368530`: 1,270,720 similarity rows, zero missing interface assets, byte-identical TSV hashes |
| 5 | unresolved | Modern/derived filter parity: 19,005 disagreements and 850 missing modern profiles |
| 6 | complete (alignment scope) | 25-, 946-, and 19,855-template GTAlign/TM-align diagnostics reconcile on AI/COSBI |
| 7 | partial | Relaxed positive canary passes; normal-threshold AI/COSBI canaries `1368555`, `1368556`, `1368560`, `1368561`, `1368567`, `1368568` are explicit `completed_no_predictions` |
| 8–9 | blocked | Frozen source policy has `status=blocked_source_authority` and `confirmatory_run_authorized=false` |
| 10 | partial | Confirmatory preflight exists and correctly returns `blocked`/exit 2; no claim ledger is published |

## Reproducibility anchors

- Frozen template list SHA-256: `4680d3eda8030861a40373cd193b0e8bef7c21a771a90553c4c49814e964b48d`
- AI/COSBI similarity TSV SHA-256: `f6ab1745944cfdc506fef9ee657a2425d744d7a48a60eb0829e64bed23cf3063`
- AI/COSBI eligible-list TSV SHA-256: `097ea0b78af75e87875d3c3d6952bfefdfe0476bea7bcdc2d88ec72f5258e6b`
- Confirmatory preflight: `tmp/agent/20260718-prism-verification/confirmatory-preflight/source-policy-blocked.json`

The next authorized action is to resolve the source-authority decision and the
modern-versus-derived filter-asset semantics. Until then, Gate 8 and Gate 9
must not be submitted.
