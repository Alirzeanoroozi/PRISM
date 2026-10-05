# Step 5 current validation map

## Maintained runtime path

```text
inputs.csv
  -> prism.py
       -> surface extraction (NACCESS by default; FreeSASA option)
       -> aligner provider branch
            -> src.alignment              (TMalign)
            -> src.alignment_gtalign      (GTalign)
            -> src.alignment_multiprot     (MultiProt)
            -> isolated candidate only: alignment_usalign (USalign)
       -> src.transformation.transformer
            -> transformed PDB candidates
            -> candidate audit records when audit output is enabled
       -> candidate_selector (optional)
            -> baseline or isolated PRODIGY candidate ranking
       -> selected candidates
       -> refinement (external Rosetta default; FiberDock/PyRosetta alternatives)
       -> optional DockQ comparison
```

Direct source inspection is authoritative for runtime edges. `prism.py` imports
the three maintained providers above and calls the transformer directly; it
does not import or instantiate `src/structural_aligner.py`. The generic
`StructuralAligner` file is therefore diagnostic/disconnected code, not a
validated maintained entry point. The isolated USalign worktree adds an
opt-in `prism.py` branch and `src/alignment_usalign.py`; it has not been
promoted to the canonical checkout.

Graphify was run on a temporary copy of the canonical `src/` tree because the
canonical Graphify output was absent. It reported 44 code files, 501 nodes,
926 edges, and 24 communities; graph hash:
`8a3f3dfc159d172f6beb9c2fc9d9d76a83668597fdeb3f11436eb6860ab1406b` at
`/tmp/prism-step5-graphify/src/graphify-out/graph.json`. Its query identified
the provider and transformation symbols but did not prove a runtime call edge
to `structural_aligner.py`; direct imports decide that question.

## Provider and score contracts

| Provider | Current record shape / score meaning | Gate behavior | Status |
|---|---|---|---|
| TMalign | match count, transform, match map, parsed `tm_score`, alignment lengths, status, aligner | TM threshold plus match/coverage checks in transformation | maintained, tested |
| GTalign | TMalign-like fields plus `tm_score_ref/query` and raw-output hash | TM threshold plus match/coverage checks | maintained, tested; no explicit score contract field |
| MultiProt current | `tm_score = max(0, 1 - kabsch_rmsd/10)`; `tm_score_contract=multiprot_kabsch_rmsd_proxy`; native match/coverage gate contract | native match count and coverage, not proxy TM | maintained but semantically non-TM |
| MultiProt legacy | legacy native score contract | legacy behavior | retained/conditional only |
| USalign candidate | 214-match fixture, first parsed USalign TM score, transform, RMSD, status, aligner | candidate adapter tests; no explicit score contract field | isolated, bounded only |

The MultiProt proxy must not be compared numerically with TMalign, GTalign,
or USalign TM-scores. Earlier retained calibration succeeded for 47/56 pairs;
the proxy-to-true-TM correlation was approximately 0.047, and only two pairs
passed the native true-TM/match/coverage criteria. This is calibration
evidence, not a fresh full-pipeline run.

## Attrition and possible silent-drop locations

The retained 1gte ledger records 76,248 potential alignment sides, 2,997
alignment JSON records, 560 orientation attempts, 560 transformation
attempts, 489 clash rejections, and 71 final passes. All 2,997 loaded records
passed the recorded transform-field and match/TM/coverage checks. There are
1,163 structural missing-partner records. The 1,877 residual after consuming
the orientation attempts is explicitly not classified as structural orphanage.

The difference between potential sides and written alignment records is
73,251 observed missing/unwritten records. The evidence does not establish
whether those records represent missing inputs, provider failures, parse
failures, or an upstream enumeration mismatch. The next run must persist a
record for each such reason rather than infer one from absence.

Potential non-observable drop points in the current code are:

- provider subprocess failure or unparseable output producing an empty
  alignment record;
- missing alignment JSON, which `transformation.py` can substitute with an
  in-memory empty record;
- residue mapping or invalid-transform failure returned from transformation;
- match/TM/coverage gate failure;
- clash rejection after transformation;
- optional selector/refiner failures that are not represented by a complete
  end-to-end stage ledger;
- external Rosetta commands launched with `os.system` without per-candidate
  return-code or score-gate records.

The existing candidate-audit status vocabulary covers generated,
alignment-failed, transformation-failed, clash-rejected, refinement-failed,
and refinement-accepted, but does not by itself distinguish all requested
`missing_input`, `alignment_unavailable`, `parse_failed`, `residue_mapping_failed`,
`invalid_transform`, `score_gate_failed`, `missing_contact`,
`score_gate_rejected`, `not_configured`, and `skipped_by_configuration` cases.

## Refiner, ranking, and evaluator boundaries

- Canonical defaults remain TMalign, NACCESS, refinement enabled with external
  Rosetta, ranking disabled, baseline ranking, and `top_k=5`.
- Retained six-pair comparison evidence reports 144 alignment records, 14
  transformation records, 22 external-Rosetta outputs, and 21 FiberDock
  outputs. It is historical and not a matched DockQ quality panel.
- FiberDock's `fiberdock_energies.ref` versus guessed `fd_params.ref` issue was
  reproduced and corrected only in an isolated replay; no canonical source
  promotion occurred.
- The isolated PRODIGY adapter produced one selected candidate on success and
  preserved both candidates when PRODIGY failed for no contacts. This validates
  bounded state/failure handling, not ranking quality.
- The current lightweight `compare.py` DockQ wrapper does not persist the full
  native hash/evaluator-version contract required for a reproducible scientific
  ranking comparison; stronger benchmark scripts exist but have not been run
  as a fresh matched panel.

## Validation boundaries

| Layer | Evidence | Status |
|---|---|---|
| Transformation helpers and thresholds | canonical focused tests, 22 passed | `TESTED` |
| USalign provider compatibility | isolated tests, 20 passed; Slurm 1656966; 214 matches; TM 0.98359 | `VALIDATED` within fixture |
| PRODIGY adapter states | isolated tests, 9 passed; Slurm 1656992/1656993 | `VALIDATED` within fixture |
| Retained ranking load | baseline 4 selected/2 refined vs ranked top-1 1/1 | `VALIDATED` as retained load evidence |
| Fresh current-tree three-arm smoke | job 1657005, no durable output; controllers unavailable | `UNKNOWN` |
| GPU/CPU/multithread GTalign comparison | no Slurm execution in this run | `UNRUN` |
| Native DockQ quality comparison | no fresh matched panel | `UNKNOWN` |

The complete machine-readable counts and provenance are in
`tool-matrix.json`, `alignment-stage-ledger.json`, `refiner-stage-ledger.json`,
`matched-candidate-ledger.json`, and `no-drop-ledger.json` in this directory.
