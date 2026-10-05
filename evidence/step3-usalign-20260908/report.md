# PRISM-prescript — Step 3 USalign report

## Outcome

The existing USalign installation is explicitly `AVAILABLE` and compatible
with the PRISM alignment JSON contract for a bounded real case. The safe
integration shape is an opt-in provider branch, not a default change and not
the currently disconnected `StructuralAligner` abstraction.

No canonical PRISM-prescript file was modified. The implementation and tests
remain in the isolated worker worktree for review/promotion only.

## Root-cause findings

### USalign CLI/output mismatch — `EVIDENCE`

USalign's default full output includes the alignment summary and TM-scores but
does not emit the transform matrix. The installed help explicitly documents
`-m -` as “print matrix to stdout.” A real Slurm run confirmed that the matrix
appears only with `-m -` and uses five tokens per row: row index, translation,
and three rotation values.

Therefore, the original minimal two-argument adapter contract was insufficient:
it could not populate `rotation_mat` or `translation`. The isolated adapter
now invokes `-outfmt -1 -m -` and parses the observed row layout.

### Entry-point disconnection — `OBSERVATION` / `INFERENCE`

`prism.py` directly selects `src.alignment`, `src.alignment_gtalign`, or
`src.alignment_multiprot`, then directly calls `transformer`. It does not
import or call `StructuralAligner` from `src/structural_aligner.py`.

Graphify’s shortest path between `StructuralAligner` and `transformer()` is a
generic shared-import path through `json`, not a maintained pipeline call edge.
The adapter was therefore tested at the maintained CLI provider boundary.

## Compatibility evidence

The live raw probe on Slurm job `1656964` produced:

- `Aligned length=214`
- `RMSD=0.71`
- `TM-score=0.98359`
- a three-row transform matrix

Parsing and mapping the same output produced 214 residue correspondences,
with sample direction `L.D.1 -> A.D.1`. Applying the parsed matrix to the 428
CA coordinates of the two probe structures produced RMSD `0.7500643` and max
CA error `2.07497`, supporting the documented equation and orientation.

The corrected live adapter smoke on Slurm job `1656966` wrote a successful
PRISM-compatible record with `status=success`, `aligner=USalign`,
`match_count=214`, `tm_score=0.98359`, and the same transform/mapping values.

## Isolated implementation state

The candidate is `IMPLEMENTED` and `TESTED` in the isolated worktree. The
focused test suite has 17 passing tests; the complete isolated test directory
has 20 passing tests. Tests cover parser formats, transform fields, mapping
direction, missing binary/input behavior, fail-closed malformed output, CLI
selection, default `tmalign`, and absence of a shared USalign alignment alias.

The bounded live probe supports `VALIDATED_BOUNDED`. Source review and
Graphify review support `REVIEWED_ISOLATED_ONLY`. This is not a claim that the
dirty canonical checkout has been promoted or that a full scientific panel is
validated.

## Acceptance gate

| Requirement | Result | Evidence |
|---|---|---|
| USalign classified | `AVAILABLE` | existing wrapper, help, version, Slurm probe |
| flags/output/matrix recorded | pass | `-outfmt -1 -m -`, real five-column rows |
| chain/residue mapping recorded | pass | 214 mappings, `L -> A` sample |
| transform convention recorded | pass | equation and 428-CA RMSD check |
| failure behavior recorded | pass | 17 focused tests and fail-closed records |
| `structural_aligner.py` connection determined | pass | direct source + Graphify query |
| default behavior preserved | pass | default parser remains `tmalign`; canonical status hash unchanged |
| canonical project files untouched | pass | 373-entry status hash unchanged |
| full real-project validation | not claimed | explicitly outside bounded Step 3 scope |

## Decision

Step 3’s bounded compatibility gate passes. Stop at this step. Do not merge or
copy the isolated candidate without explicit promotion authorization and a
fresh validation against the current dirty canonical baseline.
