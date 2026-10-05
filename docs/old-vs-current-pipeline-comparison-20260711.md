# Legacy MultiProt/FiberDock vs Current TM-align/Rosetta

## Scope and Evidence

The comparison uses the maintained legacy implementation in
`/scratch/rshadi25/GitHub/prism-oldversion` and the current implementation in
this repository. The two smoke runs use the same legacy fixture and template;
the current backend comparison uses a separate current-pipeline PDB because
the repositories do not currently expose one shared end-to-end fixture with
both template formats. Therefore, the legacy end-to-end counts are directly
paired, while the current result is a verified stage-level result.

Evidence is classified as:

- **Observed:** directly measured from source, artifacts, or Slurm logs.
- **Inferred:** follows from the observed control flow and thresholds.
- **Hypothesis:** plausible but requiring a controlled replay to quantify.

## Stage-by-Stage Comparison

| Stage | Legacy MultiProt + FiberDock | Current TM-align + Rosetta | Consequence for report count |
|---|---|---|---|
| Input/surface extraction | Configurable POPS, NACCESS, or FreeSASA; the validated NACCESS smoke produced both target surfaces. | NACCESS default, with new explicit FreeSASA backend. On Slurm job `1340649`, both produced 136 surface residues with Jaccard `1.0` for the tested target. | **Not the cause in the tested case.** Backend/radii differences remain a dataset-dependent hypothesis. |
| Structural alignment | MultiProt stores up to three parsed solutions per chain (`multiprotcount = 3`). In the validated run, each of four chain/target files contained 3 solutions. | TMalign runs once per target/template/chain and writes one JSON transform per pair. | **Primary confirmed cause:** legacy has multiple candidate transforms; current has one. |
| Alignment acceptance | Match count, match percentage, hotspot criterion, contact mapping, and clash filtering. No TM-score threshold. | Requires at least 15 matched residues, TM-score at least `0.5`, match-percentage threshold, and clash filtering. | **Confirmed additional loss mechanism:** valid-looking alignments below TM-score `0.5` are rejected by current. |
| Orientation/candidate expansion | MultiProt solutions are paired across chain sides and can generate many transformed candidates. | Two orientations are attempted, but each uses only the single alignment per chain. | Legacy candidate multiplicity can grow approximately as `solutions_chain1 x solutions_chain2`; current is bounded by two orientations. |
| Transformation filtering | Uses explicit hotspot matching, contact count, match thresholds, and clash limits from `prism.ini`. | Uses match count, TM-score, match percentage, and clash limits; `hotspot_analysis()` currently returns `True`. | Current hotspot logic is looser, so it does **not** explain fewer reports; it may increase current candidates relative to a fully enforced hotspot filter. |
| Refinement | FiberDock accepts a negative energy; the validated smoke model had energy `-45.23`. | Rosetta retains structures only when interaction score is at most `-5.0`. | **Confirmed potential loss mechanism:** Rosetta’s acceptance rule is stricter and its score is not numerically comparable to FiberDock energy. |
| Final output | Validated NACCESS and FreeSASA runs each produced 1 FiberDock model, 1 interface file, and energy `-45.23`. | Existing current source writes only accepted Rosetta structures and refinement energies; a same-input full comparison is not currently available. | Final counts cannot be compared scientifically until the same candidates are replayed through both refinement engines. |

## Measured Smoke Results

| Pipeline/backend | Slurm job | State | Models | Energy/result |
|---|---:|---|---:|---|
| Legacy MultiProt + FiberDock + NACCESS | `1340567` | COMPLETED | 1 | FiberDock energy `-45.23` |
| Legacy MultiProt + FiberDock + FreeSASA | `1340665` | COMPLETED | 1 | FiberDock energy `-45.23` |
| Current surface stage + NACCESS/FreeSASA comparison | `1340649` | COMPLETED | 136 residues each | Intersection 136; Jaccard `1.0` |

The legacy FreeSASA run initially failed to produce candidates because its RSA
writer used a different column layout from NACCESS while the legacy parser
used fixed-column offsets. The writer was corrected to emit NACCESS-compatible
columns; the rerun then completed successfully with the same FiberDock result.
This is an important pipeline-specific error distinction, not a biological
difference.

## Why TM-align Reports Fewer Predictions

### Confirmed causes

1. **Candidate multiplicity is intentionally reduced.** MultiProt retains up to three solutions per chain. TMalign retains one transform per pair. The legacy pipeline therefore explores alternatives that the current pipeline never creates.
2. **TM-score filtering removes candidates before transformation.** Current `TM_SCORE_THRESHOLD = 0.5` is an additional gate absent from the legacy MultiProt path.
3. **Rosetta acceptance is a separate stricter gate.** Current Rosetta requires interaction score `<= -5.0`; legacy FiberDock accepts any negative energy. The scores are tool-specific and cannot be compared by numeric magnitude.
4. **The output definitions differ.** A legacy prediction is accepted after FiberDock negative-energy validation; a current prediction additionally needs a Rosetta score, structure file, and accepted interaction score.

### Plausible causes requiring controlled replay

1. **Surface-set differences:** FreeSASA and NACCESS can differ because of atomic radii, residue reference values, CA main-chain/side-chain conventions, or parser formatting. The tested current target showed identical surface sets, but that does not prove equivalence for all proteins or protein-DNA cases.
2. **Template/interface set differences:** Different template lists, chain IDs, interface residue files, or template filtering can change the candidate universe before alignment.
3. **Alignment objective differences:** MultiProt multiple structural alignment and TMalign pairwise alignment optimize different objectives and can produce different residue mappings, transforms, RMSDs, and coverage even on identical structures.
4. **Hotspot/contact semantics:** Legacy hotspot criteria are explicit and configurable; current hotspot analysis is currently unconditional. Contact definitions and residue numbering can still differ downstream.
5. **Clash implementation differences:** Both use a 3 Angstrom-style clash rule, but the candidate structures, atom subsets, and exact rejection timing differ.
6. **External-tool failure handling:** Current and legacy wrappers may silently skip failed alignment/refinement outputs or produce missing score files. Missing, rejected, and never-generated candidates must be counted separately.
7. **Caching and stale artifacts:** Legacy alignment files and current JSON outputs can be reused across runs. A stale result directory can make apparent prediction counts reflect prior configuration rather than the current run.
8. **Duplicate suppression and orientation naming:** The two pipelines deduplicate and name candidates differently, so raw file counts may not equal unique biological predictions.
9. **Order and numerical behavior:** Different parsers, residue numbering, floating-point thresholds, and tie handling can move borderline candidates across acceptance boundaries.

## Controlled Experiments That Isolate the Difference

1. Run both pipelines on the same target PDBs, same template PDBs, same chain IDs, and the same surface files; record candidate counts after every stage.
2. Export every alignment candidate as a normalized record containing target, template, chain, match count, match percentage, TM-score/RMSD, transform ID, clash count, and status.
3. Replay all MultiProt candidates through the current transformation and Rosetta stages without changing their transforms. This isolates alignment multiplicity from refinement filtering.
4. Replay the single TMalign candidate through FiberDock. This isolates FiberDock versus Rosetta acceptance.
5. Run current with TM-score thresholds `0.0`, `0.4`, and `0.5`, while keeping all other thresholds fixed. The difference in accepted counts estimates the TM-score gate contribution.
6. Run current with one-to-three retained alignment hypotheses per chain, if available, to test whether multiplicity alone closes the count gap.
7. Record `generated`, `surface_passed`, `alignment_passed`, `transform_passed`, `refinement_attempted`, `refinement_accepted`, `failed`, and `rejected` counts separately. Never treat missing output as a biological negative.

## Bottom Line

The strongest evidence-backed explanation is not that TM-align is universally
worse. It is that the newer pipeline searches a smaller hypothesis space and
applies additional acceptance criteria, especially one-transform-per-chain,
TM-score `>= 0.5`, and Rosetta interaction score `<= -5.0`. The old pipeline’s
higher prediction count is therefore expected to be partly a multiplicity and
reporting-policy effect. A fair scientific comparison requires candidate-level
replay on identical inputs before attributing the remaining difference to
alignment quality.
