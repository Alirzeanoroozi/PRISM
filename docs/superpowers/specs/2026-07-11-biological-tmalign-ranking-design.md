# Biological TMalign Candidate Ranking and Learning Design

## Objective

Improve the biological usefulness of the current TM-align + Rosetta pipeline
without replacing TM-align prematurely. TM-align remains the deterministic
candidate generator and transformation source. A new, auditable ranking layer
will prioritize candidates that are more likely to form a biologically
meaningful interface. The legacy MultiProt + FiberDock pipeline remains an
independent baseline for comparison.

## Approach Options

1. **Rule-only calibration**: add interface and surface filters to TM-align.
   This is low risk and interpretable, but cannot learn interactions between
   weak signals and cannot estimate confidence.
2. **Tabular candidate reranker (recommended first stage)**: train a calibrated
   gradient-boosting model on TM-align, interface, surface, clash, and Rosetta
   features. This needs modest data, is explainable with feature importance,
   and can be evaluated without changing generated structures.
3. **Residue/contact deep model (second stage)**: train a residue graph model
   or equivariant graph model on candidate interfaces, optionally augmented
   with protein-language-model embeddings. This can capture discontinuous
   interfaces but requires substantially more labeled and diverse complexes.

## Staged Design

### Stage 0: Candidate audit and reproducibility

Create one row per candidate orientation and retain rows for failures. Required
identity fields are query proteins, template, template chains, orientation,
source pipeline, input/native complex identifier, and repository/tool versions.

Required feature and status fields include:

- TM-align match count, TM-score, alignment coverage, and matched-residue
  percentage for each chain;
- transformed interface residue count, residue mapping coverage, surface-area
  features, contact count, clash count, and transformation status;
- Rosetta prepack/docking status and interaction score;
- explicit terminal status such as `generated`, `alignment_failed`,
  `transformation_failed`, `clash_rejected`, `refinement_failed`, or
  `refinement_accepted`.

No candidate, chain, orientation, or failure may be silently dropped.

### Stage 1: Supervised candidate reranking

For benchmark complexes with native structures, compute DockQ, interface RMSD,
fraction of native contacts, and ligand RMSD where applicable. Define the
primary binary label as `native_like = DockQ >= 0.23`; retain continuous iRMSD
and DockQ for regression and ranking analysis.

The first model will be a scikit-learn gradient-boosting baseline or an
XGBoost/LightGBM model only if the environment already provides it. It will
rank candidates, not create new transformations. A probability calibration
step will be fitted on validation data only. If model confidence is low or the
candidate is out of distribution, the existing deterministic ranking remains
the fallback.

Splits must be grouped by native complex and, where possible, by sequence
clusters at approximately 40% identity. All decoys derived from one complex,
template family, or parent structure must remain in one split. Feature
normalization, calibration, and feature selection must be fitted inside each
training fold.

### Stage 2: Deep-learning model

Only after Stage 1 establishes a clean labeled table, evaluate a residue/contact
model. Residues are graph nodes; edges represent spatial contacts, with the
contact cutoff recorded as a configuration value. Initial node features are
sequence identity class, physicochemical properties, solvent exposure, and
interface membership. Optional ESM embeddings are an ablation, not a required
dependency. The model should predict DockQ/native-like probability and/or
iRMSD, with class-balanced sampling or balanced bootstrap ensembles for the
strongly imbalanced decoy set.

An equivariant architecture is preferred only if a reproducible environment,
GPU allocation, and enough diverse training complexes are available. It must be
compared against the tabular model and a TM-score-only baseline, not evaluated
in isolation.

## Evaluation and Acceptance Criteria

Every experiment reports the exact split manifest, features, labels, model
version, random seed, and excluded/failed rows. Primary metrics are:

- top-1 native-like success and DockQ;
- enrichment factor and success rate at top 1%, 5%, and 10%;
- Spearman correlation for ranking and MAE/RMSE for iRMSD;
- precision-recall AUC and probability calibration;
- coverage and performance of the deterministic fallback path.

The learned stage is accepted only if it improves held-out top-k biological
quality over the current baseline without reducing candidate generation
coverage, and if the improvement survives a complex/sequence-grouped split.
Results must include an ablation of TM-score, interface features, surface
features, clash features, and Rosetta features. The legacy MultiProt +
FiberDock pipeline is reported separately; it is not mixed into training data.

## Risks and Controls

- **Leakage**: grouped complex and sequence-cluster splits; immutable split
  manifests.
- **Too few positives**: report prevalence and confidence intervals; do not
  claim accuracy from one complex.
- **Rosetta score shortcut**: ablate Rosetta features and report pre-Rosetta
  ranking separately.
- **Rigid-alignment blind spots**: mark flexible/discontinuous interfaces and
  quantify them as an explicit subgroup.
- **Distribution shift**: retain the deterministic fallback and report an OOD
  score or training-support distance.
- **Reproducibility**: preserve candidate rows, failed statuses, tool versions,
  seeds, and Slurm job metadata.

## Implementation Boundary

The first implementation consists of the candidate audit schema, deterministic
feature extraction, split/label generation, baseline ranking report, and tests.
Deep-learning training is a separate follow-up after the audit table and Stage
1 ablations pass. No default threshold or pipeline acceptance behavior changes
until held-out evidence supports the change.
