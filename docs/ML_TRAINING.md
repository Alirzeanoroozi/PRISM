# PRISM Learning Workflow

The learning components are optional and do not change the default TM-align +
Rosetta pipeline. They require a labeled candidate table containing
`native_complex_id`, `dockq`, and the four TM-align feature columns.

## Build and Label Candidates

Use exact validated cases rather than scanning the historical artifact tree:

```bash
python benchmark/scripts/build_candidate_table.py \
  --query-pair 1FGNH,1TFHA \
  --template 1h5bAB \
  --native-complex-id benchmark-1ahw \
  --refinement-dir processed/rosetta_refinement/structures \
  --output results/candidates.csv
```

Attach DockQ/iRMSD scores produced by the benchmark scoring workflow:

```bash
python benchmark/scripts/attach_native_labels.py \
  results/candidates.csv scores.csv results/candidates_labeled.csv
```

Rows with missing or ambiguous score matches remain present and must not be
used as supervised labels.

## Stage 1 Training

Run this on a compute allocation with scikit-learn installed. The trainer uses
native-complex-grouped splits and writes both a pickle model and metrics JSON:

```bash
python benchmark/scripts/train_reranker.py \
  results/candidates_labeled.csv \
  --model results/reranker.pkl \
  --metrics results/reranker_metrics.json \
  --seed 0
```

The model is a CPU-safe `HistGradientBoostingClassifier`. It is a baseline,
not evidence of improvement, until held-out top-k DockQ/native-like metrics
beat the deterministic biological ranking baseline.

## Stage 2 Deep Model

The optional `src/residue_contact_model.py` requires PyTorch and is not wired
into production. It consumes residue node features and a contact edge list and
returns native-like logits plus an iRMSD estimate. Train it only after Stage 1
has established a leakage-safe table and after comparing against the tabular
baseline. GPU execution, sequence-cluster splits, class balancing, and model
checkpoint provenance must be recorded in the experiment output.

## Acceptance Gate

Do not enable learned ranking by default unless the held-out evaluation reports
an improvement in top-1/top-k DockQ and enrichment without reducing candidate
generation coverage. Report the TM-score-only baseline, deterministic
biological baseline, tabular model, and contact model separately.
