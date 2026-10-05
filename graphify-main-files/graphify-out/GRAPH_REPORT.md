# Graph Report - .  (2026-08-21)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 556 nodes · 1100 edges · 19 communities (14 shown, 5 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 17 edges (avg confidence: 0.5)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `1026a7fc`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- Community 0
- Community 1
- Community 2
- Community 3
- Community 4
- Community 5
- Community 6
- Community 7
- Community 8
- Community 9
- Community 10
- Community 11
- Community 12
- Community 13
- Community 14
- Community 15
- Community 16
- Community 17

## God Nodes (most connected - your core abstractions)
1. `TemplateGenerator` - 25 edges
2. `TransformFilter` - 19 edges
3. `main()` - 15 edges
4. `normalize_target_id()` - 14 edges
5. `StructuralAligner` - 14 edges
6. `SurfaceExtractor` - 12 edges
7. `refine_pairs()` - 11 edges
8. `process_pair_for_template()` - 11 edges
9. `_align_one()` - 10 edges
10. `select_top_candidates()` - 10 edges

## Surprising Connections (you probably didn't know these)
- `main()` --calls--> `run_analysis()`  [EXTRACTED]
  prism.py → src/analyse_pdbs.py
- `main()` --calls--> `select_top_candidates()`  [EXTRACTED]
  prism.py → src/candidate_selector.py
- `main()` --calls--> `compare_pairs_from_outputs()`  [EXTRACTED]
  prism.py → src/compare.py
- `main()` --calls--> `refine_pairs()`  [EXTRACTED]
  prism.py → src/fiberdock_refinement.py
- `main()` --calls--> `pdb_downloader()`  [EXTRACTED]
  prism.py → src/pdb_download.py

## Import Cycles
- None detected.

## Communities (19 total, 5 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.07
Nodes (57): current_tmalign_rosetta::bio_pdb, current_tmalign_rosetta::bio_pdb_polypeptide, current_tmalign_rosetta::gzip, current_tmalign_rosetta::json, current_tmalign_rosetta::os, current_tmalign_rosetta::pandas, analyse_pdb(), run_analysis() (+49 more)

### Community 1 - "Community 1"
Cohesion: 0.06
Nodes (56): current_tmalign_rosetta::concurrent_futures, current_tmalign_rosetta::datetime, build_parser(), main(), parse_bool(), Parse CLI booleans without Python's bool('false') trap., Append an opt-in, machine-readable pipeline stage event., record_stage_event() (+48 more)

### Community 2 - "Community 2"
Cohesion: 0.05
Nodes (49): current_tmalign_rosetta::argparse, current_tmalign_rosetta::freesasa, current_tmalign_rosetta::numpy, current_tmalign_rosetta::pathlib, main(), Namespace, Write an artifact observation to the ledger., write_artifact() (+41 more)

### Community 3 - "Community 3"
Cohesion: 0.07
Nodes (47): current_tmalign_rosetta::math, current_tmalign_rosetta::shlex, biological_baseline_score(), _bounded(), _optional_bounded(), rank_candidates(), Deterministic biological candidate baseline used before ML reranking., Score a candidate without allowing failed rows into the ranking. TM-score and… (+39 more)

### Community 4 - "Community 4"
Cohesion: 0.07
Nodes (15): old_multiprot_fiberdock::checktemplate, old_multiprot_fiberdock::flexiblerefinement, old_multiprot_fiberdock::pdbdownload, old_multiprot_fiberdock::preprocessor, TemplateChecker, FlexibleRefinement, Controller, PDBdownload (+7 more)

### Community 5 - "Community 5"
Cohesion: 0.09
Nodes (40): current_tmalign_rosetta::importlib, current_tmalign_rosetta::importlib_metadata, PathLike, current_tmalign_rosetta::platform, get_contacts_from_atom_lines(), Write contacting residue-number pairs from two PDB atom-line groups., _base_report(), _combine_partners() (+32 more)

### Community 6 - "Community 6"
Cohesion: 0.08
Nodes (39): alignment_features(), CandidateAudit, CandidateRecord, Any, Leakage-safe candidate audit records for PRISM ranking experiments. The audit…, Extract stable, numeric features from one alignment JSON object., Append candidate records without changing pipeline acceptance logic., record_alignment_pair() (+31 more)

### Community 7 - "Community 7"
Cohesion: 0.07
Nodes (23): old_multiprot_fiberdock::configparser, old_multiprot_fiberdock::fiberdockinterfaceextractor, old_multiprot_fiberdock::glob, old_multiprot_fiberdock::gzip, old_multiprot_fiberdock::maincontroller, old_multiprot_fiberdock::math, old_multiprot_fiberdock::mysqldb, old_multiprot_fiberdock::os (+15 more)

### Community 8 - "Community 8"
Cohesion: 0.10
Nodes (21): current_tmalign_rosetta::dataclasses, Enum, current_tmalign_rosetta::logging, AlignmentResult, AlignmentTool, Perform pairwise structural alignment (direct TMalign replacement). Args:…, Run GTalign (exact TMalign replacement), Run US-align (universal alignment) (+13 more)

### Community 9 - "Community 9"
Cohesion: 0.11
Nodes (27): current_tmalign_rosetta::bio_pdb_superimposer, current_tmalign_rosetta::csv, current_tmalign_rosetta::scratch_rshadi25_github_prism_prescript_graphify_main_files_current_tmalign_rosetta_src_eval_dockq_py, current_tmalign_rosetta::scratch_rshadi25_github_prism_prescript_graphify_main_files_current_tmalign_rosetta_src_eval_irmsd_backbone_py, _ca_by_key(), ca_rmsd(), compare_and_summarize(), compare_pair() (+19 more)

### Community 11 - "Community 11"
Cohesion: 0.18
Nodes (18): apply_tm_transform(), build_docked_complex(), find_staged_pdb(), get_chain_mapping(), main(), parse_multiprot_output(), parse_multiprot_transforms(), Path (+10 more)

### Community 13 - "Community 13"
Cohesion: 0.25
Nodes (8): current_tmalign_rosetta::hashlib, grouped_split(), native_like_label(), Validation and leakage-safe splitting for candidate ranking tables., Return the documented CAPRI native-like label., Validate and copy labeled rows; unlabeled rows are rejected explicitly., Split row indices by native complex, never by individual decoy rows., validate_training_rows()

### Community 14 - "Community 14"
Cohesion: 0.50
Nodes (3): build_model(), Optional residue-contact model for the Stage 2 PRISM experiment. PyTorch is…, Build a small contact-message-passing model when PyTorch is available.

## Knowledge Gaps
- **2 isolated node(s):** `prism.sh script`, `setup_cluster.sh script`
  These have ≤1 connection - possible missing edges or undocumented components.
- **5 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `TemplateGenerator` connect `Community 10` to `Community 4`, `Community 7`?**
  _High betweenness centrality (0.016) - this node is a cross-community bridge._
- **What connects `prism.sh script`, `setup_cluster.sh script` to the rest of the system?**
  _2 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.0710085933966531 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.059322033898305086 - nodes in this community are weakly interconnected._
- **Should `Community 2` be split into smaller, more focused modules?**
  _Cohesion score 0.0512987012987013 - nodes in this community are weakly interconnected._
- **Should `Community 3` be split into smaller, more focused modules?**
  _Cohesion score 0.06823529411764706 - nodes in this community are weakly interconnected._
- **Should `Community 4` be split into smaller, more focused modules?**
  _Cohesion score 0.07030527289546716 - nodes in this community are weakly interconnected._