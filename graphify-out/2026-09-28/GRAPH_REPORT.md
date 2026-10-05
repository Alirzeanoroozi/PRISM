# Graph Report - /scratch/rshadi25/GitHub/PRISM-prescript  (2026-09-17)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 584 nodes · 1113 edges · 26 communities (23 shown, 3 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 21 edges (avg confidence: 0.69)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `1026a7fc`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- run_evidence.py
- transformation.py
- normalize_target_id
- compare.py
- stepwise_analysis.py
- pyrosetta_refinement.py
- prodigy_ranker.py
- StructuralAligner
- alignment_multiprot.py
- multiprot_pyrosetta.py
- declare_contract
- alignment_gtalign.py
- irmsd.py
- refine_pairs
- validate_gate
- compute_tm_score
- extract_sequences.py
- pipeline_inputs.py
- ranking_data.py
- ranking_metrics.py
- get_asa_complex
- residue_contact_model.py
- plot_num_chains_distogram.py
- freesasa_runner.py
- eda/__init__.py

## God Nodes (most connected - your core abstractions)
1. `normalize_target_id()` - 17 edges
2. `build_execution_attempt()` - 17 edges
3. `main()` - 14 edges
4. `StructuralAligner` - 14 edges
5. `ArtifactObservation` - 13 edges
6. `validate_before_consume()` - 13 edges
7. `alignment_inventory()` - 13 edges
8. `process_pair_for_template()` - 13 edges
9. `observe_artifact()` - 11 edges
10. `append_artifact_observation()` - 11 edges

## Surprising Connections (you probably didn't know these)
- `declare_contract()` --calls--> `canonical_json()`  [INFERRED]
  run_identity.py → provenance/run_evidence.py
- `init_run()` --calls--> `canonical_json()`  [INFERRED]
  run_identity.py → provenance/run_evidence.py
- `init_run()` --calls--> `_make_run_id()`  [INFERRED]
  run_identity.py → provenance/run_evidence.py
- `declare_contract()` --calls--> `build_declared_contract()`  [INFERRED]
  run_identity.py → provenance/run_evidence.py
- `init_run()` --calls--> `build_execution_attempt()`  [INFERRED]
  run_identity.py → provenance/run_evidence.py

## Import Cycles
- None detected.

## Communities (26 total, 3 thin omitted)

### Community 0 - "run_evidence.py"
Cohesion: 0.06
Nodes (70): main(), Namespace, Write an artifact observation to the ledger., write_artifact(), Exception, Canonical Phase 1 run identity and artifact-evidence contracts., append_artifact_observation(), ArtifactLedgerError (+62 more)

### Community 1 - "transformation.py"
Cohesion: 0.06
Nodes (49): alignment_features(), CandidateAudit, CandidateRecord, Any, Leakage-safe candidate audit records for PRISM ranking experiments. The audit…, Extract stable, numeric features from one alignment JSON object., Append candidate records without changing pipeline acceptance logic., record_alignment_pair() (+41 more)

### Community 2 - "normalize_target_id"
Cohesion: 0.09
Nodes (43): analyse_pdb(), run_analysis(), get_contacts(), center_of_mass(), contacting_residues(), fetch_all_atoms_coordinates(), get_contact_potentials(), hotspot_creator() (+35 more)

### Community 3 - "compare.py"
Cohesion: 0.07
Nodes (47): _ca_by_key(), ca_rmsd(), compare_and_summarize(), compare_pair(), compare_pairs_from_outputs(), _compare_task(), get_trimmed_native_pdb(), _parse_output_pdb_path() (+39 more)

### Community 4 - "stepwise_analysis.py"
Cohesion: 0.09
Nodes (49): _alignment_contract_status(), alignment_inventory(), _asset_origins(), asset_provenance_manifest(), _atom_records(), _audit_index(), _ca_residue_count(), candidate_key() (+41 more)

### Community 5 - "pyrosetta_refinement.py"
Cohesion: 0.10
Nodes (37): get_contacts_from_atom_lines(), Write contacting residue-number pairs from two PDB atom-line groups., _base_report(), _combine_partners(), _command_metadata(), _environment_metadata(), _file_hashes(), _import_report() (+29 more)

### Community 6 - "prodigy_ranker.py"
Cohesion: 0.09
Nodes (38): biological_baseline_score(), _bounded(), _optional_bounded(), rank_candidates(), Deterministic biological candidate baseline used before ML reranking., Score a candidate without allowing failed rows into the ranking. TM-score and…, Return complete candidates sorted by baseline score, best first., _build_audit_index() (+30 more)

### Community 7 - "StructuralAligner"
Cohesion: 0.11
Nodes (19): Enum, AlignmentResult, AlignmentTool, Perform pairwise structural alignment (direct TMalign replacement). Args:…, Run GTalign (exact TMalign replacement), Run US-align (universal alignment), Run Foldseek in pairwise global alignment mode, Stage 1: Fast prefilter using SSAlign for large databases. Returns list of hit… (+11 more)

### Community 8 - "alignment_multiprot.py"
Cohesion: 0.12
Nodes (27): align(), _align_one(), extract_chain_and_res_ids(), _has_ca_atoms(), iter_bounded_results(), align_multiprot(), _align_one(), _check_seccomp() (+19 more)

### Community 9 - "multiprot_pyrosetta.py"
Cohesion: 0.18
Nodes (18): apply_tm_transform(), build_docked_complex(), find_staged_pdb(), get_chain_mapping(), main(), parse_multiprot_output(), parse_multiprot_transforms(), Path (+10 more)

### Community 10 - "declare_contract"
Cohesion: 0.20
Nodes (15): declare_contract(), _git_source_inventory(), init_run(), _load_pairs(), _load_template_inventory(), main(), Any, Namespace (+7 more)

### Community 11 - "alignment_gtalign.py"
Cohesion: 0.23
Nodes (14): align_gtalign(), build_match_dict_from_aligned_sequences(), extract_chain_and_res_ids(), _extract_floats(), _extract_gtalign_alignment_seq(), _extract_gtalign_query_path(), _parse_gtalign_hit_block(), parse_gtalign_hits() (+6 more)

### Community 12 - "irmsd.py"
Cohesion: 0.25
Nodes (13): almostIdentical(), calcRotMatAndEigenValues(), combineCoords(), getChainIndices(), getCommonResiduesList(), getInterfaceResidues(), getIRMSD(), getRMSD() (+5 more)

### Community 13 - "refine_pairs"
Cohesion: 0.19
Nodes (14): _add_hydrogens(), _build_fiberdock_params(), _check_tools(), _create_ca_pdb(), Run Normal Mode Analysis., Build FiberDock parameter file. receptor_hb/ligand_hb: Paths to hydrogenated…, Run FiberDock energy calculation. FiberDock must run from its own directory to…, FiberDock refinement entry point. Args: passed_pairs: List of (ligand_pdb,… (+6 more)

### Community 14 - "validate_gate"
Cohesion: 0.33
Nodes (9): _expected_inventory_error(), _load_expected_inventory(), main(), Namespace, Path, Load the declared artifact inventory used by the consumer gate. JSON accepts…, Emit a fail-closed result when the declared inventory is unusable., Run validation gate on completed run. (+1 more)

### Community 15 - "compute_tm_score"
Cohesion: 0.31
Nodes (8): compute_tm_score(), get_ca_coords(), parse_mp_key(), process_alignment_file(), Add true_tm_score to an alignment JSON., Parse MultiProt key like 'C.T.110' -> (chain, resnum, resname3), Return dict: (chain, resnum, resname) -> CA coord array, Compute true TM-score from MultiProt match_dict.

### Community 16 - "extract_sequences.py"
Cohesion: 0.36
Nodes (7): build_sequence_dict(), main(), parse_template(), Extract sequences from template PDBs using Biopython (ATOM records). Requires…, Split a template id like ``1a0dAB`` into (pdb_id, chain1, chain2)., Return {chain_id: one-letter sequence} for all chains (model 0, pdb-atom)., read_all_chains_from_pdb()

### Community 17 - "pipeline_inputs.py"
Cohesion: 0.36
Nodes (7): normalize_template_ids(), Input and template-panel selection helpers for the current PRISM CLI., Normalize template IDs from CLI tokens or a plain-text manifest. Current PRISM…, Read a six-character template ID per non-comment line., Select the template panel while preserving the current default.…, read_template_list(), select_templates()

### Community 18 - "ranking_data.py"
Cohesion: 0.29
Nodes (7): grouped_split(), native_like_label(), Validation and leakage-safe splitting for candidate ranking tables., Return the documented CAPRI native-like label., Validate and copy labeled rows; unlabeled rows are rejected explicitly., Split row indices by native complex, never by individual decoy rows., validate_training_rows()

### Community 19 - "ranking_metrics.py"
Cohesion: 0.48
Nodes (6): enrichment_factor(), _label(), Dependency-free ranking metrics for PRISM candidate experiments., Spearman correlation between score rank and binary biological label., spearman_score(), top_k_success()

### Community 20 - "get_asa_complex"
Cohesion: 0.43
Nodes (6): get_asa_complex(), get_asa_flat(), Compute relative ASA for the requested chains of a target. `target` is…, Return a flat dict keyed by `RESNAME_RESNUMBER_CHAINID`. Useful for the hotspot…, Resolve the path to a target's PDB file. `pdb_root` is the directory containing…, _resolve_pdb_path()

### Community 21 - "residue_contact_model.py"
Cohesion: 0.50
Nodes (3): build_model(), Optional residue-contact model for the Stage 2 PRISM experiment. PyTorch is…, Build a small contact-message-passing model when PyTorch is available.

## Knowledge Gaps
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `sha256_file()` connect `stepwise_analysis.py` to `declare_contract`?**
  _High betweenness centrality (0.166) - this node is a cross-community bridge._
- **Why does `declare_contract()` connect `declare_contract` to `run_evidence.py`?**
  _High betweenness centrality (0.136) - this node is a cross-community bridge._
- **Why does `_load_template_inventory()` connect `declare_contract` to `stepwise_analysis.py`?**
  _High betweenness centrality (0.082) - this node is a cross-community bridge._
- **Should `run_evidence.py` be split into smaller, more focused modules?**
  _Cohesion score 0.06288448393711552 - nodes in this community are weakly interconnected._
- **Should `transformation.py` be split into smaller, more focused modules?**
  _Cohesion score 0.061016949152542375 - nodes in this community are weakly interconnected._
- **Should `normalize_target_id` be split into smaller, more focused modules?**
  _Cohesion score 0.09049773755656108 - nodes in this community are weakly interconnected._
- **Should `compare.py` be split into smaller, more focused modules?**
  _Cohesion score 0.06787330316742081 - nodes in this community are weakly interconnected._