# BM55 native-chain audit

This audit compares the selected PRISM-prescript BM55 manifest with the
native-chain assignments encoded in the benchmark `Complex` fields from
`T_Rigid.csv`, `T_medium.csv`, and `T_difficult.csv`. It also inspects the
actual chain IDs present in each assembled native PDB on KUACC.

## Result

- 69,895 selected candidate rows across 257 unique BM55 cases.
- 257/257 case mappings agree with the benchmark `Complex` assignments after
  removing trailing benchmark annotation markers such as `*`.
- 233 assembled native PDBs have exactly the expected chain set.
- 20 contain extra blank or auxiliary chains, but all expected chains are
  present; explicit DockQ mappings can still restrict scoring to the intended
  receptor/ligand groups.
- 4 lack an expected native chain and require assembly repair before scoring:
  - `medium_3aad_052`: expected `A,D`, actual `A,B`.
  - `rigid_1oyv_069`: expected `B,I`, actual `A,I`.
  - `rigid_3p57_159`: expected `CD,P`, actual `AB,P`.
  - `rigid_3p57_160`: expected `IJ,P`, actual `AB,P`.

The `T_*.csv` files establish the expected labels but do not provide missing
coordinates. These four cases must therefore remain excluded from DockQ until
the correct native assemblies are recovered and hash-validated.

## Provenance

- Remote manifest: `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-large-refine-20260927/manifest/selected_candidates.csv`
- Remote manifest SHA-256: `16caa2a8e35b9c5b99e72d325540ed1901299a274faded199d89893b735dc5fc`
- Remote audit namespace: `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-large-refine-20260927-repair-v1/audit/`
- Local compact outputs: `audit_summary.json`, `native_chain_mapping_audit.csv`,
  and `native_chain_inventory.csv`.

No active jobs, raw benchmark files, or existing result checkpoints were
modified. The staged repair worker strips trailing benchmark annotation
markers from chain fields before DockQ preflight.
