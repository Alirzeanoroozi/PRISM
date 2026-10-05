# Bounded goal

Validate and document the PRISM-prescript DockQ parser and canonical benchmark
scorer corrections without overwriting historical benchmark evidence.

Success criteria:

- Complete mappings use bounded `GlobalDockQ`.
- Pairwise recovery is explicitly `scored_cross_only`.
- Valid-unscored JSON and failures remain auditable.
- No-align mapping checks protect residue correspondence.
- CPU use and raw-command provenance are explicit.
- Focused tests and syntax checks pass.
