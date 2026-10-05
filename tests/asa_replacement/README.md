# ASA Backend Replacement Test

This folder benchmarks residue-level solvent accessibility backends against
`Naccess` for PRISM-compatible structures.

Backends covered here:

- `Naccess`
- `FreeSASA`
- `RustSASA`

The harness extracts only the requested PRISM chains into a clean temporary PDB
before calling the backends. This keeps the input consistent across tools and
avoids parser issues on verbose source PDB headers.

Run with the dedicated Python 3.11 env because `rust-sasa-python` requires it:

```bash
tests/asa_tools_py311/bin/python tests/asa_replacement/compare_backends.py \
  --structure-kind target \
  --structure-id 1FGNH \
  --output tests/asa_replacement/output/1FGNH.backends.json
```

Template example:

```bash
tests/asa_tools_py311/bin/python tests/asa_replacement/compare_backends.py \
  --structure-kind template \
  --structure-id 2AI9AB \
  --output tests/asa_replacement/output/2AI9AB.backends.json
```
