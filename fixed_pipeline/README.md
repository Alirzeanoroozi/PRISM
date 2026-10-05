# Fixed Pipeline — Isolated Bug-Fix Workspace

This directory contains fixed, backportable versions of PRISM pipeline
modules, each addressing specific bugs identified in the main pipeline.
None of these files modify the original pipeline sources.

## What's Fixed

| Module | Bug | Fix |
|--------|-----|-----|
| `src/rosetta_refinement.py` | Shell injection via `os.system()` | Replaced with `subprocess.run()` + `shutil` |
| `src/rosetta_refinement.py` | Bare exception handlers | Structured error codes + logging |
| `src/rosetta_refinement.py` | Race condition on `_0001.pdb` | Unique output paths |
| `src/alignment.py` | Truncated mapping silently passed downstream | `mapping_truncated` check in transformer |
| `src/alignment_gtalign.py` | Missing files silently degrade coverage | Counted summary at end |
| `src/transformation.py` | `apply_tm_transform()` swallows all errors | Structured return with diagnostics |
| `src/naccess_utils.py` | FreeSASA subprocess fails silently | Pre-check before pipeline start |
| `src/hotspot.py` | `ATOM_DICT` undefined (dead code) | Defined + re-enabled filter |
| `src/alignment_multiprot.py` | Magic threshold `5`, seccomp exit 159 | Configurable threshold + diagnostic |

## Usage

```bash
# Deploy fixed modules by copying them over originals:
cp fixed_pipeline/src/*.py src/

# Or run regression tests:
cd fixed_pipeline
python -m pytest tests/
```

## Verification

Each fix includes:
- A test in `tests/` that validates the fix.
- Error-handling improvements that maintain backward compatibility.
- Structured logging where appropriate.
