# Bug Fix Manifest

| # | File | Bug | Severity | Fixed? |
|---|------|-----|----------|--------|
| 1 | `rosetta_refinement.py` | **Shell injection** via `os.system()` — all Rosetta commands use f-string interpolation of file paths | 🔴 CRITICAL | ✅ |
| 2 | `rosetta_refinement.py` | **Bare exception handlers** — all errors caught as generic `Exception`, return bare `"-"` string | 🟠 HIGH | ✅ |
| 3 | `rosetta_refinement.py` | **Race condition** — Rosetta output `_0001.pdb` naming causes concurrent runs to clobber each other | 🟠 HIGH | ✅ |
| 4 | `rosetta_refinement.py` | **Score parsing fragile** — column index positions assumed without validation | 🟡 MEDIUM | ✅ |
| 5 | `alignment.py` | **Mapping truncation silent** — truncated alignment passed to transformer with `"mapping_truncated"` status never checked downstream | 🟠 HIGH | ✅ |
| 6 | `alignment_gtalign.py` | **Missing file coverage invisible** — missing query/template files printed but never counted | 🟡 MEDIUM | ✅ |
| 7 | `alignment_gtalign.py` | **No hit summary** — total hits, skips, and written JSONs not reported after run | 🟢 LOW | ✅ |
| 8 | `transformation.py` | **Transform failure silent** — `apply_tm_transform()` catches all exceptions, returns `False` with no diagnostic | 🟠 HIGH | ✅ |
| 9 | `transformation.py` | **NaN coordinates undetected** — rotation/translation with invalid matrices produce NaN coords silently | 🟠 HIGH | ✅ |
| 10 | `naccess_utils.py` | **FreeSASA env mismatch** — subprocess fails if FreeSASA in different conda env, no pre-flight check | 🟡 MEDIUM | ✅ |
| 11 | `naccess_utils.py` | **NACCESS binary missing** — stale path produces unhelpful error | 🟢 LOW | ✅ |
| 12 | `hotspot.py` | **ATOM_DICT filter disabled** — commented out with TODO, all atoms included in hotspot features | 🟡 MEDIUM | ✅ |
| 13 | `alignment_multiprot.py` | **Magic number 5** — `largest_solution < 5` hardcoded, not configurable | 🟢 LOW | ✅ |
| 14 | `alignment_multiprot.py` | **Seccomp exit 159 undiagnosed** — 32-bit binary blocked by seccomp, no useful error message | 🟡 MEDIUM | ✅ |

## Verification

```bash
# Run all tests:
bash fixed_pipeline/run_tests.sh

# Deploy fixes (backup originals first):
bash fixed_pipeline/deploy_fixes.sh --backup

# Run original pipeline tests to confirm no regression:
bash benchmark/scripts/run_stable_checks.sh
```
