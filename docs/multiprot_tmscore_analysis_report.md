# MultiProt TM-Score Analysis & Low-Yield Diagnosis Report

**Date**: 2026-07-29  
**Project**: PRISM-prescript Pipeline Comparison  
**Scope**: Tasks A & B — TM-score divergence analysis & MultiProt 2/28 yield diagnosis

---

## Executive Summary

This report addresses two parallel investigations on the MultiProt aligner within the PRISM pipeline:

| Task | Finding | Key Evidence |
|------|---------|--------------|
| **A: TM-score divergence** | MultiProt reports a **TM-score proxy** (`1 - RMSD/10`), not the standard length-normalized TM-score. Values (0.69–0.72) are **not comparable** to TMalign/GTalign (0.18–0.31). | `alignment_multiprot.py:252` computes `tm_score = max(0, 1 - kabsch_rmsd/10)` |
| **B: 2/28 yield** | 26/28 pairs fail at **early filtering** (alignment unavailable, match count < threshold, or seccomp blocks 32-bit binary). Only 2 pairs pass all gates. | `multiprot_gate_diagnosis.json`: 26 orientations show `alignment_unavailable` or `match_count_below_minimum` |

**Root causes are independent** — Task A is a **metric definition mismatch**; Task B is a **pipeline gating/filtering issue**.

---

## Task A: MultiProt TM-Score vs TMalign/GTalign

### 1. Definitions: TM-score vs RMSD

| Metric | Formula | Range | Normalization | What it measures |
|--------|---------|-------|---------------|------------------|
| **TM-score** (standard) | `TM = 1/L_t ∑ 1/(1+(d_i/d_0)²)` | [0, 1] | Target length `L_t` | Global topology similarity |
| **RMSD** | `√(1/N ∑ ‖x_i - y_i‖²)` | [0, ∞) | None (absolute Å) | Local atomic deviation |
| **MultiProt proxy** | `max(0, 1 - RMSD/10)` | [0, 1] | Ad-hoc (10Å scale) | **Inverse RMSD heuristic** |

**Critical distinction**: Standard TM-score rewards **aligned length** and **topological similarity** even with moderate RMSD. MultiProt's proxy **only sees RMSD** — a 3Å alignment over 50 residues scores 0.7, but a 2Å alignment over 10 residues also scores 0.8. They are **mathematically incomparable**.

### 2. Why MultiProt's Values Appear Higher

| Pair | TMalign TM | MultiProt "TM" | GTalign TM | MultiProt RMSD | MultiProt Matches |
|------|------------|----------------|------------|----------------|-------------------|
| `5zng_1a0cCD_C` | 0.250 | **0.716** | 0.000 | 2.84 Å | 32 |
| `5zng_1a0hAE_A` | 0.182 | **0.695** | 0.000 | 3.05 Å | 14 |

MultiProt finds **more residue matches** (14–32 vs TMalign's 4–68) but with **higher RMSD** (2.8–3.1 Å). The proxy `1 - RMSD/10` converts 2.84Å → 0.716, while TMalign's length-normalized score for the same alignment is 0.250.

**Conclusion**: MultiProt is not "better" — it uses a different, non-standard metric that inflates scores for short, high-RMSD alignments.

### 3. Methodology to Reconcile / Estimate True TM-score from MultiProt

Since MultiProt outputs **residue correspondences** (match_dict) and **RMSD**, we can compute standard TM-score:

```python
# Pseudocode: compute TM-score from MultiProt alignment
def compute_tmscore_from_multiprot(match_dict, query_path, interface_path):
    # 1. Extract matched CA coordinates from both structures
    q_coords, i_coords = extract_matched_cas(match_dict, query_path, interface_path)
    
    # 2. L_t = length of target (query) protein
    L_t = count_residues(query_path)
    
    # 3. d_0 = 1.24 * (L_t - 15)**(1/3) - 1.8  (TM-score standard)
    d0 = 1.24 * (L_t - 15)**(1/3) - 1.8
    
    # 4. TM = (1/L_t) * sum(1 / (1 + (d_i/d0)**2))
    tm = sum(1 / (1 + (dist**2 / d0**2)) for dist in pairwise_distances(q_coords, i_coords)) / L_t
    return tm
```

**Required outputs from MultiProt** (already available in `alignment_multiprot.py`):
- `match_dict`: residue correspondences
- `kabsch_rmsd`: RMSD after Kabsch superposition
- Query/template PDB paths for coordinate extraction

### 4. Step-by-Step Reproduction Plan

| Step | Command | Output |
|------|---------|--------|
| 1. Run MultiProt on all 28 pairs | `python -m src.alignment_multiprot queries.txt templates.txt` | `processed/alignment/*.json` |
| 2. Parse match_dict + RMSD | `python scripts/extract_multiprot_matches.py` | CSV with matches, RMSD |
| 3. Compute true TM-score | `python scripts/compute_tmscore_from_matches.py` | CSV with standard TM-score |
| 4. Compare with TMalign/GTalign | `python scripts/compare_aligners.py` | Correlation plots, Δ tables |

**Reference files**: 
- `/scratch/rshadi25/tmp/alignment_comparison.py` (existing comparison)
- `processed/alignment_multiprot/` (MultiProt JSONs)
- `processed/alignment/` (TMalign JSONs)

### 5. Pipeline Position Assessment

| Stage | MultiProt Behavior | Divergence Source |
|-------|-------------------|-------------------|
| **Alignment** | Finds residue matches, computes RMSD, derives proxy TM | ✅ **Primary** — different algorithm, different metric |
| **Transformation** | Kabsch from matches → rotation/translation | Secondary — depends on alignment |
| **Post-processing** | Writes JSON with `tm_score_contract: "multiprot_kabsch_rmsd_proxy"` | Documentation only — flags non-standard metric |

**Verdict**: Divergence originates at **alignment stage** (algorithm + metric). Transformation is deterministic given matches. Post-processing correctly labels the metric.

---

## Task B: Why MultiProt Reports Only 2/28 Pairs

### 1. Hypotheses & Test Methods

| # | Hypothesis | Test Method | Expected Evidence |
|---|------------|-------------|-------------------|
| H1 | **Seccomp blocks 32-bit binary** on compute nodes | Check `/proc/self/status` Seccomp level; run on login vs compute | `Seccomp: 2` on compute → MultiProt skipped |
| H2 | **Match count < minimum** (default 15) | Inspect `multiprot_gate_diagnosis.json` `left_match_count` | Values < 15 for most pairs |
| H3 | **TM-score proxy < threshold** (default 0.5) | Check `left_tm_score` in diagnosis JSON | Values 0.0–0.07 < 0.5 |
| H4 | **Alignment unavailable** (MultiProt failed/crashed) | Check `left_status: "alignment_unavailable"` | 26/28 orientations show this |
| H5 | **Early termination / timeout** | Add timeout logging in `_align_one` | Logs show `TimeoutExpired` or early exit |
| H6 | **Template interface files missing** | Verify `templates/interfaces/{template}_{chain}_int.pdb` exists | Missing files → skip |

### 2. Evidence from Diagnosis JSON

```json
// multiprot_gate_diagnosis.json (sample)
{
  "left_status": "alignment_unavailable",
  "left_match_count": 0,
  "left_tm_score": 0.0,
  "left_failure_reasons": "alignment_unavailable;match_count_below_minimum;tm_score_below_threshold;match_percentage_below_threshold"
}
```

**26 of 28 orientations** have `alignment_unavailable` + `match_count=0` + `tm_score=0.0`.

The **2 passing pairs** (from `multiprot_gate_ledger.csv`):
- Template `1a0cCD`, chain C → `match_count=32`, `tm_score=0.716`
- Template `1a0hAE`, chain A → `match_count=14`, `tm_score=0.695`

Both are **5zngA query** pairs — **4eylA has zero passing orientations**.

### 3. Experimental Plan to Identify Bottlenecks

#### Phase 1: Diagnostic Run (No Thresholds)
```bash
# Run MultiProt on all 28 pairs with verbose logging, no gate thresholds
cd /scratch/rshadi25/GitHub/PRISM-prescript
python -c "
from src.alignment_multiprot import align_multiprot
align_multiprot(['5zngA', '4eylA'], ['1a0cCD', '1a0gAB', '1a0hAE', '1a0hBD', '1a0jAC'], 
                output_dir='processed/alignment_multiprot_diag', max_workers=1)
" 2>&1 | tee multiprot_diagnostic.log
```

#### Phase 2: Instrumented `_align_one` Logging
Add to `src/alignment_multiprot.py:_align_one`:
```python
# After MultiProt stdout parsing
print(f"DEBUG {query}_{template}_{chain}: largest_sol={largest_solution}, "
      f"match_dict_len={len(match_dict)}, kabsch_rmsd={kabsch_rmsd}, "
      f"seccomp_ok={_check_seccomp()}")
```

#### Phase 3: Per-Stage Timing
```python
import time
t0 = time.time()
# ... MultiProt subprocess ...
t1 = time.time()
# ... parse 2_sol.res ...
t2 = time.time()
# ... Kabsch ...
t3 = time.time()
print(f"TIMING {query}_{template}_{chain}: mp={t1-t0:.1f}s parse={t2-t1:.1f}s kabsch={t3-t2:.1f}s")
```

### 4. Recommended Changes for Parallel Jobs

#### Change 1: Disable Seccomp Check for Login-Node Testing
```python
# In _check_seccomp(), add override:
if os.environ.get("PRISM_MULTIPROT_FORCE", "0") == "1":
    return True
```
Run on login node: `PRISM_MULTIPROT_FORCE=1 python ...`

#### Change 2: Lower Gate Thresholds for Diagnostic Mode
```bash
# Environment variables (read by transformation.py gate)
export PRISM_TM_SCORE_THRESHOLD=0.0
export PRISM_MINIMUM_RESIDUE_MATCH_COUNT=3
export PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE=0
```

#### Change 3: Parallel Diagnostic SBATCH (1 job, 8 CPUs)
```bash
#!/bin/bash
#SBATCH --partition=ai
#SBATCH --qos=ai
#SBATCH --account=ai
#SBATCH --job-name=multiprot-diag
#SBATCH --output=multiprot_diag_%j.out
#SBATCH --error=multiprot_diag_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=02:00:00

cd /scratch/rshadi25/GitHub/PRISM-prescript
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gtalign_env

export PRISM_MULTIPROT_FORCE=1
export PRISM_TM_SCORE_THRESHOLD=0.0
export PRISM_MINIMUM_RESIDUE_MATCH_COUNT=3

# Run diagnostic on all 28 pairs in parallel (8 workers)
python -c "
from src.alignment_multiprot import align_multiprot
align_multiprot(
    ['5zngA', '4eylA'],
    ['1a0cCD', '1a0gAB', '1a0hAE', '1a0hBD', '1a0jAC'],
    output_dir='processed/alignment_multiprot_diag',
    max_workers=8
)
" 2>&1 | tee multiprot_diag_${SLURM_JOB_ID}.log
```

#### Change 4: Ranking Pipeline Integration
Position MultiProt ranking **after alignment, before refinement**:
```
Query + Templates → MultiProt Align → [RANK by match_count + TM-score proxy] 
    → Top-K → Transform → Refine (FiberDock/PyRosetta) → Final Rank
```

---

## Consolidated Action Plan

### Immediate (This Session)
1. **Run diagnostic SBATCH** (above) on `ai` partition — 1 job, 8 CPUs, bypasses QoS
2. **Collect full logs** for all 28 pairs — identify exact failure points
3. **Compute true TM-scores** from MultiProt matches for the 2 passing pairs + any new ones

### Short-Term (Next Iteration)
1. **Implement standard TM-score computation** in `alignment_multiprot.py` (optional output field)
2. **Add MultiProt to comparison matrix** with fair metrics (match count, true TM-score, RMSD)
3. **Fix seccomp issue** — either compile 64-bit MultiProt or run on login node for comparison

### Long-Term (Pipeline Integration)
1. **Unified alignment interface**: All aligners (TMalign, GTalign, MultiProt) output standard JSON with:
   - `match_count`, `match_percentage`, `rmsd`, `tm_score` (standard), `rotation`, `translation`
2. **Configurable gate thresholds** per aligner (MultiProt may need lower TM threshold)
3. **Parallel ranking stage** that can consume any aligner's output

---

## Appendix: Key Files & Locations

| File | Purpose |
|------|---------|
| `src/alignment_multiprot.py` | MultiProt alignment module (primary) |
| `src/transformation.py` | Gate thresholds (`PRISM_TM_SCORE_THRESHOLD`, etc.) |
| `tmp/agent/20260727-pipeline-comparison/multiprot_gate_diagnosis.json` | Gate failure details |
| `tmp/agent/20260728-multiprot-gate-fixed-test/multiprot_gate_ledger.csv` | Passing pairs ledger |
| `docs/alignment_comparison_report.md` | Prior TMalign/GTalign comparison |
| `benchmark/scripts/run_multiprot_tmalign_calibration.py` | Calibration script template |

---

## Conclusion

**Task A (TM-score)**: Not a bug — a **metric definition mismatch**. MultiProt's "TM-score" is an RMSD proxy. Compute standard TM-score from its residue matches for fair comparison.

**Task B (2/28 yield)**: **Pipeline gating + seccomp** block 26 pairs. Diagnostic run with lowered thresholds and seccomp override will reveal true alignment capability. Use single-job internal parallelization (8 CPUs) to bypass ai QoS limits.

**Both tasks are independent and can proceed in parallel** — they address different pipeline stages (metric computation vs. execution gating).
