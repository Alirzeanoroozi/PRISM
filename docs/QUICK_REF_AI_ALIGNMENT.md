# Quick Reference: AI Structural Alignment Tools for PRISM

## 🎯 One-Line Decision

| If you need... | Use this tool |
|----------------|---------------|
| **Exact TMalign replacement, 1000x faster** | `gtalign` (GPU) |
| **Universal (protein+RNA+complex)** | `US-align` |
| **CPU-only, no GPU** | `foldseek --alignment-type 1` |
| **Search millions of structures** | `ssalign` → `gtalign` (two-stage) |
| **Binding site / motif alignment** | `localign` |
| **Interpretable local alignment** | `plasma` or `softalign` |

---

## 🚀 Quick Start Commands

### GTalign (Primary Recommendation)
```bash
# Install
conda install -c minmarg gtalign_gpu  # GPU (recommended)
conda install -c minmarg gtalign_mp   # CPU only

# Replace TMalign
gtalign query.pdb target.pdb --outfmt json --gpu
# Output: {"tm_score": 0.85, "rmsd": 1.2, "aligned_length": 150, "transform": [...]}
```

### US-align (Universal)
```bash
# Build from source
git clone https://github.com/pylelab/USalign && cd USalign && make

# Run (monomeric = identical to TMalign)
USalign query.pdb target.pdb -outfmt json
```

### Foldseek (CPU, Classical)
```bash
# Install
conda install -c bioconda foldseek

# Global alignment mode (replaces TMalign)
foldseek align query.pdb target.pdb out.tsv --alignment-type 1
# Or database search
foldseek easy-search query.pdb db out.tsv tmp --alignment-type 1
```

### SSAlign (Large Scale Prefilter)
```bash
# Install
git clone https://github.com/ISYSLAB-HUST/SSAlign.git
cd SSAlign && conda env create -f env.yml
# Download SaProt_650M_AF2.pt from HuggingFace to models/

# Two-stage pipeline
ssalign_prefilter query.pdb afdb50/ --top-k 100 --gpu 4 --out hits.json
saligner hits.json --out refined.tsv
# Or refine with GTalign
```

### LocAlign (Functional Sites)
```bash
# Install
git clone https://github.com/hagairavid18/LocAlign.git
cd LocAlign && conda env create -f inference_env.yaml

# Motif alignment
localign --source query.pdb --target target.pdb --motif "10,15,20,25" --out results/
# Outputs: transform, correspondence_map, pLRMSD, ChimeraX scripts
```

---

## 📊 Performance Comparison

| Tool | Speedup vs TMalign | Accuracy | GPU | Best For |
|------|-------------------|----------|-----|----------|
| **GTalign GPU** | **1,000x - 900,000x** | **Exact match** | ✅ | Primary replacement |
| **GTalign CPU** | 100-500x | Exact match | ❌ | No GPU env |
| **US-align** | ~100x | Exact match (protein) | ❌ | Multi-molecule |
| **Foldseek** | 100-1000x | Near (misses simple folds) | ❌ | CPU-only pipelines |
| **SSAlign** | 100-770x (prefilter) | TM-align comparable | ✅ | AFDB-scale search |
| **LocAlign** | N/A (local only) | Functional sites | ✅ | Binding motifs |

---

## ⚠️ Tools That CANNOT Replace TMalign

| Tool | Reason |
|------|--------|
| **PLASMA** | Substructure only; fails on full-length folds |
| **SoftAlign** | Needs downstream Kabsch solver; lower success rate |
| **DeepBLAST** | Sequence-only input (no 3D coordinates) |
| **DeepAlign** | Legacy, 7+ years old, CPU-only, too slow |
| **Caretta/FoldMason** | Multiple structure alignment (MSTA), not pairwise |
| **All Diffusion Models** | Generative only (RFdiffusion, EvoDiff, etc.) |

---

## 🔧 PRISM Integration Snippet

```python
# In prism.py - replace TMalign call
from src.structural_aligner import replace_tmalign_call

# OLD:
# tm_score, rmsd = run_tmalign(query.pdb, target.pdb)

# NEW (one line change):
result = replace_tmalign_call("query.pdb", "target.pdb", use_gpu=True)
tm_score, rmsd = result["tm_score"], result["rmsd"]
```

---

## 📦 Installation for PRISM Environment

```bash
# Minimal (GTalign only - recommended)
conda create -n prism_align python=3.10
conda activate prism_align
conda install -c minmarg gtalign_gpu

# Full stack
conda install -c bioconda foldseek
# SSAlign, LocAlign: see individual install above
```

---

## 📄 Full Documentation

See: `docs/AI_STRUCTURAL_ALIGNMENT_ALTERNATIVES.md`

- Complete tool specifications (11 tools)
- Installation details
- Input/output formats
- License compatibility
- Benchmarking framework
- Two-stage pipeline examples