# PRISM Pipeline vs `nprot.2011.367` Protocol

This document summarizes a code-level audit of the repository pipeline against the protocol described in `references/nprot.2011.367.pdf` (PRISM protocol paper), with explicit notes on implemented, changed, missing, and likely broken steps.

## Scope and context

- Paper baseline: PRISM protocol workflow (template preparation, target preprocessing, structural matching, transformation filtering, flexible refinement/ranking).
- Repository pipeline entrypoint: `prism.py`.
- Important local substitutions already present in this repo:
  - MultiProt replaced by TM-align (default) or GTalign (optional).
  - FiberDock replaced by Rosetta docking refinement.

## High-level verdict

- The code implements the overall PRISM-like workflow shape.
- Several core steps are present (interface extraction, hotspot computation, alignment, transform generation, clash filtering, refinement stage).
- However, the current implementation is **not fully faithful** to the protocol as written.
- Most impactful gaps are in transformation filtering and in Rosetta stage correctness.

---

## 1. Paper workflow steps and code mapping

### 1.1 Input handling / target acquisition

Paper intent:
- Accept target proteins/chains and prepare structures for analysis.

Code:
- Reads receptor/ligand list from CSV and downloads PDB files.
- Files:
  - `src/pdb_download.py:33`
  - `src/pdb_download.py:39`
  - `src/pdb_download.py:41`
  - `src/pdb_download.py:49`
- Notes:
  - Input is constrained to 5-character IDs (PDB+chain assumption), not arbitrary multi-chain targets.

Status: **Implemented (simplified assumptions).**

### 1.2 Template preprocessing and filtering

Paper intent:
- Template set curation (valid complexes, structure quality/format constraints).

Code:
- Optional template analysis and filtering stage.
- Keeps templates with `num_chains > 1` and `structure_type == single`.
- Files:
  - `src/analyse_pdbs.py:61`
  - `src/analyse_pdbs.py:85`
  - `src/analyse_pdbs.py:88`
  - `prism.py:21`

Status: **Implemented (coarse filtering; not all protocol curation rules are evident).**

### 1.3 Interface residue extraction for templates

Paper intent:
- Build template interfaces for each partner chain.

Code:
- Computes interacting residues using heavy-atom vdW overlap criterion (`vdW1 + vdW2 + 0.5`), then adds nearby residues by Cα proximity.
- Writes chain-specific Cα interface PDBs and interface residue lists.
- Files:
  - `src/interface.py:37`
  - `src/interface.py:73`
  - `src/interface.py:88`
  - `src/interface.py:106`
  - `src/interface.py:145`

Status: **Implemented.**

### 1.4 Hotspot prediction for template interfaces

Paper intent:
- HotPoint-like hotspot identification using burial + contact potential thresholds.

Code:
- Uses NACCESS relative ASA and pair-potential-based contact scoring.
- Thresholds match HotPoint-style criteria in README/comments:
  - relative ASA <= 20%
  - |contact potential| >= 18
- Files:
  - `src/hotspot.py:11`
  - `src/hotspot.py:12`
  - `src/hotspot.py:24`
  - `src/hotspot.py:30`
  - `src/naccess_utils.py:10`

Status: **Implemented (local reimplementation, not web-service call).**

### 1.5 Template contact map generation

Paper intent:
- Keep interface residue-residue contacts for later filtering/scoring.

Code:
- Generates residue contact pairs using minimum heavy-atom distance <= 4.0 Å.
- File:
  - `src/contact.py:7`
  - `src/contact.py:12`

Status: **Implemented as artifact generation.**

### 1.6 Target surface extraction

Paper intent:
- Surface residue extraction using ASA and neighborhood expansion.

Code:
- Uses NACCESS on targets, keeps residues with relative ASA > 15.
- Adds residues by Cα proximity threshold `SCFFTHRESHOLD = 1.4`.
- Writes `.asa.pdb` Cα-only files for alignment.
- Files:
  - `src/surface_extract.py:8`
  - `src/surface_extract.py:9`
  - `src/surface_extract.py:20`
  - `src/surface_extract.py:56`

Status: **Implemented, but likely parameter mismatch (see gaps).**

### 1.7 Structural alignment stage

Paper intent:
- MultiProt structural matching of target surfaces to template interfaces; keep meaningful matches.

Code:
- Default backend: TM-align, parsed into transform + match dictionary JSON.
- Optional backend: GTalign with PRISM-compatible JSON output.
- Files:
  - `src/alignment.py:7`
  - `src/alignment.py:14`
  - `src/alignment.py:25`
  - `src/alignment.py:95`
  - `src/alignment_gtalign.py:12`
  - `prism.py:42`

Status: **Implemented with replacement algorithm (not MultiProt-equivalent).**

### 1.8 Transformation filtering / candidate assembly

Paper intent:
- Apply paper thresholds (match quality, hotspot support, contacts, clash checks) before refinement.

Code:
- Applies `match_count`, `tm_score`, match-percentage checks, and CA-clash filtering.
- Supports two orientations (chain mapping swapped).
- Files:
  - `src/transformation.py:54`
  - `src/transformation.py:72`
  - `src/transformation.py:80`
  - `src/transformation.py:93`
  - `src/transformation.py:108`
  - `src/transformation.py:156`

Status: **Partially implemented; major protocol filters missing/stubbed (see gaps).**

### 1.9 Refinement and final ranking

Paper intent:
- FiberDock flexible refinement + energy-based ranking.

Code:
- Replaced by Rosetta prepack + local docking refine, then score extraction and thresholding.
- Files:
  - `src/rosetta_refinement.py:14`
  - `src/rosetta_refinement.py:31`
  - `src/rosetta_refinement.py:35`
  - `src/rosetta_refinement.py:58`
  - `prism.py:67`

Status: **Implemented as a replacement stage, not protocol-equivalent.**

---

## 2. Missing or diverging protocol behavior

### 2.1 Hotspot filtering is effectively disabled in transformation stage

- `hotspot_analysis()` returns `True` unconditionally.
- File:
  - `src/transformation.py:51`
- Impact:
  - Any hotspot-based acceptance criterion from the protocol is currently bypassed.

### 2.2 Contact-based filtering threshold is not used

- Contact-related constants exist, but template contact files are not consumed in candidate filtering.
- Files:
  - `src/transformation.py:20`
  - `src/contact.py:12`
- Impact:
  - Matching-contact criterion from PRISM filtering is not enforced.

### 2.3 MultiProt replacement is a behavioral change, not a drop-in substitution

- TM-align/GTalign impose different matching semantics versus MultiProt.
- PRISM paper emphasizes MultiProt’s structural matching behavior and multiple transforms.
- Files:
  - `src/alignment.py:14`
  - `src/alignment_gtalign.py:95`
- Impact:
  - Candidate generation characteristics may differ significantly from protocol outputs.

### 2.4 Multiple-solution handling appears compressed to one JSON per pair

- Downstream expects one alignment JSON per `(query, template, chain)`.
- Files:
  - `src/alignment.py:95`
  - `src/alignment_gtalign.py:203`
- Impact:
  - The paper’s multi-solution matching/ranking behavior is likely underrepresented.

### 2.5 Target surface neighborhood threshold is suspiciously strict

- `SCFFTHRESHOLD = 1.4` Å in `src/surface_extract.py:9`.
- Impact:
  - This is unusually small for Cα-neighborhood expansion and likely under-selects nearby residues.

### 2.6 Protocol-level sequence/homology handling not clearly present

- No explicit FASTA/homology grouping step found in active pipeline.
- Impact:
  - Some protocol assumptions around homologous chain handling may be missing.

---

## 3. Likely bugs in current implementation

### 3.1 Rosetta partner chain IDs likely parsed incorrectly

- Code derives partner chains from fixed character positions in full file paths:
  - `partner_chains = f"{passed0[4]}_{passed1[4]}"`
- File:
  - `src/rosetta_refinement.py:30`
- Impact:
  - Rosetta may receive wrong chain definitions and produce invalid or failed docking runs.

### 3.2 Same fragile indexing used when splitting Rosetta output by chain

- Uses `passed0[4]` / `passed1[4]` to identify chain letters in output parsing.
- Files:
  - `src/rosetta_refinement.py:63`
  - `src/rosetta_refinement.py:65`
- Impact:
  - Chain extraction from docked structures can be incorrect.

### 3.3 Function signature mismatch for `get_contacts()` call in Rosetta stage

- `get_contacts()` is defined as single-argument (`template`) function.
- Rosetta stage calls it with 4 arguments.
- Files:
  - `src/contact.py:12`
  - `src/rosetta_refinement.py:74`
- Impact:
  - This call path should raise an exception (currently caught and logged), so intended post-refinement contact export does not execute.

### 3.4 Global mutable `passed_pairs` can accumulate across calls

- `passed_pairs` defined at module scope and appended in `transformer()` pipeline calls.
- Files:
  - `src/transformation.py:22`
  - `src/transformation.py:27`
- Impact:
  - Repeated runs within same process can leak state.

---

## 4. What is correctly aligned with protocol spirit

- NACCESS-driven ASA use is present for both template hotspots and target surfaces.
- HotPoint-like hotspot criteria are implemented (buried + contact potential).
- Structural alignment-derived transforms are used to place full targets onto template interface frames.
- Basic steric clash filtering is applied before refinement.
- A refinement/scoring stage exists (though with Rosetta, not FiberDock).

---

## 5. Practical conclusion

The repository has a functioning PRISM-style end-to-end pipeline, but currently it should be considered a **modified variant** of the paper method rather than a faithful reimplementation.

Most critical to claim paper-level parity:
1. Implement real hotspot-based transformation filtering (replace stub at `src/transformation.py:51`).
2. Implement contact-based match filtering using generated template contacts.
3. Fix Rosetta chain parsing and contact-call signature mismatch.
4. Revisit surface neighborhood threshold (`SCFFTHRESHOLD`) and validate against protocol definitions.
5. Document and, if possible, compensate for MultiProt -> TM-align/GTalign behavioral differences.

---

## 6. Audited pipeline entry sequence

- `prism.py:14` PDB download
- `prism.py:22` template filtering (optional)
- `prism.py:29` template generation (interface/hotspot/contact)
- `prism.py:38` target surface extraction
- `prism.py:42` structural alignment (TM-align / GTalign)
- `prism.py:60` transformation filtering + clash checks
- `prism.py:67` Rosetta refinement

