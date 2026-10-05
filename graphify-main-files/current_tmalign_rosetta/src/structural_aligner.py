#!/usr/bin/env python3
"""
PRISM Pipeline - Structural Alignment Module
AI/Deep Learning Alternatives to TMalign

This module provides a unified interface for multiple structural alignment tools
that can replace or complement TMalign in the PRISM pipeline.

Tools supported:
- GTalign (primary recommendation - exact TMalign replacement)
- US-align (universal: proteins + nucleic acids + complexes)
- Foldseek (CPU-only, SIMD vectorized)
- SSAlign (large-scale prefilter)
- LocAlign (functional site / motif alignment)
- PLASMA (local substructure alignment)
- SoftAlign (optimal transport + Kabsch)
"""

import subprocess
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class AlignmentTool(Enum):
    GTALIGN_GPU = "gtalign_gpu"
    GTALIGN_CPU = "gtalign_cpu"
    USALIGN = "usalign"
    FOLDSEEK = "foldseek"
    SSALIGN = "ssalign"
    LOCALIGN = "localign"
    PLASMA = "plasma"
    SOFTALIGN = "softalign"


@dataclass
class AlignmentResult:
    """Standardized alignment result across all tools"""
    tool: str
    query: str
    target: str
    tm_score: float
    rmsd: float
    aligned_length: int
    transformation_matrix: Optional[List[List[float]]] = None
    correspondence_map: Optional[Dict] = None
    confidence_score: Optional[float] = None
    raw_output: Optional[str] = None
    success: bool = True
    error: Optional[str] = None


class StructuralAligner:
    """
    Unified interface for structural alignment tools in PRISM pipeline.
    
    Usage:
        aligner = StructuralAligner()
        
        # Simple pairwise (replaces TMalign)
        result = aligner.align_pair("query.pdb", "target.pdb", tool=AlignmentTool.GTALIGN_GPU)
        
        # Large database search (two-stage)
        hits = aligner.prefilter_database("query.pdb", "afdb50/", top_k=100)
        results = aligner.refine_hits("query.pdb", hits)
        
        # Functional site alignment
        result = aligner.align_motif("query.pdb", "target.pdb", motif_residues=[10, 15, 20])
    """
    
    def __init__(self, 
                 gtalign_path: str = "gtalign",
                 usalign_path: str = "USalign",
                 foldseek_path: str = "foldseek",
                 ssalign_prefilter_path: str = "ssalign_prefilter",
                 saligner_path: str = "saligner",
                 localign_path: str = "localign"):
        self.tools = {
            AlignmentTool.GTALIGN_GPU: gtalign_path,
            AlignmentTool.GTALIGN_CPU: gtalign_path,
            AlignmentTool.USALIGN: usalign_path,
            AlignmentTool.FOLDSEEK: foldseek_path,
            AlignmentTool.SSALIGN: ssalign_prefilter_path,
            AlignmentTool.LOCALIGN: localign_path,
        }
        self._verify_tools()
    
    def _verify_tools(self):
        """Check which tools are available"""
        self.available = {}
        for tool, path in self.tools.items():
            try:
                result = subprocess.run([path, "--help"], 
                                      capture_output=True, timeout=5)
                self.available[tool] = result.returncode == 0
                if self.available[tool]:
                    logger.info(f"✓ {tool.value} available at {path}")
            except (FileNotFoundError, subprocess.TimeoutExpired):
                self.available[tool] = False
                logger.warning(f"✗ {tool.value} not found at {path}")
    
    def align_pair(self, 
                   query: str, 
                   target: str, 
                   tool: AlignmentTool = AlignmentTool.GTALIGN_GPU,
                   **kwargs) -> AlignmentResult:
        """
        Perform pairwise structural alignment (direct TMalign replacement).
        
        Args:
            query: Path to query PDB/mmCIF
            target: Path to target PDB/mmCIF
            tool: Alignment tool to use
            **kwargs: Tool-specific options
            
        Returns:
            AlignmentResult with standardized fields
        """
        if tool == AlignmentTool.GTALIGN_GPU:
            return self._run_gtalign(query, target, gpu=True, **kwargs)
        elif tool == AlignmentTool.GTALIGN_CPU:
            return self._run_gtalign(query, target, gpu=False, **kwargs)
        elif tool == AlignmentTool.USALIGN:
            return self._run_usalign(query, target, **kwargs)
        elif tool == AlignmentTool.FOLDSEEK:
            return self._run_foldseek_pair(query, target, **kwargs)
        else:
            return AlignmentResult(
                tool=tool.value, query=query, target=target,
                tm_score=0.0, rmsd=0.0, aligned_length=0,
                success=False, error=f"Tool {tool.value} not supported for pairwise"
            )
    
    def _run_gtalign(self, query: str, target: str, gpu: bool = True, **kwargs) -> AlignmentResult:
        """Run GTalign (exact TMalign replacement)"""
        cmd = [self.tools[AlignmentTool.GTALIGN_GPU if gpu else AlignmentTool.GTALIGN_CPU]]
        cmd.extend([query, target])
        cmd.extend(["--outfmt", "json"])
        if gpu:
            cmd.append("--gpu")
        else:
            cmd.append("--cpu")
        
        # Add any extra args
        for k, v in kwargs.items():
            cmd.extend([f"--{k}", str(v)])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode != 0:
                return AlignmentResult(
                    tool="gtalign", query=query, target=target,
                    tm_score=0.0, rmsd=0.0, aligned_length=0,
                    success=False, error=result.stderr, raw_output=result.stdout
                )
            
            # Parse JSON output
            data = json.loads(result.stdout)
            return AlignmentResult(
                tool="gtalign_gpu" if gpu else "gtalign_cpu",
                query=query, target=target,
                tm_score=data.get("tm_score", 0.0),
                rmsd=data.get("rmsd", 0.0),
                aligned_length=data.get("aligned_length", 0),
                transformation_matrix=data.get("transform"),
                raw_output=result.stdout,
                success=True
            )
        except json.JSONDecodeError as e:
            return AlignmentResult(
                tool="gtalign", query=query, target=target,
                tm_score=0.0, rmsd=0.0, aligned_length=0,
                success=False, error=f"JSON parse error: {e}", raw_output=result.stdout
            )
        except subprocess.TimeoutExpired:
            return AlignmentResult(
                tool="gtalign", query=query, target=target,
                tm_score=0.0, rmsd=0.0, aligned_length=0,
                success=False, error="Timeout"
            )
    
    def _run_usalign(self, query: str, target: str, **kwargs) -> AlignmentResult:
        """Run US-align (universal alignment)"""
        cmd = [self.tools[AlignmentTool.USALIGN], query, target, "-outfmt", "json"]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode != 0:
                return AlignmentResult(
                    tool="usalign", query=query, target=target,
                    tm_score=0.0, rmsd=0.0, aligned_length=0,
                    success=False, error=result.stderr
                )
            
            data = json.loads(result.stdout)
            return AlignmentResult(
                tool="usalign", query=query, target=target,
                tm_score=data.get("tm_score", 0.0),
                rmsd=data.get("rmsd", 0.0),
                aligned_length=data.get("aligned_length", 0),
                transformation_matrix=data.get("transform"),
                raw_output=result.stdout,
                success=True
            )
        except Exception as e:
            return AlignmentResult(
                tool="usalign", query=query, target=target,
                tm_score=0.0, rmsd=0.0, aligned_length=0,
                success=False, error=str(e)
            )
    
    def _run_foldseek_pair(self, query: str, target: str, **kwargs) -> AlignmentResult:
        """Run Foldseek in pairwise global alignment mode"""
        cmd = [self.tools[AlignmentTool.FOLDSEEK], "align", query, target, "stdout"]
        cmd.extend(["--alignment-type", "1", "--format-output", "query,target,tmscore,rmsd,alnlen"])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode != 0:
                return AlignmentResult(
                    tool="foldseek", query=query, target=target,
                    tm_score=0.0, rmsd=0.0, aligned_length=0,
                    success=False, error=result.stderr
                )
            
            # Parse TSV output
            parts = result.stdout.strip().split("\t")
            if len(parts) >= 5:
                return AlignmentResult(
                    tool="foldseek", query=query, target=target,
                    tm_score=float(parts[2]),
                    rmsd=float(parts[3]),
                    aligned_length=int(parts[4]),
                    raw_output=result.stdout,
                    success=True
                )
            else:
                return AlignmentResult(
                    tool="foldseek", query=query, target=target,
                    tm_score=0.0, rmsd=0.0, aligned_length=0,
                    success=False, error="Unexpected output format"
                )
        except Exception as e:
            return AlignmentResult(
                tool="foldseek", query=query, target=target,
                tm_score=0.0, rmsd=0.0, aligned_length=0,
                success=False, error=str(e)
            )
    
    def prefilter_database(self, 
                          query: str, 
                          database_dir: str, 
                          top_k: int = 100,
                          gpu_count: int = 1) -> List[Dict]:
        """
        Stage 1: Fast prefilter using SSAlign for large databases.
        
        Returns list of hit dictionaries with estimated scores.
        """
        if not self.available.get(AlignmentTool.SSALIGN):
            raise RuntimeError("SSAlign not available. Install from https://github.com/ISYSLAB-HUST/SSAlign")
        
        cmd = [
            self.tools[AlignmentTool.SSALIGN],
            query, database_dir,
            "--top-k", str(top_k),
            "--gpu", str(gpu_count),
            "--outfmt", "json"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if result.returncode != 0:
                logger.error(f"SSAlign prefilter failed: {result.stderr}")
                return []
            
            hits = json.loads(result.stdout)
            return hits
        except Exception as e:
            logger.error(f"SSAlign prefilter error: {e}")
            return []
    
    def refine_hits(self, 
                   query: str, 
                   hits: List[Dict],
                   tool: AlignmentTool = AlignmentTool.GTALIGN_GPU) -> List[AlignmentResult]:
        """
        Stage 2: Precise refinement of prefilter hits using GTalign/US-align.
        """
        results = []
        for hit in hits:
            target_path = hit.get("path") or hit.get("target_path")
            if target_path and os.path.exists(target_path):
                result = self.align_pair(query, target_path, tool=tool)
                result.confidence_score = hit.get("ss_score") or hit.get("estimated_tm")
                results.append(result)
        return results
    
    def align_motif(self, 
                   query: str, 
                   target: str, 
                   motif_residues: List[int],
                   **kwargs) -> AlignmentResult:
        """
        Functional site / motif alignment using LocAlign.
        
        Args:
            query: Query structure PDB
            target: Target structure PDB
            motif_residues: List of residue indices defining the motif/binding site
            **kwargs: Additional LocAlign options
        """
        if not self.available.get(AlignmentTool.LOCALIGN):
            raise RuntimeError("LocAlign not available. Install from https://github.com/hagairavid18/LocAlign")
        
        # Create motif mask file
        motif_file = f"/tmp/motif_{os.path.basename(query)}.txt"
        with open(motif_file, "w") as f:
            f.write(",".join(map(str, motif_residues)))
        
        cmd = [
            self.tools[AlignmentTool.LOCALIGN],
            "--source", query,
            "--target", target,
            "--motif", motif_file,
            "--outfmt", "json"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            os.remove(motif_file)
            
            if result.returncode != 0:
                return AlignmentResult(
                    tool="localign", query=query, target=target,
                    tm_score=0.0, rmsd=0.0, aligned_length=0,
                    success=False, error=result.stderr
                )
            
            data = json.loads(result.stdout)
            return AlignmentResult(
                tool="localign", query=query, target=target,
                tm_score=data.get("tm_score", 0.0),
                rmsd=data.get("rmsd", 0.0),
                aligned_length=data.get("aligned_length", 0),
                transformation_matrix=data.get("transform"),
                correspondence_map=data.get("correspondence"),
                confidence_score=data.get("plrmsd"),
                raw_output=result.stdout,
                success=True
            )
        except Exception as e:
            if os.path.exists(motif_file):
                os.remove(motif_file)
            return AlignmentResult(
                tool="localign", query=query, target=target,
                tm_score=0.0, rmsd=0.0, aligned_length=0,
                success=False, error=str(e)
            )
    
    def benchmark_tools(self, 
                       test_pairs: List[Tuple[str, str]],
                       tools: List[AlignmentTool] = None) -> Dict:
        """
        Benchmark multiple tools on test pairs for PRISM validation.
        """
        if tools is None:
            tools = [AlignmentTool.GTALIGN_GPU, AlignmentTool.USALIGN, AlignmentTool.FOLDSEEK]
        
        results = {tool.value: [] for tool in tools}
        
        for query, target in test_pairs:
            for tool in tools:
                if self.available.get(tool):
                    result = self.align_pair(query, target, tool=tool)
                    results[tool.value].append({
                        "query": query,
                        "target": target,
                        "tm_score": result.tm_score,
                        "rmsd": result.rmsd,
                        "time": result.raw_output,  # Would need timing wrapper
                        "success": result.success
                    })
        
        return results


# Convenience functions for PRISM integration
def replace_tmalign_call(query_pdb: str, target_pdb: str, use_gpu: bool = True) -> Dict:
    """
    Drop-in replacement for TMalign call in PRISM pipeline.
    
    Returns dict with same keys as TMalign parser would produce.
    """
    aligner = StructuralAligner()
    tool = AlignmentTool.GTALIGN_GPU if use_gpu else AlignmentTool.GTALIGN_CPU
    result = aligner.align_pair(query_pdb, target_pdb, tool=tool)
    
    if result.success:
        return {
            "tm_score": result.tm_score,
            "rmsd": result.rmsd,
            "aligned_length": result.aligned_length,
            "transform_matrix": result.transformation_matrix,
            "tool": "gtalign"
        }
    else:
        raise RuntimeError(f"Alignment failed: {result.error}")


def two_stage_large_scale_search(query_pdb: str, 
                                 database_path: str, 
                                 top_k: int = 100,
                                 refine_top: int = 10) -> List[Dict]:
    """
    Two-stage search for large databases (AFDB-scale).
    
    Stage 1: SSAlign prefilter (100-770x speedup)
    Stage 2: GTalign precise refinement
    """
    aligner = StructuralAligner()
    
    # Stage 1: Fast prefilter
    logger.info(f"Stage 1: SSAlign prefilter on {database_path}")
    hits = aligner.prefilter_database(query_pdb, database_path, top_k=top_k)
    
    # Stage 2: Precise refinement on top hits
    logger.info(f"Stage 2: GTalign refinement on top {refine_top} hits")
    top_hits = hits[:refine_top]
    results = aligner.refine_hits(query_pdb, top_hits, tool=AlignmentTool.GTALIGN_GPU)
    
    return [
        {
            "target": r.target,
            "tm_score": r.tm_score,
            "rmsd": r.rmsd,
            "aligned_length": r.aligned_length,
            "prefilter_score": r.confidence_score,
            "transform_matrix": r.transformation_matrix
        }
        for r in results if r.success
    ]


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) == 3:
        query, target = sys.argv[1], sys.argv[2]
        result = replace_tmalign_call(query, target)
        print(json.dumps(result, indent=2))
    else:
        print("Usage: python structural_aligner.py query.pdb target.pdb")
        print("\nFor PRISM integration, import this module and use:")
        print("  from structural_aligner import replace_tmalign_call, two_stage_large_scale_search")