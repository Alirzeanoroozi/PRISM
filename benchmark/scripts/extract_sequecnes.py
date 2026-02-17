#!/usr/bin/env python3
"""Extract chains and sequences from PDB files using Biopython."""

from pathlib import Path

from Bio.PDB import PDBParser, PPBuilder


def extract_sequences_from_pdb(pdb_path: str) -> dict[str, str]:
    """Extract chain sequences from a PDB file using Biopython.
    
    Returns dict mapping chain_id -> sequence (1-letter codes).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    ppb = PPBuilder()
    
    result = {}
    for model in structure:
        for chain in model:
            seq_parts = []
            for pp in ppb.build_peptides(chain):
                seq_parts.append(str(pp.get_sequence()))
            if seq_parts:
                result[chain.id] = "".join(seq_parts)
    return result


def main():
    sample_dir = Path(__file__).resolve().parent.parent / "sample"
    
    if not sample_dir.exists():
        print(f"Directory not found: {sample_dir}")
        return
    
    pdb_files = sorted(sample_dir.glob("*.pdb"))
    
    for pdb_path in pdb_files:
        print(f"\n{'='*60}")
        print(f"File: {pdb_path.name}")
        print("=" * 60)
        
        chains_seqs = extract_sequences_from_pdb(str(pdb_path))
        
        if not chains_seqs:
            print("  No chains found.")
            continue
            
        for chain_id, seq in sorted(chains_seqs.items()):
            chain_label = f"Chain {chain_id!r}" if chain_id != " " else "Chain ' '"
            print(f"\n  {chain_label} ({len(seq)} residues):")
            # Print sequence in 60-char lines
            for i in range(0, len(seq), 60):
                chunk = seq[i : i + 60]
                print(f"    {chunk}")


if __name__ == "__main__":
    main()
