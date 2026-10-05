"""
Count how many unique pairs from merged.txt have Rosetta outputs.

Inputs:
- merged.txt with submitted pairs.
- rosetta_output_1 joblist folders.

Output:
- Prints total unique pairs, matched count, and missing count to stdout.

Insight:
- Hit/miss coverage is computed with order-independent pair matching.
"""

import os
import re
from collections import defaultdict

# Read merged.txt
merged_pairs = set()
with open('/scratch/rshadi25/GitHub/PRISM/benchmark/prism_processed/LISTS/merged.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            pdb1, pdb2 = parts
            # Store as a sorted tuple to handle both orders
            pair = tuple(sorted([pdb1, pdb2]))
            merged_pairs.add(pair)

print(f"Total pairs in merged.txt: {len(merged_pairs)}")

# Extract pairs from rosetta output filenames
rosetta_pairs = set()
rosetta_output_dir = '/scratch/rshadi25/GitHub/PRISM/benchmark/prism_processed/rosetta_output_1/'

# Pattern to extract pdb IDs from filename
# Format: something_pdb1_0_pdb2_0.rosetta_...
pattern = r'_([A-Za-z0-9]+)_0_([A-Za-z0-9]+)_0\.rosetta'

for joblist in os.listdir(rosetta_output_dir):
    joblist_path = os.path.join(rosetta_output_dir, joblist)
    if os.path.isdir(joblist_path):
        for filename in os.listdir(joblist_path):
            match = re.search(pattern, filename)
            if match:
                pdb1, pdb2 = match.groups()
                # Store as sorted tuple to handle both orders
                pair = tuple(sorted([pdb1, pdb2]))
                rosetta_pairs.add(pair)

print(f"Total pairs with results in rosetta_output_1: {len(rosetta_pairs)}")

# Find matches
matches = merged_pairs & rosetta_pairs
print(f"Pairs that have results: {len(matches)}")
print(f"Pairs without results: {len(merged_pairs) - len(matches)}")

# Show some statistics
missing_pairs = merged_pairs - rosetta_pairs
if missing_pairs:
    print(f"\nFirst 10 pairs WITHOUT results:")
    for pair in sorted(list(missing_pairs))[:10]:
        print(f"  {pair[0]} {pair[1]}")

