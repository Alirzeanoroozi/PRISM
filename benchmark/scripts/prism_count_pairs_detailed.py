"""
Generate a detailed job-wise report for merged.txt pairs vs Rosetta outputs.

Inputs:
- merged.txt with submitted pairs.
- rosetta_output_1 joblist folders.

Output:
- Writes pair_analysis_report.txt with:
    * Overall hit/miss totals and percentages.
    * Top/bottom joblists by pairs found.
    * Full joblist ranking and missing pairs list.
    * Per-job listing of matched pairs.

Insight:
- Highlights which joblists contribute most/least to coverage.
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

# Extract pairs from rosetta output filenames with job tracking
rosetta_pairs = set()
job_pairs = defaultdict(set)  # Track pairs per job
rosetta_output_dir = '/scratch/rshadi25/GitHub/PRISM/benchmark/prism_processed/rosetta_output_1/'

# Pattern to extract pdb IDs from filename
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
                # Track which pairs are in merged.txt
                if pair in merged_pairs:
                    job_pairs[joblist].add(pair)

print(f"Total unique pairs with results in rosetta_output_1: {len(rosetta_pairs)}")

# Find matches
matches = merged_pairs & rosetta_pairs
missing_pairs = merged_pairs - rosetta_pairs

print(f"Pairs that have results: {len(matches)}")
print(f"Pairs without results: {len(missing_pairs)}")

# Job-wise analysis
print(f"\nTotal jobs analyzed: {len(job_pairs)}")

# Write results to file
output_file = '/scratch/rshadi25/GitHub/PRISM/benchmark/prism_processed/pair_analysis_report.txt'
with open(output_file, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("ROSETTA OUTPUT PAIR ANALYSIS REPORT\n")
    f.write("=" * 80 + "\n\n")
    
    f.write(f"Total pairs in merged.txt: {len(merged_pairs)}\n")
    f.write(f"Total unique pairs with results in rosetta_output_1: {len(rosetta_pairs)}\n")
    f.write(f"Pairs that have results: {len(matches)} ({len(matches)/len(merged_pairs)*100:.1f}%)\n")
    f.write(f"Pairs without results: {len(missing_pairs)} ({len(missing_pairs)/len(merged_pairs)*100:.1f}%)\n\n")
    
    f.write("=" * 80 + "\n")
    f.write("JOB-WISE ANALYSIS\n")
    f.write("=" * 80 + "\n\n")
    
    # Sort jobs by number of pairs (descending)
    sorted_jobs = sorted(job_pairs.items(), key=lambda x: len(x[1]), reverse=True)
    
    f.write(f"Total jobs with results: {len(sorted_jobs)}\n\n")
    
    # Top performers
    f.write("TOP 10 JOBS (Most pairs found):\n")
    f.write("-" * 80 + "\n")
    f.write(f"{'Job':<20} {'Pairs Found':<15} {'Percentage':<15}\n")
    f.write("-" * 80 + "\n")
    for job, pairs in sorted_jobs[:10]:
        percentage = len(pairs) / len(merged_pairs) * 100
        f.write(f"{job:<20} {len(pairs):<15} {percentage:.1f}%\n")
    
    f.write("\n")
    
    # Bottom performers
    f.write("BOTTOM 10 JOBS (Least pairs found):\n")
    f.write("-" * 80 + "\n")
    f.write(f"{'Job':<20} {'Pairs Found':<15} {'Percentage':<15}\n")
    f.write("-" * 80 + "\n")
    for job, pairs in sorted_jobs[-10:]:
        percentage = len(pairs) / len(merged_pairs) * 100
        f.write(f"{job:<20} {len(pairs):<15} {percentage:.1f}%\n")
    
    f.write("\n" + "=" * 80 + "\n")
    f.write("COMPLETE JOB LIST (sorted by pairs found, descending)\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"{'Job':<20} {'Pairs Found':<15} {'Percentage':<15}\n")
    f.write("-" * 80 + "\n")
    for job, pairs in sorted_jobs:
        percentage = len(pairs) / len(merged_pairs) * 100
        f.write(f"{job:<20} {len(pairs):<15} {percentage:.1f}%\n")
    
    f.write("\n" + "=" * 80 + "\n")
    f.write("PAIRS WITHOUT RESULTS\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Total missing: {len(missing_pairs)}\n\n")
    for pair in sorted(missing_pairs):
        f.write(f"{pair[0]:<15} {pair[1]:<15}\n")
    
    f.write("\n" + "=" * 80 + "\n")
    f.write("PAIRS WITH RESULTS (by job)\n")
    f.write("=" * 80 + "\n\n")
    for job, pairs in sorted_jobs:
        f.write(f"\n{job} ({len(pairs)} pairs):\n")
        f.write("-" * 40 + "\n")
        for pair in sorted(pairs):
            f.write(f"  {pair[0]:<15} {pair[1]:<15}\n")

print(f"\nReport saved to: {output_file}")
print("\n" + "=" * 80)
print("JOB-WISE SUMMARY")
print("=" * 80)
print(f"{'Job':<20} {'Pairs Found':<15} {'Percentage':<15}")
print("-" * 80)

# Show top and bottom 5 jobs
sorted_jobs = sorted(job_pairs.items(), key=lambda x: len(x[1]), reverse=True)
print("\nTOP 5 JOBS:")
for job, pairs in sorted_jobs[:5]:
    percentage = len(pairs) / len(merged_pairs) * 100
    print(f"{job:<20} {len(pairs):<15} {percentage:.1f}%")

print("\nBOTTOM 5 JOBS:")
for job, pairs in sorted_jobs[-5:]:
    percentage = len(pairs) / len(merged_pairs) * 100
    print(f"{job:<20} {len(pairs):<15} {percentage:.1f}%")

