"""
Convert pair analysis into CSV files for spreadsheet analysis.

Inputs:
- merged.txt with submitted pairs.
- rosetta_output_1 joblist folders.

Outputs (in prism_processed):
- overall_summary.csv: coverage totals and percentages.
- job_wise_summary.csv: pairs found per joblist and percent of total.
- pairs_with_results.csv: matched pairs and the joblist(s) containing them.
- pairs_without_results.csv: missing pairs to investigate or re-run.

Insight:
- CSVs are sorted for easy filtering, with pair order normalized.
"""

import os
import re
import csv
from collections import defaultdict

# Read merged.txt
merged_pairs = set()
merged_pairs_list = []
with open('/scratch/rshadi25/GitHub/PRISM/benchmark/prism_processed/LISTS/merged.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            pdb1, pdb2 = parts
            merged_pairs_list.append((pdb1, pdb2))
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

# Find matches
matches = merged_pairs & rosetta_pairs
missing_pairs = merged_pairs - rosetta_pairs

# Generate CSV files
output_dir = '/scratch/rshadi25/GitHub/PRISM/benchmark/prism_processed/'

# 1. Job-wise summary CSV
csv_job_summary = os.path.join(output_dir, 'job_wise_summary.csv')
with open(csv_job_summary, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Job', 'Pairs_Found', 'Percentage', 'Status'])
    
    sorted_jobs = sorted(job_pairs.items(), key=lambda x: len(x[1]), reverse=True)
    for job, pairs in sorted_jobs:
        percentage = len(pairs) / len(merged_pairs) * 100
        status = 'High' if len(pairs) >= 8 else 'Low' if len(pairs) <= 3 else 'Medium'
        writer.writerow([job, len(pairs), f"{percentage:.2f}", status])

print(f"Created: {csv_job_summary}")

# 2. Pairs with results CSV
csv_pairs_with_results = os.path.join(output_dir, 'pairs_with_results.csv')
with open(csv_pairs_with_results, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['PDB1', 'PDB2', 'Jobs'])
    
    # Group by pair and list all jobs
    pair_jobs = defaultdict(list)
    for job, pairs in job_pairs.items():
        for pair in pairs:
            pair_jobs[pair].append(job)
    
    for pair in sorted(pair_jobs.keys()):
        jobs = ';'.join(sorted(pair_jobs[pair]))
        writer.writerow([pair[0], pair[1], jobs])

print(f"Created: {csv_pairs_with_results}")

# 3. Missing pairs CSV
csv_missing_pairs = os.path.join(output_dir, 'pairs_without_results.csv')
with open(csv_missing_pairs, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['PDB1', 'PDB2'])
    
    for pair in sorted(missing_pairs):
        writer.writerow([pair[0], pair[1]])

print(f"Created: {csv_missing_pairs}")

# 4. Overall summary CSV
csv_overall = os.path.join(output_dir, 'overall_summary.csv')
with open(csv_overall, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Metric', 'Value', 'Percentage'])
    writer.writerow(['Total Pairs in merged.txt', len(merged_pairs), '100.00'])
    writer.writerow(['Pairs with Results', len(matches), f"{len(matches)/len(merged_pairs)*100:.2f}"])
    writer.writerow(['Pairs without Results', len(missing_pairs), f"{len(missing_pairs)/len(merged_pairs)*100:.2f}"])
    writer.writerow(['Total Unique Pairs in Output', len(rosetta_pairs), '-'])
    writer.writerow(['Total Jobs Analyzed', len(job_pairs), '-'])

print(f"Created: {csv_overall}")

print("\n=== CSV Generation Complete ===")
print(f"Total CSV files created: 4")
print(f"Location: {output_dir}")
