#!/usr/bin/env python3
"""
Patch transformation.py to use separate gate for MultiProt with lowered thresholds.
This should be applied after the diagnostic run to enable MultiProt candidates
to pass transformation filtering.
"""
import re

with open("src/transformation.py", "r") as f:
    content = f.read()

# Find the alignment_score_passes function and add MultiProt-specific logic
# Look for the existing special-case and enhance it

old_pattern = r'if aligner == "multiprot":\s+# MultiProt uses proxy TM-score\s+return tm_score >= 0\.5 and match_count >= 10'

new_block = '''if aligner == "multiprot":
    # MultiProt uses proxy TM-score (1 - RMSD/10), not standard TM-score
    # Lowered thresholds based on diagnostic run results
    return tm_score >= 0.3 and match_count >= 10 and match_pct >= 30'''

content = re.sub(old_pattern, new_block, content)

with open("src/transformation.py", "w") as f:
    f.write(content)

print("Updated transformation.py MultiProt gate")
