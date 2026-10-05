# Quickstart: Run Identity and Manifest Validation

**Phase**: 01-run-identity-and-manifest  
**Audience**: Implementers and reviewers  
**Goal**: Runnable scenarios that validate contracts, ledger, and gate against the spec.

---

## Prerequisites

```bash
cd /scratch/rshadi25/GitHub/PRISM-prescript
conda activate gtalign_env
pip install -e .  # if setup.py/pyproject.toml exists, else skip
```

All commands run from repository root unless noted.

---

## Scenario 1: Declare a Contract (US1)

**Goal**: Create an immutable declared contract from user selectors + template inventory + parameters.

```bash
# 1. Prepare a minimal pair list (1 row for smoke test)
cat > /tmp/test_pairs.txt <<'EOF'
5zngA_1buhAB_A
EOF

# 2. Declare contract via CLI (to be implemented: prism.py --declare-contract)
python -m src.run_identity declare-contract \
  --pairs /tmp/test_pairs.txt \
  --aligner tmalign \
  --refiner external_rosetta \
  --stages input,surface,alignment,transformation,refinement,comparison \
  --out /tmp/declared-contract.json

# 3. Verify canonical JSON + hash
cat /tmp/declared-contract.json
python -c "
import json, hashlib
with open('/tmp/declared-contract.json') as f:
    data = json.load(f)
canon = json.dumps(data, sort_keys=True, separators=(',', ':')).encode('utf-8')
h = hashlib.sha256(canon).hexdigest()
print('declared_contract_hash:', h)
"
```

**Expected**: JSON validates against `declared-contract.schema.json`; hash is 64 lowercase hex chars.

---

## Scenario 2: Initialize Run Identity (US2)

**Goal**: Create run_id and execution attempt manifest from declared contract.

```bash
# From Scenario 1's declared contract hash
DECLARED_HASH=$(python -c "
import json, hashlib
with open('/tmp/declared-contract.json') as f:
    data = json.load(f)
canon = json.dumps(data, sort_keys=True, separators=(',', ':')).encode('utf-8')
print(hashlib.sha256(canon).hexdigest())
")

# Initialize run identity
python -m src.run_identity init-run \
  --declared-contract-hash $DECLARED_HASH \
  --declared-contract-path /tmp/declared-contract.json \
  --out /tmp/run-manifest.json

# Verify run_id format
python -c "
import json, re
with open('/tmp/run-manifest.json') as f:
    m = json.load(f)
rid = m['run_id']
print('run_id:', rid)
print('matches pattern:', bool(re.match(r'^prism-\d{8}-\d{6}-\d+-[a-f0-9]{8}$', rid)))
print('declared_contract_hash match:', m['declared_contract_hash'] == '$DECLARED_HASH')
print('status:', m['status'])
"
```

**Expected**: `run_id` matches `prism-YYYYMMDD-HHMMSS-PID-HASH8`; `status=declared`; contract hash matches.

---

## Scenario 3: Record Artifacts in Ledger (US3)

**Goal**: Append rows to append-only TSV ledger with correct primary key and symlink handling.

```bash
# Create a minimal run directory structure
mkdir -p /tmp/prism-run/input /tmp/prism-run/alignment /tmp/prism-run/refinement

# Simulate stage outputs
echo "ATOM dummy" > /tmp/prism-run/input/5zngA_1buhAB_A_query.pdb
echo "ATOM dummy" > /tmp/prism-run/input/5zngA_1buhAB_A_template.pdb
echo "ALIGNMENT DUMMY" > /tmp/prism-run/alignment/5zngA_1buhAB_A.aln
ln -s ../alignment/5zngA_1buhAB_A.aln /tmp/prism-run/alignment/latest.aln  # symlink

# Record artifacts
python -m src.artifact_ledger write \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --run-root /tmp/prism-run \
  --stage input \
  --row-id 5zngA_1buhAB_A \
  --role query_structure \
  --rel-path input/5zngA_1buhAB_A_query.pdb

python -m src.artifact_ledger write \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --run-root /tmp/prism-run \
  --stage input \
  --row-id 5zngA_1buhAB_A \
  --role template_structure \
  --rel-path input/5zngA_1buhAB_A_template.pdb

python -m src.artifact_ledger write \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --run-root /tmp/prism-run \
  --stage alignment \
  --row-id 5zngA_1buhAB_A \
  --role alignment_output \
  --rel-path alignment/5zngA_1buhAB_A.aln

# Symlink entry
python -m src.artifact_ledger write \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --run-root /tmp/prism-run \
  --stage alignment \
  --row-id 5zngA_1buhAB_A \
  --role alignment_output \
  --rel-path alignment/latest.aln

# View ledger
cat /tmp/prism-run/artifact_ledger.tsv

# Verify primary key uniqueness
python -c "
import pandas as pd
df = pd.read_csv('/tmp/prism-run/artifact_ledger.tsv', sep='\t')
print('Rows:', len(df))
print('Unique PK:', df[['dataset_row_id','scientific_role','run_relative_path']].drop_duplicates().shape[0])
print('Symlink row:')
print(df[df['run_relative_path']=='alignment/latest.aln'].to_string())
"
```

**Expected**:
- TSV has header matching schema
- 4 rows (3 files + 1 symlink)
- Symlink row: `path_kind=symlink`, `target_sha256` = hash of target file, `link_target=../alignment/5zngA_1buhAB_A.aln`
- Primary key `(dataset_row_id, scientific_role, run_relative_path)` is unique

---

## Scenario 4: Ledger Collision Rejection (US3 - ADR-0002 Enforcement)

**Goal**: Verify duplicate primary key is rejected (not silently overwritten).

```bash
# Attempt duplicate PK: same row_id + role + rel_path
python -m src.artifact_ledger write \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --run-root /tmp/prism-run \
  --stage alignment \
  --row-id 5zngA_1buhAB_A \
  --role alignment_output \
  --rel-path alignment/5zngA_1buhAB_A.aln  # SAME PATH as earlier

# Should exit non-zero with clear error about duplicate PK
echo "Exit code: $?"
```

**Expected**: Non-zero exit; error message mentions duplicate primary key `(5zngA_1buhAB_A, alignment_output, alignment/5zngA_1buhAB_A.aln)`.

---

## Scenario 5: Missing/Unavailable Artifacts (US3 - ADR-0002)

**Goal**: Verify explicit records for missing/unavailable artifacts with `sha256=null`.

```bash
# Record a missing artifact (stage output expected but absent)
python -m src.artifact_ledger write \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --run-root /tmp/prism-run \
  --stage refinement \
  --row-id 5zngA_1buhAB_A \
  --role refined_model \
  --rel-path refinement/5zngA_1buhAB_A_refined.pdb \
  --status missing

# Record unavailable (e.g., native PDB not found)
python -m src.artifact_ledger write \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --run-root /tmp/prism-run \
  --stage comparison \
  --row-id 5zngA_1buhAB_A \
  --role native_structure \
  --rel-path comparison/5zngA_1buhAB_A_native.pdb \
  --status unavailable

# Verify null hashes
cat /tmp/prism-run/artifact_ledger.tsv | grep -E 'missing|unavailable'
python -c "
import pandas as pd
df = pd.read_csv('/tmp/prism-run/artifact_ledger.tsv', sep='\t')
missing = df[df['status'].isin(['missing','unavailable'])]
print('Missing/Unavailable rows:', len(missing))
print('All sha256 null:', missing['sha256'].isna().all())
"
```

**Expected**: Rows present with `sha256=NA` (or empty), `status=missing`/`unavailable`.

---

## Scenario 6: Validation Gate Pre-Consumption (US4)

**Goal**: Run validation gate on completed run; detect mismatches/missing.

```bash
# Run validation gate
python -m src.validation_gate validate \
  --run-root /tmp/prism-run \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --ledger /tmp/prism-run/artifact_ledger.tsv \
  --out /tmp/validation_gate.json

# Inspect
cat /tmp/validation_gate.json | python -m json.tool

# Verify structure
python -c "
import json
with open('/tmp/validation_gate.json') as f:
    v = json.load(f)
print('Overall:', v['overall_status'])
print('Artifact checks:', len(v['artifact_checks']))
print('Row identity checks:', len(v['row_identity_checks']))
print('Secret scan clean:', v['secret_scan']['clean'])
"
```

**Expected**: `overall_status` in `pass|fail|warn`; `secret_scan.clean=true`; artifact checks cover all ledger rows.

---

## Scenario 7: Artifact Tampering Detection (US4 - Mutation Detection)

**Goal**: Verify gate detects content hash mismatch after file modification.

```bash
# Modify a file after ledger recorded
echo "TAMPERED" > /tmp/prism-run/input/5zngA_1buhAB_A_query.pdb

# Re-run validation
python -m src.validation_gate validate \
  --run-root /tmp/prism-run \
  --run-id prism-20260729-120000-12345-abcd1234 \
  --ledger /tmp/prism-run/artifact_ledger.tsv \
  --out /tmp/validation_gate2.json

# Check detection
python -c "
import json
with open('/tmp/validation_gate2.json') as f:
    v = json.load(f)
mismatches = [c for c in v['artifact_checks'] if c['status']=='mismatch']
print('Mismatches detected:', len(mismatches))
for m in mismatches:
    print(f\"  {m['dataset_row_id']} | {m['scientific_role']} | expected={m['expected_sha256']} actual={m['actual_sha256']}\")
"
```

**Expected**: At least 1 mismatch detected; `expected_sha256` ≠ `actual_sha256`.

---

## Scenario 8: Retry with Superseded Attempt (US2)

**Goal**: New run_id links to failed parent via `parent_run_id`.

```bash
# First attempt (simulate failure)
python -m src.run_identity init-run \
  --declared-contract-hash $DECLARED_HASH \
  --declared-contract-path /tmp/declared-contract.json \
  --out /tmp/run-attempt1.json
ATTEMPT1_ID=$(python -c "import json; print(json.load(open('/tmp/run-attempt1.json'))['run_id'])")

# Mark as failed (manual for test)
python -c "
import json
with open('/tmp/run-attempt1.json') as f: m = json.load(f)
m['status'] = 'failed'
with open('/tmp/run-attempt1.json', 'w') as f: json.dump(m, f, indent=2)
"

# Retry: new run_id, links to parent
python -m src.run_identity init-run \
  --declared-contract-hash $DECLARED_HASH \
  --declared-contract-path /tmp/declared-contract.json \
  --parent-run-id $ATTEMPT1_ID \
  --supersedes-reason 'stage refinement failed; retry with relaxed gate' \
  --out /tmp/run-attempt2.json

python -c "
import json
with open('/tmp/run-attempt2.json') as f: m = json.load(f)
print('Attempt 2 run_id:', m['run_id'])
print('Parent:', m['parent_run_id'])
print('Reason:', m['supersedes_reason'])
print('Status:', m['status'])
"
```

**Expected**: New `run_id`; `parent_run_id` = attempt1; `supersedes_reason` set; `status=declared`.

---

## Scenario 9: Secret Redaction in Manifest (US5)

**Goal**: Verify Git provenance captures diff + untracked; env vars redacted except allowlist.

```bash
# Create a test repo state
cd /tmp
mkdir test_prism_repo && cd test_prism_repo
git init -q
echo "PRISM_API_KEY=secret123" > .env
echo "PRISM_MULTIPROT_FORCE=1" >> .env
git add .env && git commit -m "init" -q
echo "untracked_config.txt" > untracked_config.txt
echo "PRISM_API_KEY=secret123" >> .env
git add -N .  # stage for diff

# Run provenance capture (from investigation_provenance.py pattern)
cd /scratch/rshadi25/GitHub/PRISM-prescript
python -c "
from benchmark.scripts.investigation_provenance import capture_git_provenance, redact_secrets
import os
os.chdir('/tmp/test_prism_repo')
prov = capture_git_provenance()
print('Git provenance:')
for k, v in prov.items():
    print(f'  {k}: {v}')

# Test secret redaction
test_env = {'PRISM_API_KEY': 'secret', 'PATH': '/usr/bin', 'SLURM_JOB_ID': '123', 'MY_SECRET': 'leak'}
redacted = redact_secrets(test_env)
print('Redacted env:', redacted)
"

# Expected: .env diff captured, untracked declared; PRISM_API_KEY → [REDACTED], PATH/SLURM_JOB_ID kept
```

**Expected**: `git_diff_hash` non-empty; `declared_untracked_files` includes `untracked_config.txt`; allowlisted vars preserved; secret vars redacted.

---

## Scenario 10: Canonical JSON Determinism (US1)

**Goal**: Verify identical DeclaredContract → identical hash regardless of input key order.

```bash
python -c "
import json, hashlib

d1 = {'b': 2, 'a': 1}
d2 = {'a': 1, 'b': 2}

def canon(d):
    return json.dumps(d, sort_keys=True, separators=(',', ':')).encode('utf-8')

h1 = hashlib.sha256(canon(d1)).hexdigest()
h2 = hashlib.sha256(canon(d2)).hexdigest()
print('Hash 1:', h1)
print('Hash 2:', h2)
print('Equal:', h1 == h2)
"

# Also test nan handling
python -c "
import json, hashlib, math
d = {'x': float('nan')}
try:
    canon = json.dumps(d, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')
    print('Should have failed')
except ValueError as e:
    print('Correctly rejects NaN:', e)
"
```

**Expected**: Hashes equal; NaN rejected with `ValueError`.

---

## Validation Checklist

After all scenarios, confirm:

| Check | Pass? |
|-------|-------|
| DeclaredContract canonical hash stable | ☐ |
| RunIdentity run_id pattern correct | ☐ |
| Ledger TSV primary key enforced | ☐ |
| Symlink rows hash target content | ☐ |
| Missing/unavailable = null hash | ☐ |
| Validation gate detects tampering | ☐ |
| Validation gate rejects secret leaks | ☐ |
| Retry creates new run_id with parent link | ☐ |
| All JSON validates against schemas | ☐ |

---

## Notes for Implementers

- Scenarios 1-2 require `src/run_identity.py` with `declare-contract` and `init-run` subcommands
- Scenarios 3-5 require `src/artifact_ledger.py` with `write` subcommand (append-only)
- Scenarios 6-7 require `src/validation_gate.py` with `validate` subcommand
- Scenario 8 uses `src/run_identity.py` with `--parent-run-id` and `--supersedes-reason`
- Scenario 9 uses `benchmark/scripts/investigation_provenance.py` patterns (extend)
- All JSON outputs must validate against contracts in `contracts/`
- TSV ledger must be append-only; never rewrite history
- Run `python -m pytest tests/contract/ -v` after implementation
