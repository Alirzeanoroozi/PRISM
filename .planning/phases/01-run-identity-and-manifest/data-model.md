# Data Model: Run Identity and Manifest

**Phase**: 01-run-identity-and-manifest  
**Generated**: 2026-07-29

---

## Entities

### 1. DeclaredContract

**Purpose**: Immutable scientific declaration. Hash of canonical JSON = `declared_contract_hash`. No runtime/attempt fields.

| Field | Type | Constraints | Description |
|-------|------|-------------|-------------|
| `contract_version` | string | const: "1.0" | Schema version |
| `pipeline_version` | string | SemVer from `prism.py` | Pipeline code version |
| `stages_enabled` | array[string] | Subset: `input`, `surface`, `alignment`, `transformation`, `refinement`, `comparison` | Declared active stages |
| `aligner` | string | Enum: `tmalign`, `gtalign`, `multiprot` | Primary aligner |
| `refiner` | string | Enum: `external_rosetta`, `pyrosetta`, `fiberdock`, `none` | Primary refiner |
| `input_selectors` | InputSelectors | Required | Raw + normalized dataset selectors |
| `template_inventory` | array[TemplateAsset] | Required | Template assets with hashes |
| `parameters` | object | JSON-serializable | All user-facing config knobs |
| `resource_request` | ResourceRequest | Required | CPU, memory, time, GPU, partition |
| `source_inventory` | SourceInventory | Required | Git HEAD + diff + declared untracked |
| `tool_fingerprints` | array[ToolFingerprint] | Required | External tool versions + hashes |

**Canonicalization**: JSON with `sort_keys=True`, `separators=(',', ':')`, `allow_nan=False`, UTF-8 encoded before SHA256.

---

### 2. RunIdentity (Execution Attempt)

**Purpose**: Immutable identity for a single operational pipeline attempt.

| Field | Type | Constraints | Description |
|-------|------|-------------|-------------|
| `run_id` | string | Pattern: `^prism-\d{8}-\d{6}-\d+-[a-f0-9]{8}$` | Human-readable unique ID: `prism-YYYYMMDD-HHMMSS-PID-HASH8` |
| `declared_contract_hash` | string | SHA256 (64 hex chars) | Links to DeclaredContract |
| `created_at` | datetime | ISO 8601 UTC | When run was initialized |
| `host` | string | | Hostname where run initiated |
| `user` | string | | Unix user who launched the run |
| `status` | enum | `declared`, `running`, `completed`, `failed`, `partial`, `validated` | Attempt status |
| `parent_run_id` | string? | Same pattern as run_id | For retries: links to superseded attempt |
| `supersedes_reason` | string? | Required if parent_run_id present | Why this attempt supersedes parent |

**Validation Rules**:
- `run_id` generated once at pipeline entry (`prism.py` init)
- `declared_contract_hash` computed from DeclaredContract canonical JSON
- `parent_run_id` + `supersedes_reason` required together; never mutate parent ledger

---

### 3. InputSelectors

| Field | Type | Description |
|-------|------|-------------|
| `raw` | array[string] | User-provided pair selectors exactly as declared |
| `normalized` | array[NormalizedSelector] | Parsed, validated selector tuples consumed by execution |

### 4. NormalizedSelector

| Field | Type | Description |
|-------|------|-------------|
| `row_id` | string | Durable benchmark row identity (e.g., `5zngA_1buhAB_A`) |
| `query_chain` | string | Query PDB code + chain(s) |
| `template_chain` | string | Template PDB code + chain(s) |

---

### 5. TemplateAsset

| Field | Type | Description |
|-------|------|-------------|
| `asset_id` | string | Logical identifier (e.g., `1a0cCD_prepared.pdb`) |
| `role` | enum | `template_pdb`, `template_interface`, `native_pdb`, `surface_file`, `other` |
| `logical_path` | string | Path as declared in run |
| `resolved_path` | string | Absolute resolved path |
| `sha256` | string | SHA256 of file content (64 hex chars) |
| `size_bytes` | integer | File size ≥ 0 |

---

### 6. ResourceRequest

| Field | Type | Description |
|-------|------|-------------|
| `cpus` | integer | CPU cores requested |
| `memory_gb` | integer | Memory in GB |
| `time_hours` | number | Wall time limit in hours |
| `partition` | string | Slurm partition (e.g., `ai`) |
| `qos` | string? | Slurm QoS (e.g., `ai`) |
| `gpu` | boolean | GPU requested |

---

### 7. SourceInventory

| Field | Type | Description |
|-------|------|-------------|
| `git_head` | string | Full commit SHA (40 hex chars) |
| `git_diff_hash` | string | SHA256 of `git diff HEAD` (empty if clean) |
| `declared_untracked` | array[string] | Declared untracked files that affect run |

---

### 8. ToolFingerprint

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Tool identifier (e.g., `tmalign`, `rosetta_scripts`) |
| `version` | string | Version string from probe |
| `sha256` | string? | File hash if readable |

---

### 9. RuntimeEnvironment (in Execution Attempt)

| Field | Type | Description |
|-------|------|-------------|
| `command_argv` | array[string] | Structured argv (not shell string) |
| `working_directory` | string | Absolute path of run root |
| `environment` | object | Allowlisted env vars only (PRISM_*, PATH, HOME, SLURM_*, CONDA_*, CUDA_VISIBLE_DEVICES, OMP_NUM_THREADS, MKL_NUM_THREADS) |
| `python_version` | string | `sys.version` |
| `python_executable` | string | `sys.executable` |
| `conda_env` | string? | Conda env name if active |
| `seeds` | object | Explicit seeds: `PYTHONHASHSEED`, `numpy`, `random`, etc. |
| `packages` | array[PackageVersion] | Key packages with versions |
| `external_tools` | array[ExecutableProbe] | Resolved external tools |
| `slurm` | SlurmContext? | Null if local run |

---

### 10. PackageVersion

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Package name |
| `version` | string | `__version__` or importlib.metadata |
| `location` | string? | Install path |

---

### 11. ExecutableProbe

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Tool identifier (e.g., `tmalign`, `gtalign`, `rosetta`) |
| `resolved_path` | string | Absolute path from `which`/`shutil.which` |
| `version_output` | string | Stdout of `--version` or probe command |
| `sha256` | string? | File hash if readable |
| `probe_ok` | boolean | Whether version probe succeeded |

---

### 12. SlurmContext

| Field | Type | Description |
|-------|------|-------------|
| `job_id` | string | `SLURM_JOB_ID` |
| `job_name` | string | `SLURM_JOB_NAME` |
| `partition` | string | `SLURM_JOB_PARTITION` |
| `qos` | string | `SLURM_JOB_QOS` |
| `account` | string | `SLURM_JOB_ACCOUNT` |
| `nodes` | array[string] | `SLURM_JOB_NODELIST` expanded |
| `cpus_per_task` | integer | `SLURM_CPUS_PER_TASK` |
| `mem_per_cpu` | string | `SLURM_MEM_PER_CPU` |
| `time_limit` | string | `SLURM_TIMELIMIT` |
| `submit_dir` | string | `SLURM_SUBMIT_DIR` |

---

### 13. GitProvenance (in RuntimeEnvironment)

| Field | Type | Description |
|-------|------|-------------|
| `head_commit` | string | Full SHA |
| `head_short` | string | First 8 chars |
| `branch` | string? | Current branch name |
| `is_dirty` | boolean | True if worktree has changes |
| `status_short` | string | `git status --short` output |
| `diff_hash` | string? | SHA256 of `git diff HEAD` (empty if clean) |
| `declared_untracked_files` | array[string] | Declared untracked files that affect run |

---

## Artifact Ledger Entities

### 14. ArtifactRecord

**Purpose**: Single row in append-only TSV ledger. One per materialized scientific artifact.

| Field | Type | Constraints | Description |
|-------|------|-------------|-------------|
| `run_id` | string | FK → RunIdentity.run_id | Parent run |
| `stage` | string | Enum: `input`, `surface`, `alignment`, `transformation`, `refinement`, `comparison`, `scoring`, `other` | Pipeline stage |
| `dataset_row_id` | string | | Durable benchmark row ID (e.g., `5zngA_1buhAB_A`) |
| `scientific_role` | string | Enum: `query_structure`, `template_structure`, `template_interface`, `native_structure`, `surface_file`, `alignment_output`, `transformation_output`, `refined_model`, `scored_model`, `evaluation_output`, `other` | Role in pipeline |
| `run_relative_path` | string | Relative to run root | Logical path |
| `path_kind` | string | Enum: `file`, `symlink`, `directory`, `missing`, `unavailable`, `broken_link` | FS entry type |
| `size_bytes` | integer | ≥0 | Target content size |
| `sha256` | string? | 64 hex chars or null | Content hash of consumed bytes (symlink = target hash; missing/unavailable = null) |
| `target_sha256` | string? | 64 hex chars or null | Hash of symlink target file; null for non-symlinks |
| `link_target` | string? | | If symlink: raw link target text |
| `resolved_path` | string? | | If symlink: absolute resolved path |
| `is_broken_link` | boolean | default: false | True if symlink target missing |
| `materialized_at` | datetime | ISO 8601 UTC | When observed |
| `produced_by` | string | | Stage/module that produced this (e.g., `alignment`, `tmalign`, `pyrosetta_refinement`) |
| `status` | string | Enum: `ok`, `missing`, `unavailable`, `failed`, `partial`, `not_scoreable`, `error` | Observation status |
| `error_detail` | string? | Required if status ≠ `ok` | Failure reason |

**Uniqueness (Primary Key)**: `(dataset_row_id, scientific_role, run_relative_path)` — **duplicate = error**

---

### 15. ValidationResult

**Purpose**: Output of pre-consumption validation gate.

| Field | Type | Description |
|-------|------|-------------|
| `run_id` | string | Validated run |
| `validated_at` | datetime | When validation ran |
| `overall_status` | enum | `pass`, `fail`, `warn` |
| `artifact_checks` | array[ArtifactCheck] | Per-artifact hash/integrity |
| `row_identity_checks` | array[RowIdentityCheck] | Duplicate/mismatch detection |
| `secret_scan` | SecretScanResult | Provenance secret leakage check |

---

### 16. ArtifactCheck

| Field | Type | Description |
|-------|------|-------------|
| `dataset_row_id` | string | |
| `scientific_role` | string | |
| `run_relative_path` | string | |
| `expected_sha256` | string? | From ledger (null if missing) |
| `actual_sha256` | string? | Current file hash (null if missing) |
| `match` | boolean | True if hashes equal |
| `status` | enum | `ok`, `mismatch`, `missing`, `unavailable` |

---

### 17. RowIdentityCheck

| Field | Type | Description |
|-------|------|-------------|
| `dataset_row_id` | string | Benchmark row |
| `expected_roles` | array[string] | Roles declared in ledger for this row |
| `found_roles` | array[string] | Roles with status=`ok` |
| `missing_roles` | array[string] | Roles with status=`missing`/`unavailable` |
| `duplicate_keys` | array[string] | Any duplicate primary keys found |
| `status` | enum | `complete`, `incomplete`, `duplicate`, `mismatch` |

---

### 18. SecretScanResult

| Field | Type | Description |
|-------|------|-------------|
| `clean` | boolean | No sentinels found |
| `scanned_fields` | array[string] | Fields checked |
| `sentinel` | string | Test value used |

---

## State Transitions

### RunIdentity (Execution Attempt)
```
CREATED → (manifest declared) → DECLARED_CONTRACT_HASHED → RUNNING → (stages complete) → CLOSED
                                                    ↓
                                              (retry/correction)
                                                    ↓
                                               SUPERSEDED (new RunIdentity with parent_run_id)
```

### ArtifactRecord
```
PENDING (stage not yet reached) → MATERIALIZED (ok) → VALIDATED (gate pass)
         ↓                              ↓
    SKIPPED (stage disabled)        MISMATCH (gate fail → block consumption)
```

---

## Relationships

```
DeclaredContract (1) ───< (many) RunIdentity (execution attempts)
RunIdentity (1) ───< (many) ArtifactRecord
RunIdentity (1) ───< (1) ValidationResult
ValidationResult (1) ───< (many) ArtifactCheck
ValidationResult (1) ───< (many) RowIdentityCheck
ArtifactRecord.record_id = concat(dataset_row_id, "|", scientific_role, "|", run_relative_path)
```

---

## Serialization

| Entity | Primary Format | Secondary |
|--------|----------------|-----------|
| DeclaredContract | JSON (`declared-contract.json`) | Hashed for declared_contract_hash |
| RunIdentity | JSON (`run-manifest.json`) | — |
| ArtifactRecord | TSV (`artifact_ledger.tsv`) | JSON Lines for programmatic |
| ValidationResult | JSON (`validation_gate.json`) | — |
| GitProvenance | Embedded in RunIdentity | — |

All timestamps: ISO 8601 UTC (`YYYY-MM-DDTHH:MM:SSZ`).
All hashes: lowercase hex SHA256.
