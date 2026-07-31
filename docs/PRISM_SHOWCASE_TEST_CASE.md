# PRISM Showcase Test Case: Two Orientations, One Ranked Candidate

## What this case demonstrates

This retained test case shows one receptor–ligand query with two transformed
candidate orientations. PRODIGY scores both candidates successfully, and the
PRISM ranking adapter forwards the single candidate with the more favorable
predicted affinity when `top_k=1`.

It demonstrates:

- candidate identity encoded in transformed filenames;
- construction of a combined complex for an external scorer;
- command, input hash, stdout/stderr, return-code, and score retention;
- top-K selection before an expensive refinement stage;
- the difference between mechanical ranking behavior and biological quality.

It does **not** demonstrate that PRODIGY improves DockQ, that orientation `o1`
is native-like, or that ranking should be enabled by default. No native DockQ
label is part of this fixture.

## Case contract

| Field | Value |
| --- | --- |
| Query receptor | `5zngA` |
| Query ligand | `4eylA` |
| Template/interface | `1a0cCD` |
| Candidates before ranking | `o1`, `o2` |
| Ranking method | PRODIGY predicted affinity |
| Selection rule | Lower/more negative affinity first |
| `top_k` | `1` |
| Expected selected orientation | `o1` |
| Evidence root | `tmp/agent/20260730-prodigy-ranking-paired-test/` |
| Authoritative summary | `summary-corrected.json` |

Expected evidence:

| Orientation | Status | Return code | Affinity (kcal/mol) | Combined-input SHA256 | Forwarded? |
| --- | --- | ---: | ---: | --- | --- |
| `o1` | `scored` | 0 | -65.827 | `8cbf99f996fb65094cfd065f053bb51ad6157de52ca3831285f99f4c55de5274` | Yes |
| `o2` | `scored` | 0 | -65.274 | `a222886f4be249c263240cf684ed91a1d41f5f4d6cd088ad2452613b3c3f4fd0` | No |

The score difference is `0.553 kcal/mol`. The adapter chooses `o1` because
`-65.827 < -65.274`; it does not choose by orientation name or file order.

## Presentation route: 8–10 minutes

### Step 1: Frame the hypothesis

Tell the audience:

> We have two geometrically accepted orientations. If the PRODIGY adapter is
> working and `top_k=1`, both candidates must be scored, but only the candidate
> with the more negative affinity should reach refinement.

Ask before showing the result:

> Which fields would you require to distinguish a legitimate top-1 selection
> from silently dropping a failed candidate?

Look for: status for both candidates, commands, input hashes, return codes,
scores, group membership, and the complete pre-ranking denominator.

### Step 2: Verify the retained fixture

From the repository root, define the evidence path:

```bash
case_root=tmp/agent/20260730-prodigy-ranking-paired-test
test -f "$case_root/summary-corrected.json"
find "$case_root/scores-corrected" -maxdepth 1 -type f -printf '%f\n' | sort
```

This is read-only and safe for a live presentation. Expect two sets of:

- combined `.pdb` input;
- `.json` structured score record;
- `.stdout.txt` raw scorer output;
- `.stderr.txt` raw scorer diagnostics.

Do not use `summary.json` as the result of record. It preserves the earlier
pre-fix behavior in which the PRODIGY input PDB was placed after `--selection`
and both candidates were retained. `summary-corrected.json` records the fixed
argument order and the correct top-1 behavior.

### Step 3: Show the two transformed candidate identities

```bash
ls -lh \
  processed/transformation/1a0cCD_5zngA_4eylA_o1_L.pdb \
  processed/transformation/1a0cCD_5zngA_4eylA_o1_R.pdb \
  processed/transformation/1a0cCD_5zngA_4eylA_o2_L.pdb \
  processed/transformation/1a0cCD_5zngA_4eylA_o2_R.pdb
```

Decode one filename with the audience:

```text
1a0cCD_5zngA_4eylA_o1_L.pdb
|      |     |     |  |
|      |     |     |  +-- left transformed partner
|      |     |     +----- orientation 1
|      |     +----------- query ligand
|      +----------------- query receptor
+------------------------ template/interface
```

Ask: **Why must `o1_L` and `o1_R` be kept together?**

Expected answer: they are the paired partners of one candidate orientation;
mixing partners across orientations would create a different, undeclared
complex.

### Step 4: Inspect one structured score record

```bash
python -m json.tool \
  "$case_root/scores-corrected/1a0cCD_5zngA_4eylA_o1_L__1a0cCD_5zngA_4eylA_o1_R.json"
```

Point out:

- `status: scored`;
- `return_code: 0`;
- `affinity_kcal_mol: -65.827`;
- combined-input SHA256;
- original left and right transformed paths;
- exact PRODIGY argument order;
- raw stdout and stderr paths.

The input PDB appears before `--selection`. This ordering matters because
PRODIGY's selection option consumes following tokens as chain selections.

### Step 5: Assert the expected result automatically

The following read-only check fails loudly if the retained showcase contract
drifts:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python - <<'PY'
import json
from pathlib import Path

summary_path = Path(
    "tmp/agent/20260730-prodigy-ranking-paired-test/summary-corrected.json"
)
data = json.loads(summary_path.read_text())

assert data["case"] == "5zngA,4eylA"
assert data["template"] == "1a0cCD"
assert data["without_ranking_forwarded_count"] == 2
assert data["with_prodigy_top_k"] == 1
assert data["with_prodigy_forwarded_count"] == 1
assert data["changed_selection_count"] is True

scores = {row["left_pdb"].split("_")[-2]: row for row in data["scores"]}
assert set(scores) == {"o1", "o2"}
assert all(row["status"] == "scored" for row in scores.values())
assert all(row["return_code"] == 0 for row in scores.values())
assert scores["o1"]["affinity_kcal_mol"] == -65.827
assert scores["o2"]["affinity_kcal_mol"] == -65.274
assert scores["o1"]["affinity_kcal_mol"] < scores["o2"]["affinity_kcal_mol"]
assert "_o1_" in data["with_prodigy_forwarded"][0][0]

for row in data["scores"]:
    for key in ("input_pdb", "stdout_path", "stderr_path"):
        assert Path(row[key]).is_file(), (key, row[key])

print("PASS: both candidates scored; top-1 correctly forwards orientation o1")
PY
```

Expected terminal line:

```text
PASS: both candidates scored; top-1 correctly forwards orientation o1
```

### Step 6: Connect the result to source code

```bash
rg -n "def (select_top_candidates|select_top_candidates_with_prodigy|score_candidate|combine_pair_pdbs|parse_affinity)" \
  src/candidate_selector.py src/prodigy_ranker.py
```

Trace in this order:

1. `src/candidate_selector.py:select_top_candidates()` dispatches on
   `rank_method="prodigy"`.
2. `src/prodigy_ranker.py:select_top_candidates_with_prodigy()` groups and
   scores candidates.
3. `combine_pair_pdbs()` constructs one chain-renamed scorer input.
4. `score_candidate()` runs PRODIGY and retains execution evidence.
5. `parse_affinity()` extracts the numeric prediction.
6. The group is sorted by increasing affinity and truncated to `top_k` only
   after all candidates score successfully.

### Step 7: State the conclusion precisely

Observation:

- Both candidates received successful scores.
- `o1` had the more negative predicted affinity.
- Top-1 selection forwarded `o1` and reduced the next-stage candidate count
  from two to one.

Inference:

- The PRODIGY adapter and top-K forwarding logic worked mechanically for this
  retained pair.

Not established:

- `o1` is closer to the native complex;
- fewer candidates reduced wall time in a controlled refinement run;
- PRODIGY generalizes across independent complexes;
- PRODIGY should replace the stable unranked reference path.

## Optional code-level regression check

This test uses a fake scorer for deterministic unit coverage; it does not run
the external PRODIGY executable:

```bash
/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q \
  tests/test_prodigy_ranker.py
```

Expected behaviors covered include combined-chain selection, correct CLI
argument ordering, lower-affinity ranking, and preservation of a group when a
candidate score fails.

## Full pipeline extension: Slurm only

If the session needs to show how this ranking mode is requested from
`prism.py`, present the command without launching it interactively:

```bash
PRISM_INPUTS_CSV=<isolated-case-inputs.csv> \
PRISM_PRODIGY_EXECUTABLE=<verified-prodigy-executable> \
PRISM_STAGE_STATUS_PATH=<run-root>/stage-status.jsonl \
/home/rshadi25/.conda/envs/gtalign_env/bin/python prism.py \
  --template-limit 1 \
  --rank true \
  --rank-method prodigy \
  --top-k 1 \
  --prodigy-executable <verified-prodigy-executable>
```

This is a template, not a ready-to-run reproduction contract: the first
calculated template is not guaranteed to be `1a0cCD`, and a complete isolated
workspace must stage the intended template assets, source identity, input
structures, runtime, and refiner. Prepare that workspace explicitly and use
Slurm for alignment/refinement. Do not mutate the repository's shared
`inputs.csv`, templates, or validated outputs for a presentation.

## Failure variants to discuss

| Failure | Expected safe behavior | What not to do |
| --- | --- | --- |
| One candidate has nonzero return code | Preserve the full group rather than rank incomplete evidence. | Do not silently promote the only successful candidate. |
| Combined-input hash changes | Treat it as a different artifact and investigate provenance. | Do not reuse the old score as if input bytes were unchanged. |
| Score cannot be parsed | Record failure, stdout, stderr, and command. | Do not substitute `0` or an invented affinity. |
| Only one transformed partner exists | Candidate is incomplete and not scoreable. | Do not combine it with the other orientation's partner. |
| `o2` scores better in a rerun | Verify tool/version/input hashes before interpreting stochastic or environmental drift. | Do not overwrite retained evidence. |

## Cleanup policy

The showcase is read-only. Do not delete or rewrite
`tmp/agent/20260730-prodigy-ranking-paired-test/`. If a live replay is approved,
write to a new directory under `tmp/agent/<timestamp>-prodigy-showcase-replay/`
and retain its command, environment, source hash, inputs, outputs, and status.

## Learning recap

- **Concept:** ranking is a selection operation over a complete candidate
  group, not a substitute for native-complex evaluation.
- **Pitfall:** a successful score is meaningless without exact candidate and
  input identity.
- **Exercise:** change the hypothetical `top_k` to 2 and predict what should be
  forwarded without running anything. Then explain why the scientific quality
  claim remains unchanged.
