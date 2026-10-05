#!/usr/bin/env python3
"""Write the frozen paper/dataset/reference crosswalk for an investigation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = (
    "claim_id",
    "claim",
    "arm_or_cohort",
    "observed_value",
    "expected_or_reference_value",
    "reference_file",
    "reference_anchor",
    "status",
    "evidence_required",
)


ROWS = (
    {
        "claim_id": "paper-bm3-cohort",
        "claim": "PRISM paper validation cohort",
        "arm_or_cohort": "paper",
        "observed_value": "88 rigid-body cases from 165 chains",
        "expected_or_reference_value": "88 Benchmark 3.0 interfaces",
        "reference_file": "references/Proteins - 2011 - Tuncbag - Fast and accurate modeling of protein protein interactions by combining.md",
        "reference_anchor": "line 260; lines 641-642",
        "status": "blocked-until-source-freeze",
        "evidence_required": "authoritative 88-row list with raw selectors and hashes",
    },
    {
        "claim_id": "paper-unbiased-templates",
        "claim": "Paper unbiased template set",
        "arm_or_cohort": "paper",
        "observed_value": "not reconstructed from repository BM5/5.5 table",
        "expected_or_reference_value": "7,922 protein interfaces with self-hit removal",
        "reference_file": "references/Proteins - 2011 - Tuncbag - Fast and accurate modeling of protein protein interactions by combining.md",
        "reference_anchor": "lines 270-275; lines 350-375",
        "status": "blocked",
        "evidence_required": "authoritative template list and exclusion policy",
    },
    {
        "claim_id": "archive-aliases",
        "claim": "Synthetic archive prefixes identify repeated complexes",
        "arm_or_cohort": "repository BM5/5.5 extension",
        "observed_value": "BAAD, BOYV, BP57, CP57 resolved by chain-set matching",
        "expected_or_reference_value": "README-defined synthetic prefixes",
        "reference_file": "benchmark/originals/benchmark5.5/README",
        "reference_anchor": "FILE NAMING FORMAT; synthetic-prefix paragraph",
        "status": "supported",
        "evidence_required": "four role-file hashes per dataset_row_id",
    },
    {
        "claim_id": "naccess-surface",
        "claim": "Surface accessibility tool and threshold",
        "arm_or_cohort": "historical/current implementation",
        "observed_value": "NACCESS; RSA threshold and scaffold are implementation parameters",
        "expected_or_reference_value": "NACCESS surface extraction; protocol threshold crosswalk required",
        "reference_file": "references/nprot.2011.367.md",
        "reference_anchor": "NACCESS/surface-extraction methods section",
        "status": "crosswalk-required",
        "evidence_required": "effective configuration and executable hash",
    },
    {
        "claim_id": "fiberdock-output",
        "claim": "Flexible refinement and ranking output",
        "arm_or_cohort": "working_version/multiprot",
        "observed_value": "FiberDock output and within-arm energy ranking",
        "expected_or_reference_value": "FiberDock structures/energies; do not compare energy magnitudes across tools",
        "reference_file": "references/Proteins - 2011 - Tuncbag - Fast and accurate modeling of protein protein interactions by combining.md",
        "reference_anchor": "lines 198-205; lines 350-364",
        "status": "supported-with-boundary",
        "evidence_required": "adapter-derived coordinate and energy hashes",
    },
    {
        "claim_id": "modern-dockq",
        "claim": "Modern common evaluation metric",
        "arm_or_cohort": "all confirmatory arms",
        "observed_value": "DockQ 2.1.3 evaluator contract",
        "expected_or_reference_value": "implementation-specific modern secondary metric",
        "reference_file": "references/Proteins - 2011 - Gao - New benchmark metrics for protein‐protein docking methods - IScore.md",
        "reference_anchor": "paper metrics crosswalk; DockQ is not claimed as paper reproduction",
        "status": "implementation-specific",
        "evidence_required": "evaluator version, mapping fixture results, raw JSON hashes",
    },
    {
        "claim_id": "repository-cohort-size",
        "claim": "Primary repository benchmark cohort",
        "arm_or_cohort": "repository BM5/5.5 extension",
        "observed_value": "257 rows: 162 rigid, 60 medium, 35 difficult",
        "expected_or_reference_value": "the three benchmark CSVs currently present in benchmark/data",
        "reference_file": "benchmark/data/T_Rigid.csv; benchmark/data/T_medium.csv; benchmark/data/T_difficult.csv",
        "reference_anchor": "header and complete row counts",
        "status": "source-locked",
        "evidence_required": "dataset_row_id manifest with raw CSV row numbers and hashes",
    },
    {
        "claim_id": "curated-role-contract",
        "claim": "Pipeline inputs and native truth use curated archive roles",
        "arm_or_cohort": "repository BM5/5.5 extension",
        "observed_value": "pipeline_receptor=r_u; pipeline_ligand=l_u; native_receptor=r_b; native_ligand=l_b",
        "expected_or_reference_value": "benchmark file naming format defines r/l and u/b role suffixes",
        "reference_file": "benchmark/originals/benchmark5.5/README",
        "reference_anchor": "FILE NAMING FORMAT; characters 6-8",
        "status": "source-locked",
        "evidence_required": "archive member names plus four role hashes per dataset_row_id",
    },
    {
        "claim_id": "paper-matching-threshold",
        "claim": "Paper template matching and candidate filters",
        "arm_or_cohort": "paper method crosswalk",
        "observed_value": "implementation must record threshold values per arm",
        "expected_or_reference_value": "paper reports 40% residue matching, 60% for short template chains, minimum 15, hotspot, clash and accessibility filters",
        "reference_file": "references/Proteins - 2011 - Tuncbag - Fast and accurate modeling of protein protein interactions by combining.md",
        "reference_anchor": "methods: template matching, NACCESS, collision and hotspot filters",
        "status": "crosswalk-required",
        "evidence_required": "effective arm configuration and filter decision records",
    },
    {
        "claim_id": "protocol-matching-threshold",
        "claim": "Historical protocol threshold differs from paper implementation",
        "arm_or_cohort": "working_version/multiprot",
        "observed_value": "historical protocol reference uses 50%/30% matching and DIFF_PERCENTAGE 20",
        "expected_or_reference_value": "protocol-specific values must not be silently applied to the paper arm",
        "reference_file": "references/nprot.2011.367.md",
        "reference_anchor": "template matching and DIFF_PERCENTAGE sections",
        "status": "difference-to-test",
        "evidence_required": "historical/current filter replay on identical candidate streams",
    },
    {
        "claim_id": "template-panel",
        "claim": "Matched template panel for causal comparison",
        "arm_or_cohort": "current/historical diagnostic arms",
        "observed_value": "946 current templates are fully resolvable and are a subset of the historical inventory",
        "expected_or_reference_value": "candidate generation and filter replay use the same 946-template panel",
        "reference_file": "tmp/agent/20260713-historical-current-investigation/final-preflight/template_preflight_arms.json",
        "reference_anchor": "generated preflight JSON; source hashes recorded in provenance manifest",
        "status": "preflight-supported",
        "evidence_required": "template ID list and per-asset hashes for each arm",
    },
    {
        "claim_id": "tmalign-binary-identity",
        "claim": "TM-align executable identity",
        "arm_or_cohort": "current/historical diagnostic arms",
        "observed_value": "root and working_version/TMalign must be compared by hash before attribution",
        "expected_or_reference_value": "identical executable with different filters is not a TM-align-only intervention",
        "reference_file": "benchmark/scripts/investigation_provenance.py; working_version/TMalign",
        "reference_anchor": "provenance hash manifest and executable bytes",
        "status": "implementation-contract",
        "evidence_required": "executable SHA-256 and filter/config hashes",
    },
    {
        "claim_id": "top20-ranking",
        "claim": "Equal top-20 budget across confirmatory arms",
        "arm_or_cohort": "all confirmatory arms",
        "observed_value": "top-20 selection is native-independent; oracle GlobalDockQ is evaluated only after selection",
        "expected_or_reference_value": "no native score can enter candidate ranking or budget allocation",
        "reference_file": "benchmark/scripts/investigation_contracts.py",
        "reference_anchor": "native-independent ranking and top-k contract functions",
        "status": "implementation-contract",
        "evidence_required": "ranking-key records and top-1/3/5/20 selected-pose manifests",
    },
    {
        "claim_id": "no-model-null-semantics",
        "claim": "Failure and no-model score semantics",
        "arm_or_cohort": "all confirmatory arms",
        "observed_value": "structural metrics remain null; unconditional best_GlobalDockQ_at_20 is zero only for no-model utility",
        "expected_or_reference_value": "failed pairs cannot become zero-valued structural measurements",
        "reference_file": "benchmark/scripts/investigation_contracts.py; benchmark/scripts/standardized_evaluator.py",
        "reference_anchor": "score contract and normalized DockQ records",
        "status": "implementation-contract",
        "evidence_required": "adversarial no-model and malformed-mapping fixtures",
    },
    {
        "claim_id": "kutem-profile",
        "claim": "Isolated execution resource profile",
        "arm_or_cohort": "all batch tasks",
        "observed_value": "partition/account/qos kutem, one node/task, two CPUs, 2G, five minutes, array 1-10",
        "expected_or_reference_value": "same profile for smoke/calibration and at most ten tasks per array",
        "reference_file": "benchmark/jobs/isolated_kutem_array.sbatch",
        "reference_anchor": "SBATCH resource directives",
        "status": "execution-contract",
        "evidence_required": "exit.json requested/observed Slurm resources and sacct reconciliation",
    },
    {
        "claim_id": "multiprot-side-effects",
        "claim": "MultiProt compatibility boundary",
        "arm_or_cohort": "working_version/multiprot",
        "observed_value": "database, mail, HTML, destructive cleanup, and Python 2 behavior require an explicit adapter gate",
        "expected_or_reference_value": "algorithmic alignment/filtering behavior is preserved while side effects are disabled",
        "reference_file": "working_version/multiprot/mainController.py; working_version/multiprot/prism.ini",
        "reference_anchor": "controller imports/calls and local path configuration",
        "status": "blocked-until-smoke",
        "evidence_required": "pristine/effective source hashes and side-effect sentinel logs",
    },
)


def write_crosswalk(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(ROWS)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    write_crosswalk(args.output)
    print(f"wrote {len(ROWS)} crosswalk rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
