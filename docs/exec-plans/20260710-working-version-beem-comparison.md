# Diagnose the working-version and current PRISM output gap

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` up to date as work proceeds.

## Purpose / Big Picture

Establish why `working_version/` produces a final model on another HPC while the current root pipeline can finish without writing a model or `.intRes.txt` file. Install the upstream BeEM `master` implementation at the project-local path expected by `working_version`, run bounded previews through Slurm, and change code only after a stage boundary and root cause are demonstrated. Success means both variants have reproducible stage manifests, explicit failure reasons, and either a generated final output or an evidence-backed explanation of the exact gate that prevented it.

## Progress

- [x] Read repository guidance, project memory, and prior smoke diagnoses.
- [x] Confirm execution host and establish that compute previews require Slurm.
- [x] Inventory the `working_version` entry point and its BeEM path contract.
- [x] Pin and install BeEM `master` under `working_version/external_tools/BeEM-master`.
- [x] Build a non-destructive preview harness with per-stage output checks.
- [x] Run matched minimal previews for `working_version` and the current pipeline.
- [x] Confirm root-cause hypotheses with minimal tests and controlled threshold previews.
- [x] Implement and verify only the confirmed corrections.
- [x] Update project memory with durable verified conclusions.

## Surprises & Discoveries

- Observation: `working_version` actively imports `structuralAlignmentTM.StructuralAligner` and `flexibleRefinementRosetta.FlexibleRefinement`; despite legacy names in `prism.ini`, this copy is a TMalign plus Rosetta variant, not the pure MultiProt plus FiberDock implementation.
  Evidence: `working_version/run_files/mainController.py` imports and invokes those classes.
- Observation: BeEM is only on the mmCIF fallback/template-download path and is resolved as `../../external_tools/BeEM-master/BeEM` from a job directory.
  Evidence: `working_version/run_files/pdbDownload.py` and `working_version/run_files/templateGenerator.py`.
- Observation: the current root pipeline can exit successfully with `Passed pairs 0`, creating an empty refinement-energy file but no model or `.intRes.txt` file.
  Evidence: `tmp/agent/20260702-prism-main-simple-rigid-diagnosis/DIAGNOSIS.md` and `tmp/agent/20260702-prism-version-smoke-validation/RESULTS.md`.
- Observation: network access to GitHub is blocked in the default sandbox, so fetching BeEM requires an approved network-enabled command.
  Evidence: `git ls-remote https://github.com/kad-ecoli/BeEM.git` failed with `Could not resolve host: github.com`.
- Observation: the working parser crashed on a blank matrix separator and did not support the bundled TMalign's tokenized `0,1,2` rows or modern TM-score labels.
  Evidence: Slurm job `1336973` and `working-preview-1336973.log`.
- Observation: chain-file repair alone did not recover the current rigid positive; `SCFFTHRESHOLD=1.4` still produced zero passed pairs, while the otherwise matched `5.0` run produced one pair and a final model.
  Evidence: Slurm jobs `1337002`, `1337003`, and `1337008`.

## Decision Log

- Decision: install BeEM project-locally at the exact path expected by the working code, and record its commit rather than modifying a global environment.
  Rationale: this preserves portability and avoids changing unrelated environments.
  Date/Author: 2026-07-10 / Codex
- Decision: treat the model PDB and `.intRes.txt` as final outputs; an empty `refinement_energies.txt` is diagnostic metadata, not successful model generation.
  Rationale: this matches the user's observed failure and the actual refinement writer.
  Date/Author: 2026-07-10 / Codex
- Decision: separate download/conversion validation from alignment/refinement validation.
  Rationale: BeEM cannot explain no-output cases where an ordinary PDB file already exists and the pipeline reaches alignment; conflating these stages would obscure causality.
  Date/Author: 2026-07-10 / Codex
- Decision: restore the current surface scaffold compatibility default to `5.0` and expose `PRISM_SCFF_THRESHOLD` for controlled alternatives.
  Rationale: the matched rigid-positive experiment isolated this threshold as the candidate-generation gate after chain handling was fixed.
  Date/Author: 2026-07-10 / Codex

## Outcomes & Retrospective

BeEM was installed and validated at upstream commit `4d71e6bf120312669859fc08b8e0c97918913779`. The working pipeline and the repaired current default both produced final model and `.intRes.txt` files for the rigid positive. The working model exported with `I_sc=-24.91`; the final current-default validation exported with `I_sc=-23.245`. The main root causes were packaging/path assumptions, an incompatible working TMalign parser, inconsistent current chain-qualified paths, and the current `1.4` scaffold threshold. The remaining limitation is that the current Python 3 downloader still lacks a BeEM-backed mmCIF fallback; this is distinct from the fixed rigid-case failure.

## Context and Orientation

The current pipeline starts at `prism.py`, writes intermediates under `processed/`, aligns target surfaces with template interfaces using `src/alignment.py`, filters candidates in `src/transformation.py`, and refines accepted pairs in `src/rosetta_refinement.py`. The copied working pipeline starts at `working_version/run_files/prism.py`, creates `working_version/jobs/<job_id>/`, and coordinates legacy-style modules through `working_version/run_files/mainController.py`. Its relative configuration expects shared assets under `working_version/pdb`, `working_version/template`, `working_version/external_tools`, and final Rosetta outputs under `working_version/rosetta_output_1/<job_id>`.

BeEM converts mmCIF files that cannot be represented directly as one legacy PDB file into PDB bundles. The working code then merges `*-bundle*.pdb` files. This is an input-representation dependency, not the aligner or refiner itself.

## Plan of Work

First pin and smoke-test BeEM independently with a small local mmCIF fixture or a staged download, recording generated bundle names and chain mappings. Next construct matched preview workspaces under `tmp/agent/20260710-working-version-beem-comparison/`, using one known rigid positive and the same template where both implementations can consume equivalent chain-qualified inputs. Add boundary checks for downloaded PDBs, extracted surfaces, alignment payloads, transformed pairs, Rosetta score files, final model files, and `.intRes.txt` files. Run lightweight setup on the current host and submit structural alignment/refinement to a minimal CPU Slurm job. Compare all thresholds, target identifiers, template interfaces, executable paths, command return codes, and output path assumptions. Only after the first divergent stage is shown, test one isolated hypothesis and implement a correction if necessary.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, fetch `https://github.com/kad-ecoli/BeEM.git` into `working_version/external_tools/BeEM-master`, then record `git rev-parse HEAD`, executable metadata, and `BeEM --help` or a bounded fixture conversion result.
2. From the same repository root, inspect and stage the known rigid positive inputs in `tmp/agent/20260710-working-version-beem-comparison/` without modifying source datasets.
3. Submit a CPU-only preview with `sbatch --partition=kutem --account=kutem --qos=kutem` after verifying live queue access; write Slurm logs beneath the preview artifact directory.
4. Summarize each stage in a TSV with `pipeline`, `case_id`, `stage`, `status`, `expected_path`, `observed_count`, and `reason`.
5. Run the smallest regression check that demonstrates the confirmed discrepancy before and after any source correction.

## Validation and Acceptance

BeEM acceptance requires a pinned upstream commit, an executable `working_version/external_tools/BeEM-master/BeEM`, and a successful bounded conversion that produces the bundle files expected by the merge code. Pipeline acceptance requires exact commands and logs for both variants, non-empty alignment artifacts, explicit accepted/rejected candidate counts, command return codes for Rosetta, and a final status that distinguishes `model_written`, `filtered_no_candidates`, `refinement_failed`, and `score_rejected`. A successful end-to-end preview must produce both a model PDB and its `.intRes.txt` file.

## Idempotence and Recovery

All preview work uses a dedicated `tmp/agent/20260710-working-version-beem-comparison/` tree. Reruns use new case/run subdirectories or explicit resume checks and do not overwrite benchmark data. If BeEM already exists, verify its origin and commit before reusing it; do not replace an unknown installation silently. Slurm jobs write isolated logs and can be resubmitted without deleting prior evidence. No cleanup occurs until outputs are classified, and retained diagnostic artifacts are listed in the final report.

## Artifacts and Notes

- Plan: `docs/exec-plans/20260710-working-version-beem-comparison.md`
- Prior current-pipeline diagnosis: `tmp/agent/20260702-prism-main-simple-rigid-diagnosis/DIAGNOSIS.md`
- Prior version smoke report: `tmp/agent/20260702-prism-version-smoke-validation/RESULTS.md`
- New preview root: `tmp/agent/20260710-working-version-beem-comparison/`
- Detailed results: `tmp/agent/20260710-working-version-beem-comparison/RESULTS.md`
- Stage table: `tmp/agent/20260710-working-version-beem-comparison/stage_summary.tsv`

## Interfaces and Dependencies

- BeEM upstream: `https://github.com/kad-ecoli/BeEM/`, branch `master`, commit `4d71e6bf120312669859fc08b8e0c97918913779`.
- Working BeEM executable contract: `working_version/external_tools/BeEM-master/BeEM <input.cif>`.
- Working pipeline interpreter: Python 2-compatible runtime; exact environment must be discovered before execution.
- Current pipeline interpreter: a Python 3 environment providing `pandas`, `numpy`, and Biopython; prior smoke used `/home/rshadi25/.conda/envs/gtalign_env/bin/python3.11`.
- Structural tools: TMalign, NACCESS, and Rosetta 2022.42 paths must be verified at run time.
- Compute placement: setup/network on the current host; alignment and refinement inside a Slurm CPU allocation.
