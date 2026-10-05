# VALAR local-Qwen Copilot instructions

These instructions define Copilot behavior for the VALAR framework and for
autonomous sessions backed by a local Qwen server. Follow the repository
`AGENTS.md`, selected-project instructions, active plan, and bounded objective.
The current Copilot CLI dispatches that objective with `/autopilot` plus
`--experimental`; `/goal` is the legacy VALAR name.
Higher-level safety and scope constraints cannot be weakened by project notes.

## Detect the operating mode

- Local-Qwen worker mode is active when VALAR supplies a run directory, persisted session ID, and localhost provider metadata.
- Framework-maintenance mode applies when editing this repository outside a worker run. Preserve existing architecture, use tests, and do not launch project work unless explicitly requested.
- Never infer or select a project. The manifest and prompt identify the only approved project and workstream.
- The main orchestrator, not the local-Qwen worker, owns hourly supervision and workflow changes. The orchestrator uses GPT-5.6 Luna XHigh for operational review and GPT-5.6 Sol Medium for architecture, decomposition, or adjudication when warranted.

## Main-orchestrator supervision

During a long-running local-Qwen session, the main orchestrator inspects durable state at launch, hourly, at checkpoint/status changes, and at walltime warnings. Review:

- `manifest.json`, external provider/job/session metadata, and reservation state;
- `checkpoint.md`, `handoff.md`, `report.md`, `evidence/`, and `validation/`;
- bounded tails or structured records from `logs/slurm.log`, `logs/qwen-server.log`, `logs/provider-smoke.json`, `logs/copilot.log`, `logs/copilot-output.jsonl`, and `logs/copilot-help-providers.txt`, plus `evidence/copilot-runtime.json`;
- relevant selected-project Git status and protected-path changes.

Record a timestamped `CONTINUE`, `NO_CHANGE`, `ADJUST_PROMPT`, `ADJUST_WORKFLOW`, `CHECKPOINT_AND_RESUME`, `BLOCKED`, or `STOP` decision. The orchestrator may
clarify/narrow the objective, reorder approved steps, select a relevant
skill/MCP, adjust a bounded budget, or resume the same session. It must not
switch projects, broaden scope, overwrite `prompt.md`, change validation
criteria silently, or duplicate jobs.

After restart, reconcile manifests, reservations, Slurm status, sessions, checkpoints, handoffs, and evidence, then perform overdue reviews before new scheduling decisions. Hourly supervision must be evidenced by durable timestamps and decisions; it does not imply a background daemon exists.

## Project-scoped local reports and logs

- The controller supplies `RUN_DIR` as
  `/home/rshadi25/valar-agent-framework/evidence/<project-slug>/<run-id>/`.
- Use only `RUN_DIR/report.md`, `RUN_DIR/logs/`, `RUN_DIR/evidence/`,
  `RUN_DIR/validation/`, and `RUN_DIR/evidence/index.json` for project reports,
  logs, evidence, validation, and canonical path metadata.
- Never write project outputs to the shared framework `evidence/`, `logs/`, or
  repository root. Provider logs stay in `RUN_DIR/logs/`; summarize them in the
  report and do not duplicate them across projects.
- Reference these paths in the durable handoff/report. If daily-organizer sync
  fails, keep the fallback summary and blocker in `RUN_DIR`; the organizer is
  only a compact index.

## Project-specific orchestrator audit log

- The main orchestrator owns a concise daily audit entry for every selected project after meaningful execution, hourly review, checkpoint/status change, and closeout. The worker should provide the paths and facts needed for this entry; it should not write to a global organizer checkout unless explicitly assigned that path.
- Append to `/home/rshadi25/GitHub/daily-organizer-life-os/projects/research/<project-slug>/daily-log/YY-MM-DD-HH.md`, using the established project slug and local project timezone. For example, `prism-histone` uses `prism_histone`.
- Include timestamp, project/run/workstream, backend/model and job/session metadata, processed files or bounded directories, actions/commands, achieved results, evidence/validation paths, decisions, blockers, next action, and branch/worktree. Append within the same hour; never overwrite an existing entry.
- Keep this as a compact tracking index. Do not copy full prompts/logs, raw data, secrets, or unverified claims. The durable run directory remains authoritative. If the organizer checkout is unavailable, preserve the summary under the selected project's local evidence namespace and record the sync blocker rather than claiming completion.
- After explicit user authorization, the main orchestrator may run `LOCAL_ORGANIZER_ROOT=/home/rshadi25/GitHub/daily-organizer-life-os /home/rshadi25/bin/sync_research_memory.sh --commit-message "Sync <project-slug> daily log"`; verify commit/push output. Never push the selected research repository.

## Explicit no-Copilot fallback

If the user explicitly requests that Copilot/local Qwen not be used:

- Do not start `qwen-llm` or Copilot and do not silently substitute this path because a Qwen worker is slow or failed.
- The main orchestrator may initiate GPT-5.6 Luna XHigh Codex subagents for approved work. Use `repo_analyst`, `engineer`, `debugger`, and `science_researcher` according to task type.
- Use GPT-5.6 Sol Medium for decomposition, architecture, conflict resolution, or adjudication, not routine implementation or repetitive log review.
- Preserve the selected project, bounded goal, plan/memory recovery, evidence contract, checkpoint/handoff, validation, and protected-path rules.
- Independent subagents must have separate durable run records and isolated worktrees when they may write; use no more than two concurrent workers globally and never merge or promote automatically.
- Record the backend/model change and rationale in durable run metadata and `decisions.jsonl`.

## Fail-closed local provider

Before project work, the worker/controller must establish this sequence:

```text
Slurm allocation
-> qwen-llm run in the same allocation
-> localhost health succeeds
-> /v1/models returns an actual model ID
-> one bounded /v1/chat/completions provider smoke succeeds
-> Copilot executable/version and provider help are recorded
-> Copilot BYOK variables point to that endpoint/model
-> provider metadata is persisted
-> non-interactive objective dispatch succeeds
```

- Require `COPILOT_PROVIDER_TYPE=openai` and a `COPILOT_PROVIDER_BASE_URL` under `http://127.0.0.1:<port>/v1`.
- Obtain the wire model from `/v1/models`; do not guess it. Set the Copilot model variables to the verified local model.
- Refuse to start or continue project work if health, model discovery, resource inventory, provider binding, or objective activation evidence is missing.
- Never fall back to GitHub-hosted or another remote model. Do not use `--model auto` in local-Qwen worker mode.
- Never expose the server beyond `127.0.0.1`; do not tunnel it.
- Read GPU type/count, quantization, context, build, profile, port, job ID, and limits from the manifest/registry. Do not hard-code dynamic hardware policy here.

Useful provider and resource checks:

```bash
qwen-llm health --port "$PORT"
qwen-llm models --port "$PORT" --json
curl --fail --silent --show-error -H 'Content-Type: application/json' \
  -d '{"model":"<discovered-id>","messages":[{"role":"user","content":"Reply with exactly VALAR_PROVIDER_SMOKE_OK."}],"temperature":0,"max_tokens":32,"stream":false,"chat_template_kwargs":{"enable_thinking":false}}' \
  "http://127.0.0.1:$PORT/v1/chat/completions"
copilot help providers
copilot plugins list --kind mcp --kind skill --kind instruction --json
```

Do not print environment secrets or credentials. A local provider should not
require copying a real API key into prompts, logs, evidence, or manifests.

## Persistent objective and session

- Generate and persist a UUID `copilot_session_id` before first launch. A different concurrent run must have a different session ID.
- Start the persistent session with the objective command advertised by
  `copilot help commands` through Copilot's non-interactive `--prompt` mode.
  The current runtime requires `/autopilot` and `--experimental`; older
  runtimes may advertise `/goal`. This is required for offline batch workers because
  interactive `-i` dispatch is gated by Copilot login state in the installed
  CLI. The session ID remains durable and is reused for recovery.
- Use the exact durable `prompt.md` as the initial objective and preserve its SHA-256/provenance.
- Resume the same session and objective after interruption when safe. If replacement is unavoidable, preserve previous session IDs and explain why.
- Record durable evidence that the objective activated before considering the worker operational.

Conceptual launch syntax:

```bash
copilot \
  --experimental \
  --session-id "$COPILOT_SESSION_ID" \
  --add-dir "$PROJECT_ROOT" \
  --add-dir "$RUN_DIR" \
  --log-dir "$RUN_DIR/logs/copilot" \
  --no-ask-user \
  --mode autopilot \
  --max-autopilot-continues "$VALAR_MAX_AUTOPILOT_CONTINUES" \
  --no-remote \
  --no-remote-export \
  --log-level debug \
  --output-format json \
  --prompt "/autopilot $(cat "$RUN_DIR/prompt.md")"
```

The worker captures the provider probe in `logs/provider-smoke.json`, the
resolved Copilot executable/version in `evidence/copilot-runtime.json`, JSONL
stdout in `logs/copilot-output.jsonl`, and Copilot diagnostics in
`logs/copilot.log`. It must record the persisted active/completed objective
artifact (the runtime-specific `*objective*.json` under the session), matching
the exact run prompt and run ID, plus a corroborating
prompt-dispatch event tied to the persisted session ID and prompt SHA-256 in
`evidence/copilot-objective-activation.json` before reporting that the
objective activated. A successful process exit without that evidence is a
worker failure, not validation.

Read-only workers run Copilot inside an approved OS-enforced compute-node
boundary, with Bubblewrap preferred: the selected project and
manifest/prompt/artifact files are read-only, while only the assigned run
evidence/log/validation/checkpoint/handoff/report/decision paths and run-local
Copilot home are writable. An unavailable boundary is a fail-closed worker
error. The model must never edit `manifest.json` or the canonical project
checkout.

### Compute-node recovery contract

- The launcher sources `/etc/bashrc` before strict shell mode, loads required
  modules explicitly, and records the runtime and sandbox artifacts.
- If `bwrap` is unavailable, use the validated `singularity/4.3.2` fallback;
  keep the project read-only and bind only the approved run paths writable.
- Singularity warnings for absent `/etc/hosts`, `/etc/localtime`, or
  `/etc/resolv.conf` are expected with the minimal sandbox. Confirm the
  boundary artifact and continue only if the worker/provider checks pass.
- A `squeue` Munge/authentication failure is a host/session visibility problem:
  retry from `valar-login`, or use `sacct` and durable run state while marking
  live status unknown. Never infer job completion from missing `squeue` output.
- For project-specific scientific tools, use the project's documented absolute
  environment path and run its lightweight version/path preflight inside Slurm;
  do not depend on interactive `conda activate` state.

The controller supplies explicit tool/path/URL/MCP permissions appropriate to
the goal. Never broaden access merely to avoid a prompt. `--allow-all-paths`,
`--allow-all-urls`, and unrestricted destructive shell access are prohibited.

## Recover durable context first

Before substantial action:

1. Read the exact run `manifest.json`, `prompt.md`, `checkpoint.md`, `handoff.md`, `decisions.jsonl`, and available evidence/validation records.
2. Read selected-project `AGENTS.md` and Copilot instructions.
3. Use `project-memory` to recover canonical project state without treating memory as a raw log.
4. Read active Speckit `spec.md`, `plan.md`, `tasks.md`, and checklist when present; do not create a competing plan.
5. Inspect relevant existing outputs and Git status.
6. State the last validated step, current task, blockers, open hypotheses, expected visible output, and next concrete action.

Do not restart solved work, scan sibling run directories, mix project memories,
or replace an approved plan with an improvised one.

## Skill and MCP routing

Inventory the current installation instead of assuming names or availability.
Use the smallest relevant set from `/home/rshadi25/skill-details.md`:

- Recovery/orientation: `project-router`, `project-memory`, `graphify`.
- Planning: existing Speckit workflow; `exec-plan` for staged/risky work; `skills:brainstorming` only for unresolved design ambiguity.
- Implementation/debugging: `skills:test-driven-development`, then `skills:systematic-debugging` for failures.
- HPC/GPU: `general`; use `optimize-for-gpu` only for GPU-code optimization, not ordinary Qwen launching.
- Parallel work: `skills:dispatching-parallel-agents` and `skills:using-git-worktrees` only for genuinely independent, isolated tasks.
- Research/evidence: `skills:research-router-skill`, `scientific-critical-thinking`, `scientific-brainstorming`, `academic-research-suite`, `paper-lookup`, `database-lookup` as relevant.
- Source-grounded notebooks: `nlm-mcp` only when NotebookLM materially helps; enable the discovered server name rather than hard-coding one.
- Review/finish: `review-agent`, `skills:requesting-code-review`, `skills:receiving-code-review`, `skills:verification-before-completion`.
- Durable memory/cleanup: `memorize`, `pipeline-storage-hygiene`.

Do not invoke every skill or MCP. Record discovered resources, selected
capabilities, failures, and omitted resources in run evidence/provenance. If a
network-backed MCP is unavailable on a compute node, preserve local work and
write a precise dependency handoff; never fabricate research results.

## Autonomous execution loop

Continue within the approved workstream:

```text
orient -> inspect evidence -> choose next concrete step
-> compare alternatives/hypotheses -> select skills/tools/MCPs
-> execute -> inspect results -> test/experiment
-> interpret -> challenge -> fix/refine -> retest
-> validate evidence -> update approved task state -> continue
```

Do not stop merely because code was written, one command/test passed, an output
file exists, or an explanation sounds plausible. Stop only when success criteria
are evidenced, a genuine blocker/approval boundary is documented, or walltime
requires closeout.

## Evidence and output contract

- Tie claims to reproducible commands, tests, logs, files, figures, calculations, databases, or primary sources as appropriate.
- Separate `OBSERVATION`, `EVIDENCE`, `INFERENCE`, `HYPOTHESIS`, and `UNKNOWN` in scientific work.
- Report `IMPLEMENTED`, `TESTED`, `VALIDATED`, and `REVIEWED` separately. Never infer one from another.
- Preserve failures and negative results. Slurm `COMPLETED` proves only process exit, not goal success.
- Incomplete validation is an evidence gap, not an automatic instruction to stop an otherwise approved bounded pipeline. Record validation needs, evidence gaps, uncertainty, alternative explanations, and next informative checks in the checkpoint/handoff/report. Use `scientific-critical-thinking` for evidence assessment and `scientific-brainstorming` for testable alternatives; neither skill authorizes a `VALIDATED` claim.
- Write visible outputs only to approved project/worktree paths or the run's `artifacts/`, `evidence/`, `validation/`, and `logs/` directories.
- Keep large logs out of the manifest. Store paths, checksums, versions, commands, inputs/outputs, seeds, Git commit, and Slurm job ID.

## Files, Git, and project boundaries

- Never switch projects or broaden the approved goal.
- Treat raw datasets, curated references, validated results, model files, Qwen builds, environment/lock files, and canonical project memory as protected unless explicitly authorized.
- Editing workers use their assigned isolated worktree/branch. Read-only workers do not modify the canonical checkout.
- Do not merge, push, promote candidate results, resolve cross-worker conflicts, or delete worktrees automatically.
- Use atomic writes for important JSON; append decision records; never overwrite the initial `prompt.md`.
- Use `prompt-continuation-NN.md` for continuation/refocus and preserve source provenance.
- Do not expose secrets in commands, output, logs, commits, handoffs, or MCP requests.

## VALAR execution and time limits

- Login nodes are for lightweight inspection, editing, setup, and submission. Heavy tests, inference, scientific analysis, and GPU work require Slurm.
- Inspect `hostname`, `SLURM_JOB_ID`, `squeue`, relevant `sinfo`, and persisted reservations before resource decisions.
- One worker owns one Slurm job, one GPU allocation, one Qwen server, one localhost port, and one Copilot session. Use `--parallel 1`.
- The hard user ceiling is 8 GPUs; preserve the configured interactive reserve and admit at most two simultaneous GPU-backed Qwen workers globally across all projects. A plan may contain up to three useful independent workstreams, but a third Qwen lane remains deferred until a global slot is available. Request the least validated node footprint that supports the goal—prefer one A40 GPU or two V100 GPUs where the selected context is supported—and the minimum validated memory that reliably starts the model/provider. Never guess below a validated profile's memory requirement or request a larger footprint merely to increase concurrency.
- The two-worker limit is one global semaphore, not one limit per project. Production launches must use the framework-owned reservation root so concurrent project drivers cannot bypass it; caller-selected project-local roots are for tests or explicitly isolated simulations only.
- Never call `qwen-llm submit` from a framework worker allocation. The framework submits once; the job calls `qwen-llm run`.
- Before project work, persist `evidence/copilot-runtime.json`, the secret-free runtime-capabilities artifact containing Copilot version/path, provider and objective-command help, selected command/flags, Qwen model/endpoint, profile/context/build, Python environment, Git commit, Slurm allocation, and sandbox result. Select from observed runtime capabilities; do not assume `/goal`, `/autopilot`, or any other CLI feature.
- Read-only workers require an approved OS-enforced compute-node boundary. Bubblewrap is one supported implementation; an explicitly validated alternative is acceptable. If no approved boundary exists, fail closed and record an environment/escalation blocker. Never run unattended work unsandboxed, and use an allowlisted filesystem view that does not expose unrelated home-directory secrets through a broad read-only mount.
- For multi-case or scientifically consequential work, freeze and hash the deterministic script, configuration, input manifest, tool versions, seeds, resource limits, and expected output locations before Slurm execution. Use explicit `afterok`/`afterany` dependencies or bounded arrays with `%N`; do not create one LLM worker per scientific case or improvise a broad pipeline through repeated shell commands.
- Treat profile walltime as an explicit registry value and configure the walltime warning independently from status reporting, normally 300–600 seconds before walltime. Do not use a periodic 30-minute termination policy or claim a 2–4 hour runtime without a validated, authorized profile.

On walltime warning or termination:

1. Stop beginning major actions.
2. Preserve useful partial evidence and finish only small safe in-flight work.
3. Run the most informative bounded verification still possible.
4. Update approved task/plan state and project memory only with validated findings.
5. Write `checkpoint.md` and `handoff.md` with exact status, blocker, and next command/action.
6. Persist session/provider/job/profile/context/port metadata.
7. Stop Copilot and Qwen cleanly; do not report walltime as completion. If the closeout is signal-confirmed and resumable, record the run outcome as `CHECKPOINTED`, release the reservation, and resume the same run/session only within the approved retry budget; do not classify it as ordinary failure.

## Completion behavior

Before claiming success, invoke `skills:verification-before-completion`, run
fresh checks, inspect expected artifacts and statuses, and confirm no unintended
files changed. State unrun checks and residual risks. Do not begin another goal,
project, merge, or promotion automatically.
