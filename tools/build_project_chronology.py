#!/usr/bin/env python3
"""Build a small, auditable, directed chronology graph for PRISM-prescript.

The source trees are read-only inputs.  Only the requested derived output
directory is rewritten, after checking its marker file.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import Path


DATE_PREFIX = re.compile(r"^(20\d{6})(?:[-_]|$)")
DATE_IN_TEXT = re.compile(r"\b(20\d{2})[-_/]?(\d{2})[-_/]?(\d{2})\b")
EXCLUDED_PARTS = {"site-packages", "python-site", "node_modules", "__pycache__", ".pytest_cache", "graphify-out"}
EVIDENCE_NAMES = {
    "exit.json", "run.log", "stdout.log", "stderr.log", "parameters.tsv",
    "inputs.csv", "runtime_manifest.json", "results.tsv",
}
EVIDENCE_SUFFIXES = {".md", ".json", ".jsonl", ".log", ".out", ".err", ".tsv", ".sbatch", ".sh"}


@dataclass(frozen=True)
class Event:
    event_id: str
    event_date: str
    date_source: str
    kind: str
    source: Path
    status: str
    evidence: tuple[str, ...]
    excluded: int


def iso_from_prefix(name: str) -> str | None:
    match = DATE_PREFIX.match(name)
    if not match:
        return None
    raw = match.group(1)
    try:
        return date(int(raw[:4]), int(raw[4:6]), int(raw[6:])).isoformat()
    except ValueError:
        return None


def iso_from_text(text: str) -> str | None:
    match = DATE_IN_TEXT.search(text)
    if not match:
        return None
    try:
        return date(*map(int, match.groups())).isoformat()
    except ValueError:
        return None


def safe_id(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def is_excluded(path: Path) -> bool:
    return any(part in EXCLUDED_PARTS or part.startswith("asa_tools_") for part in path.parts)


def json_status(path: Path) -> str | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return "unreadable_status"
    if not isinstance(value, dict):
        return "unstructured_status"
    for key in ("status", "state", "outcome", "result"):
        observed = value.get(key)
        if isinstance(observed, str) and observed.strip():
            return f"recorded:{observed.strip().lower()}"
    for key in ("exit_code", "returncode", "return_code"):
        observed = value.get(key)
        if isinstance(observed, int):
            return "recorded:success" if observed == 0 else f"recorded:exit_{observed}"
    return "unstructured_status"


def results_tsv_status(path: Path) -> str | None:
    """Return terminal status from a paired-run results table."""
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeDecodeError, csv.Error):
        return "unreadable_status"
    if not rows or "return_code" not in rows[0]:
        return "unstructured_status"
    try:
        codes = [int(row["return_code"]) for row in rows]
    except (KeyError, TypeError, ValueError):
        return "unstructured_status"
    if all(code == 0 for code in codes):
        return "recorded:success"
    nonzero = sorted({code for code in codes if code != 0})
    if len(nonzero) == 1:
        return f"recorded:exit_{nonzero[0]}"
    return "mixed_recorded_status"


def run_evidence(run_dir: Path, root: Path) -> tuple[str, tuple[str, ...], int]:
    kept: list[str] = []
    states: list[str] = []
    excluded = 0
    # Only inspect the run root and named scheduler/status roots.  Some runs
    # contain cloned repos or full environments, so recursive discovery would
    # turn a chronology index into another dependency scan.
    roots = [run_dir] + [run_dir / name for name in ("status", "ai", "cosbi", "kutem", "slurm")]
    # Paired experiments retain each immutable submission beneath
    # runs/<job-id>/. Inspect that one bounded level for its aggregate table,
    # without traversing copied source trees or environments below it.
    nested_runs = run_dir / "runs"
    if nested_runs.is_dir():
        roots.extend(item for item in sorted(nested_runs.iterdir()) if item.is_dir())
    for candidate_root in roots:
        if not candidate_root.is_dir() or is_excluded(candidate_root.relative_to(root)):
            continue
        for item in sorted(candidate_root.iterdir()):
            if item.is_dir():
                if is_excluded(item.relative_to(root)):
                    excluded += 1
                continue
            rel = item.relative_to(root)
            if item.name == "exit.json":
                state = json_status(item)
                if state:
                    states.append(state)
            elif item.name == "results.tsv":
                state = results_tsv_status(item)
                if state:
                    states.append(state)
            if item.name in EVIDENCE_NAMES or item.suffix.lower() in EVIDENCE_SUFFIXES:
                if len(kept) < 40:
                    kept.append(rel.as_posix())
    status_set = set(states)
    if not states:
        status = "no_terminal_status"
    elif len(status_set) == 1:
        status = states[-1]
    else:
        status = "mixed_recorded_status"
    return status, tuple(kept), excluded


def add_doc_event(
    events: list[Event], root: Path, path: Path, kind: str,
    fixed_date: str | None = None, fixed_date_source: str | None = None,
) -> None:
    rel = path.relative_to(root)
    observed = fixed_date or iso_from_prefix(path.name) or iso_from_prefix(path.parent.name)
    if not observed:
        try:
            observed = iso_from_text(path.read_text(encoding="utf-8", errors="ignore")[:8000])
        except OSError:
            observed = None
    if not observed:
        return
    events.append(Event(
        event_id=f"{kind}_{safe_id(rel.as_posix())}", event_date=observed,
        date_source=fixed_date_source or "filename_or_heading", kind=kind, source=rel,
        status="documented", evidence=(rel.as_posix(),), excluded=0,
    ))


def collect_events(root: Path) -> list[Event]:
    events: list[Event] = []
    agent_root = root / "tmp" / "agent"
    if agent_root.is_dir():
        for run_dir in sorted(agent_root.iterdir()):
            if not run_dir.is_dir() or is_excluded(run_dir.relative_to(root)):
                continue
            observed = iso_from_prefix(run_dir.name)
            if not observed:
                continue
            status, evidence, excluded = run_evidence(run_dir, root)
            events.append(Event(
                event_id=f"run_{safe_id(run_dir.name)}", event_date=observed,
                date_source="tmp/agent directory prefix", kind="run",
                source=run_dir.relative_to(root), status=status,
                evidence=evidence, excluded=excluded,
            ))
    for path in sorted((root / "docs" / "exec-plans").glob("*.md")):
        add_doc_event(events, root, path, "exec_plan")
    for path in sorted((root / "docs" / "memory-archive").rglob("*.md")):
        add_doc_event(events, root, path, "memory_archive")
    for name in ("summary.md", "decisions.md", "open_questions.md"):
        path = root / ".agents" / "skills" / "project-memory" / "references" / name
        if path.exists():
            add_doc_event(
                events, root, path, "active_memory",
                fixed_date=date.fromtimestamp(path.stat().st_mtime).isoformat(),
                fixed_date_source="active-memory file mtime",
            )
    return sorted(events, key=lambda event: (event.event_date, event.kind, event.event_id))


def write_corpus(root: Path, out: Path, events: list[Event]) -> None:
    marker = out / ".chronology-generated"
    if out.exists() and not marker.exists():
        raise RuntimeError(f"refusing to overwrite unmarked directory: {out}")
    if out.exists():
        for child in (out / "events", out / "graphify-out"):
            if child.exists():
                shutil.rmtree(child)
    (out / "events").mkdir(parents=True, exist_ok=True)
    marker.write_text("Derived by tools/build_project_chronology.py; safe to regenerate.\n", encoding="utf-8")
    with (out / "manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["event_id", "date", "date_source", "kind", "source", "status", "evidence_files", "excluded_files"])
        for event in events:
            writer.writerow([event.event_id, event.event_date, event.date_source, event.kind, event.source.as_posix(), event.status, len(event.evidence), event.excluded])
            note = [
                f"# {event.event_date}: {event.source.name}", "",
                f"- Event ID: `{event.event_id}`",
                f"- Kind: `{event.kind}`",
                f"- Source: `{event.source.as_posix()}`",
                f"- Date evidence: `{event.date_source}`",
                f"- Recorded status: `{event.status}`",
                f"- Excluded vendor/cache files under this source: {event.excluded}", "",
                "## Selected evidence paths", "",
            ]
            note.extend(f"- `{path}`" for path in event.evidence) if event.evidence else note.append("- No selected evidence files.")
            (out / "events" / f"{event.event_id}.md").write_text("\n".join(note) + "\n", encoding="utf-8")
    index = ["# PRISM-prescript chronology corpus", "", "This derived corpus preserves source paths; it does not copy or alter source evidence.", "", "## Event order", ""]
    index.extend(f"- {event.event_date} — [{event.kind}]({event.source.as_posix()}) — `{event.status}`" for event in events)
    (out / "README.md").write_text("\n".join(index) + "\n", encoding="utf-8")


def build_graph(root: Path, out: Path, events: list[Event]) -> tuple[int, int]:
    from graphify.analyze import god_nodes, surprising_connections, suggest_questions
    from graphify.build import build_from_json
    from graphify.cluster import cluster, score_all
    from graphify.export import to_json
    from graphify.report import generate

    nodes, edges = [], []
    seen_dates, seen_statuses = set(), set()
    for event in events:
        event_file = f"events/{event.event_id}.md"
        event_node = f"events_{event.event_id}"
        nodes.append({"id": event_node, "label": event.source.name, "file_type": "rationale", "source_file": event_file, "source_location": "L1", "_origin": "ast", "metadata": {"date": event.event_date, "kind": event.kind, "status": event.status, "source": event.source.as_posix()}})
        date_id = f"date_{event.event_date.replace('-', '_')}"
        if date_id not in seen_dates:
            nodes.append({"id": date_id, "label": event.event_date, "file_type": "concept", "source_file": "README.md", "source_location": "L6", "_origin": "ast", "metadata": {"date": event.event_date}})
            seen_dates.add(date_id)
        edges.append({"source": date_id, "target": event_node, "relation": "occurred_on", "confidence": "EXTRACTED", "confidence_score": 1.0, "source_file": event_file, "source_location": "L1"})
        status_id = f"status_{safe_id(event.status)}"
        if status_id not in seen_statuses:
            nodes.append({"id": status_id, "label": event.status, "file_type": "concept", "source_file": "README.md", "source_location": "L6", "_origin": "ast"})
            seen_statuses.add(status_id)
        edges.append({"source": event_node, "target": status_id, "relation": "has_status", "confidence": "EXTRACTED", "confidence_score": 1.0, "source_file": event_file, "source_location": "L1"})
    latest_date = max(event.event_date for event in events)
    latest_date_id = f"date_{latest_date.replace('-', '_')}"
    nodes.append({"id": "readme_current_project_chronology", "label": f"Current project chronology ({latest_date})", "file_type": "rationale", "source_file": "README.md", "source_location": "L5", "_origin": "ast", "metadata": {"as_of": latest_date}})
    edges.append({"source": "readme_current_project_chronology", "target": latest_date_id, "relation": "current_as_of", "confidence": "EXTRACTED", "confidence_score": 1.0, "source_file": "README.md", "source_location": "L5"})
    # A directory prefix proves a calendar day, not an order among same-day
    # events.  Chronology is therefore encoded between date nodes only.
    date_ids = [f"date_{value.replace('-', '_')}" for value in sorted({event.event_date for event in events})]
    for earlier, later in zip(date_ids, date_ids[1:]):
        edges.append({"source": earlier, "target": later, "relation": "before", "confidence": "EXTRACTED", "confidence_score": 1.0, "source_file": "README.md", "source_location": "L6"})
    extraction = {"nodes": nodes, "edges": edges, "hyperedges": [], "input_tokens": 0, "output_tokens": 0}
    graph_dir = out / "graphify-out"
    graph_dir.mkdir(parents=True, exist_ok=True)
    (graph_dir / "chronology_extraction.json").write_text(json.dumps(extraction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    graph = build_from_json(extraction, root=out, directed=True)
    communities = cluster(graph)
    cohesion = score_all(graph, communities)
    labels = {community: "Chronology events" for community in communities}
    analysis = {"communities": {str(k): v for k, v in communities.items()}, "cohesion": {str(k): v for k, v in cohesion.items()}, "gods": god_nodes(graph), "surprises": surprising_connections(graph, communities)}
    analysis["questions"] = suggest_questions(graph, communities, labels)
    to_json(graph, communities, str(graph_dir / "graph.json"), force=True, community_labels=labels)
    report = generate(graph, communities, cohesion, labels, analysis["gods"], analysis["surprises"], {"total_files": len(events), "total_words": 0, "files": {"document": []}}, {"input": 0, "output": 0}, str(out), suggested_questions=analysis["questions"])
    (graph_dir / "GRAPH_REPORT.md").write_text(report, encoding="utf-8")
    (graph_dir / ".graphify_analysis.json").write_text(json.dumps(analysis, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return graph.number_of_nodes(), graph.number_of_edges()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."), help="project root")
    parser.add_argument("--out", type=Path, default=Path("docs/chronology"), help="derived output directory")
    args = parser.parse_args()
    root = args.root.resolve()
    out = (root / args.out).resolve() if not args.out.is_absolute() else args.out.resolve()
    try:
        import graphify  # noqa: F401
    except ModuleNotFoundError:
        parser.error(
            "the selected interpreter lacks Graphify; rerun with "
            "$(cat graphify-out/.graphify_python)"
        )
    events = collect_events(root)
    if not events:
        raise RuntimeError("no dated project events found")
    write_corpus(root, out, events)
    nodes, edges = build_graph(root, out, events)
    print(json.dumps({"events": len(events), "nodes": nodes, "edges": edges, "output": str(out)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
