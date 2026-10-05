from pathlib import Path

from tools.build_project_chronology import run_evidence


def test_run_evidence_reads_bounded_nested_results(tmp_path: Path) -> None:
    run_dir = tmp_path / "tmp" / "agent" / "20260726-ranked-smoke"
    job_dir = run_dir / "runs" / "1392725"
    copied_source = job_dir / "baseline" / "src"
    copied_source.mkdir(parents=True)
    (copied_source / "ignored.py").write_text("pass\n", encoding="utf-8")
    (job_dir / "results.tsv").write_text(
        "arm\treturn_code\tselected_pairs\n"
        "baseline\t0\t2\n"
        "ranked\t0\t1\n",
        encoding="utf-8",
    )

    status, evidence, excluded = run_evidence(run_dir, tmp_path)

    assert status == "recorded:success"
    assert evidence == (
        "tmp/agent/20260726-ranked-smoke/runs/1392725/results.tsv",
    )
    assert excluded == 0


def test_run_evidence_reports_nested_nonzero_return(tmp_path: Path) -> None:
    run_dir = tmp_path / "tmp" / "agent" / "20260726-ranked-smoke"
    job_dir = run_dir / "runs" / "1392726"
    job_dir.mkdir(parents=True)
    (job_dir / "results.tsv").write_text(
        "arm\treturn_code\n"
        "baseline\t0\n"
        "ranked\t1\n",
        encoding="utf-8",
    )

    status, _, _ = run_evidence(run_dir, tmp_path)

    assert status == "recorded:exit_1"
