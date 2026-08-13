import pytest
import sys

import prism
from src import prodigy_ranker
from src import transformation


def test_transformer_preserves_best_only_default(monkeypatch):
    candidates = [
        {"template": "best", "output_pdb": "best.pdb"},
        {"template": "second", "output_pdb": "second.pdb"},
    ]
    monkeypatch.setattr(transformation, "process_pair_for_template", lambda *args: candidates)

    assert transformation.transformer(["rec"], ["lig"]) == [
        ("rec", "lig", "best", "best.pdb")
    ]
    assert len(transformation.transformer(["rec"], ["lig"], return_all=True)) == 2


def test_parse_affinity():
    assert prodigy_ranker.parse_affinity("candidate  -9.373\n") == pytest.approx(-9.373)


def test_prodigy_selects_lowest_affinity_per_merged_candidate(monkeypatch, tmp_path):
    executable = tmp_path / "prodigy"
    executable.write_text("placeholder")
    candidates = [
        ("1abcA", "1abcB", "tpl1", "one.pdb"),
        ("1abcA", "1abcB", "tpl2", "two.pdb"),
    ]

    def score(path, receptor, ligand, **kwargs):
        affinity = -8.0 if path == "one.pdb" else -9.0
        return prodigy_ranker.ProdigyScore(path, "scored", affinity, 0, (), "", "", None)

    monkeypatch.setattr(prodigy_ranker, "score_candidate", score)
    selected = prodigy_ranker.select_top_candidates(
        candidates, top_k=1, executable=str(executable), output_dir=str(tmp_path)
    )

    assert selected == [candidates[1]]


def test_requested_ranking_fails_when_a_candidate_cannot_be_scored(monkeypatch, tmp_path):
    executable = tmp_path / "prodigy"
    executable.write_text("placeholder")
    candidates = [("1abcA", "1abcB", "tpl1", "one.pdb")]
    monkeypatch.setattr(
        prodigy_ranker,
        "score_candidate",
        lambda *args, **kwargs: prodigy_ranker.ProdigyScore(
            "one.pdb", "failed", None, None, (), "", "broken", 2
        ),
    )

    with pytest.raises(RuntimeError, match="PRODIGY failed"):
        prodigy_ranker.select_top_candidates(
            candidates, top_k=1, executable=str(executable), output_dir=str(tmp_path)
        )


def test_score_candidate_passes_merged_pdb_and_chain_selections(tmp_path):
    candidate = tmp_path / "candidate.pdb"
    candidate.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      2  CA  ALA B   1       4.000   0.000   0.000  1.00 20.00           C  \nEND\n"
    )
    fake = tmp_path / "fake_prodigy.py"
    fake.write_text(
        "import pathlib, sys\n"
        "args = sys.argv[1:]\n"
        "assert pathlib.Path(args[-4]).is_file(), args\n"
        "assert args[-3:] == ['--selection', 'A', 'B'], args\n"
        "print('-7.5')\n"
    )

    score = prodigy_ranker.score_candidate(
        str(candidate), "1abcA", "1abcB",
        executable=f"{sys.executable} {fake}", output_dir=str(tmp_path / "scores"),
    )

    assert score.status == "scored"
    assert score.affinity_kcal_mol == pytest.approx(-7.5)


def test_cli_defaults_keep_ranking_disabled():
    args = prism.build_parser().parse_args([])
    assert args.rank is False
    assert args.rank_method == "baseline"
    assert args.aligner == "tmalign"


def test_stage_events_are_opt_in_and_structured(monkeypatch, tmp_path):
    path = tmp_path / "stages.jsonl"
    monkeypatch.setenv("PRISM_STAGE_STATUS_PATH", str(path))

    assert prism.run_stage("ranking", lambda: "done") == "done"

    records = [__import__("json").loads(line) for line in path.read_text().splitlines()]
    assert [(row["stage"], row["event"]) for row in records] == [
        ("ranking", "started"),
        ("ranking", "completed"),
    ]
    assert records[-1]["return_code"] == 0
