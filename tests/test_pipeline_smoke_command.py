from pathlib import Path


def test_smoke_uses_current_prism_cli_defaults():
    text = Path("benchmark/scripts/run_prism_pipeline_smoke.sh").read_text()
    assert "--template_list" not in text
    assert "--inputs_csv" not in text
    assert "--generate_templates" not in text
    assert "timeout \"$timeout_sec\" \"$python_bin\" prism.py" in text
    assert 'pipeline_args+=(--aligner "${PRISM_SMOKE_ALIGNER}")' in text
    assert 'pipeline_args+=(--no-refine)' in text
    assert 'status/pipeline_returned.json' in text
