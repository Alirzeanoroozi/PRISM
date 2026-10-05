from pathlib import Path


def test_launcher_records_signals_and_delegates_status_classification():
    script = Path("benchmark/scripts/submit_comparison_batches.sbatch").read_text(encoding="utf-8")

    assert 'trap on_term TERM' in script
    assert 'trap on_int INT' in script
    assert 'PRISM_STAGE_STATUS_PATH="$run_root/status/stages.jsonl"' in script
    assert 'pipeline_returned.json' in script
    assert 'pipeline_completion_contract.py' in script
    assert 'TEMPLATE_LIST' in script
    assert 'cp "$TEMPLATE_LIST" templates/calculated_templates.txt' in script
    assert "printf 'filter_mode\\t%s\\n'" in script
    assert "printf 'filter_asset_root\\t%s\\n'" in script
