from pathlib import Path


def test_reranker_slurm_has_no_node_constraint():
    text = Path("benchmark/jobs/train_reranker.sbatch").read_text()
    assert "--nodelist" not in text
    assert "--cpus-per-task=2" in text
    assert "train_reranker.py" in text
