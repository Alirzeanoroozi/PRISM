import csv
import hashlib

from benchmark.scripts.attach_native_labels import attach_labels


def test_attach_native_labels_uses_dockq_threshold(tmp_path):
    candidates = tmp_path / "candidates.csv"
    candidates.write_text(
        "model_complex,model_left,model_right,status\n"
        f"{tmp_path / 'complex.pdb'},{tmp_path / 'left.pdb'},{tmp_path / 'right.pdb'},generated\n"
    )
    scores = tmp_path / "scores.csv"
    scores.write_text(f"model_pdb,dockq,irmsd\n{tmp_path / 'complex.pdb'},0.25,1.5\n")
    output = tmp_path / "labeled.csv"
    assert attach_labels(candidates, scores, output) == 1
    row = next(csv.DictReader(output.open()))
    assert row["native_like"] == "1"
    assert row["label_status"] == "labeled"


def test_attach_native_labels_keeps_unlabeled_rows(tmp_path):
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("model_left,model_right,status\n,,alignment_failed\n")
    scores = tmp_path / "scores.csv"
    scores.write_text("model_pdb,dockq,irmsd\n")
    output = tmp_path / "labeled.csv"
    assert attach_labels(candidates, scores, output) == 0
    row = next(csv.DictReader(output.open()))
    assert row["label_status"] == "unlabeled"


def test_attach_native_labels_uses_cross_interface_metric_and_durable_identity(tmp_path):
    model = tmp_path / "complex.pdb"
    model.write_text("MODEL\n", encoding="ascii")
    model_hash = hashlib.sha256(model.read_bytes()).hexdigest()
    candidates = tmp_path / "candidates.csv"
    candidates.write_text(
        "dataset_row_id,model_complex,source_model_sha256,status\n"
        f"rigid:000001,{model},{model_hash},generated\n"
    )
    scores = tmp_path / "scores.tsv"
    scores.write_text(
        "dataset_row_id\tsource_model_path\tsource_model_sha256\tscore_status\tsource_gate_status\t"
        "dockq_global\tdockq_cross_mean\tirmsd_grouped_min\n"
        f"rigid:000001\t{model}\t{model_hash}\tscored\tstrict_clean\t0.95\t0.05\t12.5\n"
    )
    output = tmp_path / "labeled.csv"

    assert attach_labels(candidates, scores, output) == 1
    row = next(csv.DictReader(output.open()))
    assert row["dockq"] == "0.05000000"
    assert row["native_like"] == "0"
    assert row["label_metric"] == "dockq_cross_mean"
    assert row["label_model_sha256"] == model_hash


def test_attach_native_labels_accepts_explicit_cross_only_score(tmp_path):
    model = tmp_path / "complex.pdb"
    model.write_text("MODEL\n", encoding="ascii")
    model_hash = hashlib.sha256(model.read_bytes()).hexdigest()
    candidates = tmp_path / "candidates.csv"
    candidates.write_text(
        "dataset_row_id,model_complex,source_model_sha256,status\n"
        f"rigid:000001,{model},{model_hash},generated\n"
    )
    scores = tmp_path / "scores.tsv"
    scores.write_text(
        "dataset_row_id\tsource_model_path\tsource_model_sha256\tscore_status\t"
        "score_scope\tsource_gate_status\tdockq_cross_mean\tirmsd_grouped_min\n"
        f"rigid:000001\t{model}\t{model_hash}\tscored_cross_only\t"
        "requested_cross_interfaces_only\tstrict_clean\t0.30\t2.0\n"
    )
    output = tmp_path / "labeled.csv"

    assert attach_labels(candidates, scores, output) == 1
    row = next(csv.DictReader(output.open()))
    assert row["label_status"] == "labeled"
    assert row["label_score_scope"] == "requested_cross_interfaces_only"
