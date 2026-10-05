import csv
import json

from benchmark.scripts.audit_bijective_scores import audit, sha256_file


def write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_cross_only_score_allows_missing_native_interface_and_auxiliary_irmsd_failure(tmp_path):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    raw = tmp_path / "combined.json"
    component = tmp_path / "component.json"
    for path in (model, native):
        path.write_text("ATOM\n", encoding="ascii")
    component.write_text('{"GlobalDockQ": 0.4}\n', encoding="utf-8")
    raw.write_text('{"evaluation_mode": "pairwise_cross_fallback", "GlobalDockQ": null}\n', encoding="utf-8")
    components = [
        {"interface": "AB", "status": "scored", "raw_json": str(component), "raw_json_sha256": sha256_file(component)},
        {"interface": "AC", "status": "no_native_interface", "error": "no native interface"},
    ]
    models = tmp_path / "models.tsv"
    interfaces = tmp_path / "interfaces.tsv"
    write_tsv(models, [{
        "dataset_row_id": "difficult:1", "score_status": "scored",
        "score_scope": "requested_cross_interfaces_only", "source_gate_status": "strict_clean",
        "source_model_path": str(model), "staged_model_path": str(model),
        "source_model_sha256": sha256_file(model), "native_pdb_path": str(native),
        "native_pdb_sha256": sha256_file(native), "raw_dockq_json": str(raw),
        "raw_dockq_json_sha256": sha256_file(raw), "model_receptor_chains": "A",
        "model_ligand_chains": "BC", "native_receptor_chains": "A",
        "native_ligand_chains": "BC", "dockq_mapping": "ABC:ABC",
        "dockq_global": "", "dockq_global_status": "unavailable_cross_only",
        "dockq_component_runs": json.dumps(components), "irmsd_status": "failed_auxiliary",
        "dockq_version": "2.1.3", "mapping_validation_status": "aligned_default",
    }])
    write_tsv(interfaces, [{
        "staged_model_path": str(model), "record_type": "interface",
        "interface": "AB", "requested_cross_interface": "True",
    }])

    report, failures = audit(models, interfaces)

    assert failures == []
    assert report["audit_status"] == "passed"
    assert report["failed_auxiliary_irmsd_rows"] == 1
    assert report["score_scope"] == {"requested_cross_interfaces_only": 1}
