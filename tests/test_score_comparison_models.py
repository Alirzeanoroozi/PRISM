import csv
import sys
from pathlib import Path

import benchmark.scripts.score_comparison_models as scorer


def _write_table(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_score_manifest_propagates_cpu_and_json_status(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("ATOM\n", encoding="ascii")
    native.write_text("ATOM\n", encoding="ascii")
    manifest = tmp_path / "models.csv"
    _write_table(
        manifest,
        [
            {
                "pipeline": "tmalign_rosetta",
                "status": "ready",
                "model_path": str(model),
                "native_path": str(native),
                "model_receptor": "A",
                "model_ligand": "B",
                "native_receptor": "A",
                "native_ligand": "B",
            }
        ],
    )
    output = tmp_path / "scored.csv"
    commands = []

    def fake_run(command, **kwargs):
        commands.append(command)
        out_csv = Path(command[command.index("--out-csv") + 1])
        _write_table(
            out_csv,
            [
                {
                    "irmsd": "1.0",
                    "dockq": "",
                    "dockq_global": "",
                    "dockq_interface_count": "2",
                    "dockq_raw_json_path": "raw.json",
                    "dockq_raw_json_sha256": "hash",
                    "dockq_mapping": "AB:AB",
                    "mapping_validation_status": "aligned_default",
                    "mapping_validation_errors": "",
                    "dockq_json_status": "valid_unscored",
                    "error_irmsd": "",
                    "error_dockq": "missing GlobalDockQ",
                }
            ],
        )
        return type("Completed", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(scorer.subprocess, "run", fake_run)
    scorer.score_manifest(
        manifest,
        output,
        Path("python"),
        tmp_path,
        tmp_path / "raw",
        dockq_no_align=False,
        n_cpu=2,
    )

    assert commands[0][commands[0].index("--n-cpu") + 1] == "2"
    row = next(csv.DictReader(output.open(newline="", encoding="utf-8")))
    assert row["dockq_json_status"] == "valid_unscored"
