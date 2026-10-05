import csv


def _write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_merge_baseline_shards_sorts_claims_and_artifacts(tmp_path):
    from benchmark.scripts.collect_pipeline_verification_baseline import merge_baseline_shards

    shard_a = tmp_path / "baseline-ai"
    shard_b = tmp_path / "baseline-cosbi"
    _write_tsv(
        shard_a / "claims.tsv",
        ["claim_id", "classification"],
        [{"claim_id": "v3", "classification": "unsupported"}],
    )
    _write_tsv(
        shard_a / "artifact_manifest.tsv",
        ["run_root", "relative_path", "sha256", "bytes"],
        [{"run_root": "v3", "relative_path": "status/exit.json", "sha256": "a", "bytes": "1"}],
    )
    _write_tsv(
        shard_b / "claims.tsv",
        ["claim_id", "classification"],
        [{"claim_id": "smoke", "classification": "supported"}],
    )
    _write_tsv(
        shard_b / "artifact_manifest.tsv",
        ["run_root", "relative_path", "sha256", "bytes"],
        [{"run_root": "smoke", "relative_path": "status/exit.json", "sha256": "b", "bytes": "2"}],
    )
    output = tmp_path / "baseline"

    merge_baseline_shards([shard_b, shard_a], output)

    with (output / "claims.tsv").open(newline="", encoding="utf-8") as handle:
        claims = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["claim_id"] for row in claims] == ["smoke", "v3"]
    with (output / "artifact_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        artifacts = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["run_root"] for row in artifacts] == ["smoke", "v3"]


def test_merge_baseline_shards_rejects_duplicate_claim_ids(tmp_path):
    from benchmark.scripts.collect_pipeline_verification_baseline import merge_baseline_shards

    for name in ("one", "two"):
        shard = tmp_path / name
        _write_tsv(shard / "claims.tsv", ["claim_id"], [{"claim_id": "duplicate"}])
        _write_tsv(shard / "artifact_manifest.tsv", ["run_root", "relative_path"], [])

    try:
        merge_baseline_shards([tmp_path / "one", tmp_path / "two"], tmp_path / "baseline")
    except ValueError as exc:
        assert "duplicate claim_id" in str(exc)
    else:
        raise AssertionError("duplicate baseline claim identifiers were accepted")
