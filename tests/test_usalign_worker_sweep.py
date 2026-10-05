from benchmark.scripts.run_usalign_worker_sweep import (
    build_usalign_command,
    load_tasks,
    parse_worker_list,
)


def test_build_usalign_command_keeps_explicit_fast_option_and_matrix_output():
    command = build_usalign_command(
        "/opt/USalign",
        "/tmp/query.pdb",
        "/tmp/reference.pdb",
        "/tmp/matrix.out",
        fast=True,
    )

    assert command == [
        "/opt/USalign",
        "/tmp/query.pdb",
        "/tmp/reference.pdb",
        "-fast",
        "-outfmt",
        "-1",
        "-m",
        "/tmp/matrix.out",
    ]


def test_parse_worker_list_deduplicates_and_preserves_numeric_order():
    assert parse_worker_list("8,1,4,8,16") == [1, 4, 8, 16]


def test_load_tasks_limits_templates_before_expanding_chain_sides(tmp_path):
    template_list = tmp_path / "templates.txt"
    template_list.write_text("3i6eEF\n3hpgAF\n")
    interface_root = tmp_path / "interfaces"
    interface_root.mkdir()
    query = tmp_path / "query.pdb"
    query.write_text("query")

    tasks = load_tasks(template_list, interface_root, query, template_limit=1)

    assert [(task["template_id"], task["chain"]) for task in tasks] == [
        ("3i6eEF", "E"),
        ("3i6eEF", "F"),
    ]
