import prism


def test_prism_compatible_options_are_accepted():
    args = prism.build_parser().parse_args([
        "--inputs_csv", "pairs.csv",
        "--generate_templates",
        "--template_limit", "100",
        "--surface-backend", "freesasa",
        "--gtalign-dev-min-length", "4",
        "--gtalign-pre-score", "0.2",
        "--gtalign-speed", "2",
        "--gtalign-refinement", "5",
        "--multiprot-workers", "3",
        "--multiprot-path", "/tools/multiprot",
        "--refine",
        "--pyrosetta-output-dir", "pyro-out",
        "--pyrosetta-init-options", "-mute all",
        "--fiberdock-dir", "/tools/fiberdock",
        "--dockq-no-align",
    ])

    assert args.inputs_csv == "pairs.csv"
    assert args.generate_templates is True
    assert args.template_limit == 100
    assert args.surface_backend == "freesasa"
    assert args.gtalign_dev_min_length == 4
    assert args.multiprot_workers == 3
    assert args.multiprot_path == "/tools/multiprot"
    assert args.refine is True
    assert args.pyrosetta_output_dir == "pyro-out"
    assert args.fiberdock_dir == "/tools/fiberdock"
    assert args.dockq_no_align is True


def test_legacy_spellings_and_defaults_remain_supported():
    args = prism.build_parser().parse_args([
        "--generate_templates", "false",
        "--template-limit", "7",
        "--surface_backend", "naccess",
        "--gtalign_dev_min_length", "3",
        "--dockq-no-align", "false",
    ])

    assert args.generate_templates is False
    assert args.template_limit == 7
    assert args.surface_backend == "naccess"
    assert args.dockq_no_align is False
    assert args.refine is True


def test_no_refine_explicitly_skips_prescript_default():
    assert prism.build_parser().parse_args(["--no-refine"]).refine is False


def test_usalign_provider_options_are_explicit():
    args = prism.build_parser().parse_args([
        "--aligner", "usalign",
        "--usalign-path", "/tools/USalign",
        "--usalign-fast",
    ])

    assert args.aligner == "usalign"
    assert args.usalign_path == "/tools/USalign"
    assert args.usalign_fast is True
