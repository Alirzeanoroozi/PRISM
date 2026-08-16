import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from src.pdb_download import pdb_downloader
from src.eda.analyse_pdbs import run_analysis
from src.template_generate import template_generator
from src.surface_extract import extract_surfaces
from src.alignment import align
from src.alignment_gtalign import align_gtalign
from src.transformation import transformer
from src.rosetta_refinement import refiner
from src.compare import compare_and_summarize


def parse_bool(value):
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y", "on"}:
        return True
    if normalized in {"false", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"invalid boolean value: {value}")


def record_stage_event(stage, event, return_code=None, detail=""):
    """Append an event only when PRISM_STAGE_STATUS_PATH is configured."""
    raw_path = os.environ.get("PRISM_STAGE_STATUS_PATH")
    if not raw_path:
        return
    record = {
        "stage": stage,
        "event": event,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "detail": detail,
    }
    if return_code is not None:
        record["return_code"] = return_code
    path = Path(raw_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def run_stage(stage, operation):
    record_stage_event(stage, "started")
    try:
        result = operation()
    except Exception as exc:
        record_stage_event(stage, "failed", 1, f"{type(exc).__name__}: {exc}")
        raise
    record_stage_event(stage, "completed", 0)
    return result


def _select_baseline_top_k(candidates, top_k):
    """Keep transformer score order within each receptor-ligand group."""
    selected = []
    counts = {}
    for candidate in candidates:
        key = (candidate[0], candidate[1])
        if counts.get(key, 0) < top_k:
            selected.append(candidate)
            counts[key] = counts.get(key, 0) + 1
    return selected


def run_alignment_stage(args, targets, templates):
    if args.aligner == "tmalign":
        return run_stage("alignment", lambda: align(targets, templates))
    if args.aligner == "gtalign":
        return run_stage("alignment", lambda: align_gtalign(
            targets,
            templates,
            gtalign_path=args.gtalign_path,
            dev_min_length=args.gtalign_dev_min_length,
            pre_score=args.gtalign_pre_score,
            speed=args.gtalign_speed,
            refinement=args.gtalign_refinement,
        ))

    os.environ["PRISM_MULTIPROT"] = args.multiprot_path
    from src.alignment_multiprot import align_multiprot

    return run_stage("alignment", lambda: align_multiprot(
        targets,
        templates,
        output_dir="processed/alignment",
        max_workers=args.multiprot_workers,
        multiprot_path=args.multiprot_path,
        multiprot_mode=args.multiprot_mode,
        multiprot_params=args.multiprot_params,
        multiprot_solutions=args.multiprot_solutions,
    ))


def run_ranking_stage(args, candidates):
    if not args.rank:
        return candidates
    if args.top_k < 1:
        raise ValueError("--top-k must be positive when --rank is enabled")

    print(f"[5a/6] Candidate ranking ({args.rank_method})")
    if args.rank_method == "baseline":
        selected = run_stage(
            "ranking", lambda: _select_baseline_top_k(candidates, args.top_k),
        )
    else:
        from src.prodigy_ranker import select_top_candidates

        selected = run_stage("ranking", lambda: select_top_candidates(
            candidates,
            top_k=args.top_k,
            executable=args.prodigy_executable,
            output_dir=args.prodigy_output_dir,
            distance_cutoff=args.prodigy_distance_cutoff,
            acc_threshold=args.prodigy_acc_threshold,
            temperature=args.prodigy_temperature,
            timeout=args.prodigy_timeout,
        ))
    print(f"  selected {len(selected)} candidate(s)")
    return selected


def run_refinement_stage(args, candidates):
    compare_pairs = list(candidates)
    if not args.refine or not candidates:
        return compare_pairs

    print(f"[6a/6] {args.refiner} refinement")
    if args.refiner == "external_rosetta":
        refined = run_stage("refinement", lambda: refiner(candidates))
        for rec, lig, isc, tsc, output_path in refined:
            print(
                f"  refined {rec} + {lig}: "
                f"int={isc} total={tsc} -> {output_path}"
            )
        if refined:
            compare_pairs = [
                (rec, lig, template, output_pdb)
                for (rec, lig, template, _), (_, _, _, _, output_pdb)
                in zip(candidates, refined)
            ]
        return compare_pairs

    if args.refiner == "pyrosetta":
        from src.pyrosetta_refinement import refine_merged_candidates

        return run_stage("refinement", lambda: refine_merged_candidates(
            candidates,
            output_root=args.pyrosetta_output_dir,
            init_options=args.pyrosetta_init_options,
        ))

    os.environ["PRISM_FIBERDOCK_DIR"] = args.fiberdock_dir
    from src.fiberdock_refinement import refine_merged_candidates

    return run_stage("refinement", lambda: refine_merged_candidates(candidates))


def main(args):
    if args.surface_backend != "freesasa":
        raise ValueError(f"Unsupported surface backend: {args.surface_backend}")
    print("[1/6] PDB download")
    receptor_targets, ligand_targets = run_stage("input", lambda: pdb_downloader(args))
    targets = sorted(set(receptor_targets + ligand_targets))
    for r, l in zip(receptor_targets, ligand_targets):
        print(f"  {r} -> {l}")

    if args.generate_templates:
        print("[2a/6] Template analysis & filtering")
        results, filtered_count, failed_count = run_analysis()
        print(f"  analysed={len(results)} filtered={filtered_count} failed={failed_count}")

        print("[2b/6] Template artifact generation")
        templates = template_generator()
        print(f"  generated {len(templates)} templates")
    else:
        with open("templates/calculated_templates.txt") as f:
            templates = [line.strip() for line in f if line.strip()]
        print(f"[2/6] Loaded {len(templates)} pre-computed templates")

    if args.template_limit > 0:
        templates = templates[: args.template_limit]
        print(f"  using top {len(templates)} templates")

    print("[3/6] Target surface extraction")
    failed_surfaces = extract_surfaces(targets)
    if failed_surfaces:
        print(f"  WARNING: surface extraction failed for {failed_surfaces} target(s)")

    print(f"[4/6] Structural alignment ({args.aligner})")
    run_alignment_stage(args, targets, templates)

    print("[5/6] Transformation + filtering")
    passed = run_stage(
        "transformation",
        lambda: transformer(receptor_targets, ligand_targets, return_all=args.rank),
    )
    print(f"  accepted {len(passed)} receptor-ligand candidates")
    for receptor, ligand, template, output in passed:
        print(f"  {receptor} + {ligand} via {template} -> {output}")

    passed = run_ranking_stage(args, passed)
    compare_pairs = run_refinement_stage(args, passed)

    print("[6/6] Compare outputs vs native + DockQ")
    if compare_pairs:
        summary_csv, rows = compare_and_summarize(
            compare_pairs,
            dockq_no_align=args.dockq_no_align,
            n_jobs=args.compare_jobs,
        )
        print(f"  summary -> {summary_csv}")
        for row in rows:
            dockq = row["dockq"] if row["dockq"] is not None else "NA"
            rmsd = row["receptor_ca_rmsd"] if row["receptor_ca_rmsd"] is not None else "NA"
            print(
                f"  {row['template']}: output={row['output_pdb']} "
                f"receptor_RMSD={rmsd} DockQ={dockq}"
            )
    else:
        print("  no accepted outputs to compare")


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs_csv", type=str, default="inputs.csv")
    parser.add_argument("--generate_templates", action="store_true", help="Run template analysis & artifact generation before alignment.")
    parser.add_argument("--aligner", choices=["tmalign", "gtalign", "multiprot"], default="tmalign")
    parser.add_argument("--gtalign_path", default="gtalign")
    parser.add_argument("--gtalign-dev-min-length", type=int, default=3)
    parser.add_argument("--gtalign-pre-score", type=float, default=0.0)
    parser.add_argument("--gtalign-speed", type=int, default=0)
    parser.add_argument("--gtalign-refinement", type=int, default=3)
    parser.add_argument("--multiprot-workers", type=int, default=8)
    parser.add_argument(
        "--multiprot-path",
        default=os.environ.get("PRISM_MULTIPROT", "external_tools/multiprot.Linux"),
    )
    parser.add_argument(
        "--multiprot-mode",
        choices=["current", "legacy_compatible"],
        default="current",
        help="Use the current Kabsch adapter or preserve legacy solver transforms.",
    )
    parser.add_argument(
        "--multiprot-params",
        help="Optional MultiProt params.txt copied into each isolated worker directory.",
    )
    parser.add_argument(
        "--multiprot-solutions",
        type=int,
        default=3,
        help="Maximum legacy-compatible solutions retained per query/interface alignment.",
    )
    parser.add_argument(
        "--surface-backend",
        choices=["freesasa"],
        default="freesasa",
        help="Surface backend used by this PRISM implementation.",
    )
    parser.add_argument("--template_limit", type=int, default=100, help="Limit number of templates aligned (0 = all)")
    parser.add_argument("--refine", action="store_true", help="Run Rosetta refinement on accepted candidates")
    parser.add_argument(
        "--refiner",
        choices=["external_rosetta", "pyrosetta", "fiberdock"],
        default="external_rosetta",
    )
    parser.add_argument("--pyrosetta-output-dir", default="processed/pyrosetta_refinement")
    parser.add_argument(
        "--pyrosetta-init-options",
        default="-mute all -constant_seed -jran 12345",
    )
    parser.add_argument(
        "--fiberdock-dir",
        default=os.environ.get("PRISM_FIBERDOCK_DIR", "external_tools/fiberdock"),
    )
    parser.add_argument(
        "--dockq-no-align",
        action="store_true",
        help="Pass --no_align to DockQ in the final compare step (default: align)",
    )
    parser.add_argument(
        "--compare-jobs",
        type=int,
        default=1,
        help="Parallel workers for final compare/DockQ step",
    )
    parser.add_argument("--rank", type=parse_bool, default=False)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--rank-method", choices=["baseline", "prodigy"], default="baseline")
    parser.add_argument("--prodigy-executable", default="prodigy")
    parser.add_argument("--prodigy-output-dir", default="processed/ranking/prodigy")
    parser.add_argument("--prodigy-distance-cutoff", type=float, default=5.5)
    parser.add_argument("--prodigy-acc-threshold", type=float, default=0.05)
    parser.add_argument("--prodigy-temperature", type=float, default=25.0)
    parser.add_argument("--prodigy-timeout", type=float, default=120.0)
    return parser


if __name__ == "__main__":
    main(build_parser().parse_args())
