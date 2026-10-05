import argparse
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

from src.pdb_download import pdb_downloader
from src.analyse_pdbs import run_analysis
from src.template_generate import template_generator
from src.surface_extract import extract_surfaces
from src.alignment import align
from src.alignment_gtalign import align_gtalign
from src.alignment_multiprot import align_multiprot
from src.transformation import transformer
from src.rosetta_refinement import refiner
from src.pyrosetta_refinement import refine_pairs as pyrosetta_refiner
from src import fiberdock_refinement
from src.compare import compare_pairs_from_outputs
from src.candidate_selector import select_top_candidates
from src.pipeline_inputs import select_templates
from src.transformation_config import TransformationThresholds


def parse_bool(value):
    """Parse CLI booleans without Python's bool('false') trap."""
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y", "on"}:
        return True
    if normalized in {"false", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"invalid boolean value: {value}")


def record_stage_event(stage, event, return_code=None, detail=""):
    """Append an opt-in, machine-readable pipeline stage event."""
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
        record_stage_event(stage, "failed", return_code=1, detail=f"{type(exc).__name__}: {exc}")
        raise
    record_stage_event(stage, "completed", return_code=0)
    return result


_THRESHOLD_ARGUMENTS = (
    ("minimum_residue_match_count", "minimum_residue_match_count"),
    ("minimum_residue_match_percentage", "minimum_residue_match_percentage"),
    ("minimum_hotspot_match_number", "minimum_hotspot_match_number"),
    ("diff_percentage", "diff_percentage"),
    ("template_residue_count", "template_residue_count"),
    ("contact_count_threshold", "contact_count_threshold"),
    ("clashing_distance", "clashing_distance"),
    ("max_clashing_count", "max_clashing_count"),
    ("scaffold_threshold", "scaffold_threshold"),
    ("tm_score_threshold", "tm_score_threshold"),
    ("multiprot_minimum_residue_match_count", "multiprot_minimum_residue_match_count"),
    (
        "multiprot_minimum_residue_match_percentage",
        "multiprot_minimum_residue_match_percentage",
    ),
    ("alignment_gate_mode", "alignment_gate_mode"),
)


def threshold_overrides_from_args(args):
    """Return only explicitly supplied transformation threshold options."""

    return {
        field: value
        for field, attribute in _THRESHOLD_ARGUMENTS
        if (value := getattr(args, attribute, None)) is not None
    }

def main(args):
    os.environ["PRISM_SURFACE_BACKEND"] = args.surface_backend
    if args.freesasa_python:
        os.environ["PRISM_FREESASA_PYTHON"] = args.freesasa_python
    print("PDB download stage started...")
    inputs_csv = getattr(args, "inputs_csv", None)
    download_operation = (
        (lambda: pdb_downloader(inputs_csv)) if inputs_csv else pdb_downloader
    )
    receptor_targets, ligand_targets = run_stage("input", download_operation)
    targets = receptor_targets + ligand_targets
    for r, l in zip(receptor_targets, ligand_targets):
        print(r, "->", l)
    print("PDB download stage finished...")

    threshold_overrides = threshold_overrides_from_args(args)
    configured_thresholds = (
        TransformationThresholds.from_environment().with_overrides(threshold_overrides)
        if threshold_overrides
        else None
    )

    if args.generate_templates:
        print("Filtering templates stage started...")
        results, filtered_count, failed_count = run_analysis()
        print(f"Analysed {len(results)} template files.")
        print(f"Filtered {filtered_count} templates")
        print(f"Failed {failed_count} templates")
        print("Filtering templates stage finished...")

        print("Template generation stage started...")
        templates = template_generator()
        print("Templates generated, templates length", len(templates))
        print("Template generation stage finished...")
    elif not getattr(args, "templates", None) and not getattr(args, "template_list", None):
        with open("templates/calculated_templates.txt", "r") as f:
            templates = [line.strip() for line in f.readlines()]
    else:
        templates = []
    templates = select_templates(
        templates,
        explicit_templates=getattr(args, "templates", None),
        template_list_path=getattr(args, "template_list", None),
        template_limit=getattr(args, "template_limit", None),
    )
    print("Templates selected, templates length", len(templates))

    print("Surface extraction stage started...")
    scaffold_threshold = threshold_overrides.get("scaffold_threshold")
    if scaffold_threshold is None:
        extract_surfaces(targets)
    else:
        extract_surfaces(targets, scaffold_threshold=scaffold_threshold)
    print("Surface extraction stage finished...")

    print("Structural alignment stage started...")
    # Each aligner writes to its own run-scoped directory so tools never share
    # or clobber one another's output (processed/alignment_<tool>/<run_id>/).
    run_id = os.environ.get("PRISM_RUN_ID") or f"{time.strftime('%Y%m%d%H%M%S')}-{os.getpid()}"

    if args.aligner in {"tmalign", "usalign"}:
        alignment_backend_dir = "alignment_usalign" if args.aligner == "usalign" else "alignment_tmalign"
        alignment_output_dir = os.path.join("processed", alignment_backend_dir, run_id)
        prior_alignment_name = os.environ.get("PRISM_ALIGNMENT_NAME")
        prior_tmalign = os.environ.get("PRISM_TMALIGN")
        prior_usalign_fast = os.environ.get("PRISM_USALIGN_FAST")
        try:
            if args.aligner == "usalign":
                os.environ["PRISM_ALIGNMENT_NAME"] = "USalign"
                os.environ["PRISM_TMALIGN"] = args.usalign_path
                os.environ["PRISM_USALIGN_FAST"] = str(args.usalign_fast).lower()
            else:
                os.environ["PRISM_ALIGNMENT_NAME"] = "TMalign"
            run_stage("alignment", lambda: align(targets, templates, output_dir=alignment_output_dir))
        finally:
            if prior_alignment_name is None:
                os.environ.pop("PRISM_ALIGNMENT_NAME", None)
            else:
                os.environ["PRISM_ALIGNMENT_NAME"] = prior_alignment_name
            if prior_tmalign is None:
                os.environ.pop("PRISM_TMALIGN", None)
            else:
                os.environ["PRISM_TMALIGN"] = prior_tmalign
            if prior_usalign_fast is None:
                os.environ.pop("PRISM_USALIGN_FAST", None)
            else:
                os.environ["PRISM_USALIGN_FAST"] = prior_usalign_fast
    elif args.aligner == "multiprot":
        alignment_output_dir = os.path.join("processed", "alignment_multiprot", run_id)
        run_stage("alignment", lambda: align_multiprot(
            targets, templates,
            output_dir=alignment_output_dir,
            max_workers=getattr(
                args,
                "multiprot_workers",
                int(os.environ.get("PRISM_MULTIPROT_WORKERS", "8")),
            ),
            multiprot_path=getattr(args, "multiprot_path", None),
            multiprot_mode=getattr(args, "multiprot_mode", "current"),
            multiprot_params=getattr(args, "multiprot_params", None),
            multiprot_solutions=getattr(args, "multiprot_solutions", 3),
        ))
    else:
        alignment_output_dir = os.path.join("processed", "alignment_gtalign", run_id)
        run_stage(
            "alignment",
            lambda: align_gtalign(
                targets,
                templates,
                gtalign_path=args.gtalign_path,
                output_dir=alignment_output_dir,
                dev_min_length=args.gtalign_dev_min_length,
                pre_score=args.gtalign_pre_score,
                speed=args.gtalign_speed,
                refinement=args.gtalign_refinement,
                **(
                    {"min_match_count": getattr(args, "minimum_residue_match_count")}
                    if getattr(args, "minimum_residue_match_count", None) is not None
                    else {}
                ),
            ),
        )
    print("Structural alignment stage finished...")

    rank_enabled = getattr(args, "rank", False)
    if rank_enabled and args.top_k < 1:
        raise ValueError("--top-k must be positive when --rank is enabled")
    audit_path = getattr(args, "candidate_audit_path", None) or os.environ.get("PRISM_CANDIDATE_AUDIT_PATH")
    if rank_enabled and not audit_path:
        run_id = os.environ.get("PRISM_RUN_ID") or f"{time.strftime('%Y%m%d%H%M%S')}-{os.getpid()}"
        audit_path = os.path.join("processed", "candidate_audit", f"{run_id}.jsonl")

    print("Transformation filtering stage started...")
    orientation = getattr(args, "orientation", "native")

    def transform_operation():
        transformer_kwargs = {
            "alignment_dir": alignment_output_dir,
            "inputs_csv": inputs_csv or "inputs.csv",
        }
        if audit_path:
            transformer_kwargs["audit_path"] = audit_path
        if orientation != "native":
            transformer_kwargs["orientation"] = orientation
        if configured_thresholds is not None:
            transformer_kwargs["thresholds"] = configured_thresholds
        return transformer(templates, **transformer_kwargs)

    passed_pairs = run_stage("transformation", transform_operation)
    print("Passed pairs", len(passed_pairs))
    for pair in passed_pairs:
        print(pair)
    print("Transformation filtering stage finished...")

    # Optional ranking stage: select top-K candidates per pair before refinement
    if rank_enabled:
        print("Candidate ranking stage started...")
        print(f"  Candidate audit: {audit_path}")
        rank_method = getattr(args, "rank_method", "baseline")
        print(f"  Ranking method: {rank_method}")
        passed_pairs = run_stage(
            "ranking",
            lambda: select_top_candidates(
                passed_pairs,
                audit_path=audit_path,
                top_k=args.top_k,
                min_score=args.rank_min_score,
                rank_method=rank_method,
                prodigy_executable=getattr(args, "prodigy_executable", "prodigy"),
                prodigy_output_dir=getattr(args, "prodigy_output_dir", "processed/ranking/prodigy"),
                prodigy_distance_cutoff=getattr(args, "prodigy_distance_cutoff", 5.5),
                prodigy_acc_threshold=getattr(args, "prodigy_acc_threshold", 0.05),
                prodigy_temperature=getattr(args, "prodigy_temperature", 25.0),
                prodigy_timeout=getattr(args, "prodigy_timeout", 120.0),
            ),
        )
        print(f"  Selected {len(passed_pairs)} pairs after ranking (top {args.top_k} per receptor-ligand)")
        print("Candidate ranking stage finished...")

    if getattr(args, "refine", True):
        print(f"{args.refiner} refinement stage started...")
        if not passed_pairs:
            record_stage_event(
                "refinement",
                "skipped",
                return_code=0,
                detail="no candidates passed transformation filtering",
            )
        elif args.refiner == "external_rosetta":
            run_stage("refinement", lambda: refiner(passed_pairs))
        elif args.refiner == "pyrosetta":
            run_stage("refinement", lambda: pyrosetta_refiner(
                passed_pairs,
                output_root=getattr(args, "pyrosetta_output_dir", "processed/pyrosetta_refinement"),
                init_options=getattr(args, "pyrosetta_init_options", None),
            ))
        elif args.refiner == "fiberdock":
            fiberdock_refinement.FIBERDOCK_DIR = os.path.abspath(
                getattr(args, "fiberdock_dir", os.environ.get("PRISM_FIBERDOCK_DIR", "external_tools/fiberdock"))
            )
            run_stage("refinement", lambda: fiberdock_refinement.refine_pairs(passed_pairs))
        else:
            raise ValueError(f"Unknown refiner: {args.refiner}")
        print(f"{args.refiner} refinement stage finished...")
    else:
        record_stage_event(
            "refinement",
            "skipped",
            return_code=0,
            detail="refinement disabled by --no-refine",
        )

    if getattr(args, "compare", False):
        print("Compare (DockQ evaluation) stage started...")
        summary_csv, rows = run_stage(
            "compare",
            lambda: compare_pairs_from_outputs(
                passed_pairs,
                dockq_no_align=args.dockq_no_align,
                n_jobs=args.compare_jobs,
            ),
        )
        print(f"  compared {len(rows)} pairs")
        print(f"  summary written to {summary_csv}")
        print("Compare (DockQ evaluation) stage finished...")

def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--inputs_csv", "--inputs-csv", dest="inputs_csv",
        default=os.environ.get("PRISM_INPUTS_CSV", "inputs.csv"),
        help="CSV containing Receptor and Ligand target IDs.",
    )
    parser.add_argument(
        "--generate_templates",
        "--generate-templates",
        dest="generate_templates",
        nargs="?",
        const=True,
        type=parse_bool,
        default=False,
    )
    # Local change: keep TMalign as default; allow GTalign on user request.
    parser.add_argument("--aligner", choices=["tmalign", "usalign", "gtalign", "multiprot"], default="tmalign")
    parser.add_argument(
        "--orientation",
        choices=["native", "o1", "o2"],
        default="native",
        help=(
            "Transformation orientation mode: native evaluates both implicit "
            "MultiProt assignments; o1/o2 constrain the comparison run."
        ),
    )
    parser.add_argument(
        "--refiner",
        choices=["external_rosetta", "pyrosetta", "fiberdock"],
        default=os.environ.get("PRISM_REFINER", "external_rosetta"),
        help="Final refinement backend; PyRosetta is explicit and never a fallback.",
    )
    refinement = parser.add_mutually_exclusive_group()
    refinement.add_argument("--refine", dest="refine", action="store_true")
    refinement.add_argument("--no-refine", dest="refine", action="store_false")
    parser.set_defaults(refine=True)
    parser.add_argument(
        "--template-limit", "--template_limit",
        dest="template_limit",
        type=int,
        default=None,
        help="Bounded template panel for smoke tests; omit for the full manifest panel.",
    )
    template_source = parser.add_mutually_exclusive_group()
    template_source.add_argument(
        "--template-list", "--template_list",
        dest="template_list",
        default=None,
        help="Plain-text template manifest; one six-character template ID per line.",
    )
    template_source.add_argument(
        "--templates",
        nargs="+",
        dest="templates",
        default=None,
        help="Explicit template IDs (space- or comma-separated); overrides the generated/default panel.",
    )
    parser.add_argument(
        "--minimum-residue-match-count", "--minimum_residue_match_count",
        type=int,
        default=None,
        help="Minimum aligned residue count for TMalign/GTalign and the transformation gate.",
    )
    parser.add_argument(
        "--minimum-residue-match-percentage", "--minimum_residue_match_percentage",
        type=float,
        default=None,
        help="Minimum interface coverage percentage before the large-template adjustment.",
    )
    parser.add_argument(
        "--minimum-hotspot-match-number", "--minimum_hotspot_match_number",
        type=int,
        default=None,
        help="Minimum hotspot matches in published-protocol filtering mode.",
    )
    parser.add_argument(
        "--diff-percentage", "--diff_percentage",
        type=float,
        default=None,
        help="Coverage reduction applied to templates above --template-residue-count.",
    )
    parser.add_argument(
        "--template-residue-count", "--template_residue_count",
        type=float,
        default=None,
        help="Template-size boundary at which --diff-percentage is applied.",
    )
    parser.add_argument(
        "--contact-count-threshold", "--contact_count_threshold",
        "--contact-count", "--contact_count",
        type=int,
        default=None,
        help="Minimum complementary contacts in published-protocol filtering mode.",
    )
    parser.add_argument(
        "--clashing-distance", "--clashing_distance",
        type=float,
        default=None,
        help="C-alpha distance in Angstroms below which a cross-partner clash is counted.",
    )
    parser.add_argument(
        "--max-clashing-count", "--max_clashing_count",
        type=int,
        default=None,
        help="Reject a transformed pair when the clash count reaches this value.",
    )
    parser.add_argument(
        "--scaffold-threshold", "--scaffold_threshold",
        type=float,
        default=None,
        help="C-alpha radius in Angstroms used to build the surface scaffold (default: 5.0).",
    )
    parser.add_argument(
        "--tm-score-threshold", "--tm_score_threshold",
        type=float,
        default=None,
        help="Minimum TM-score for TMalign/GTalign transformation filtering.",
    )
    parser.add_argument(
        "--alignment-gate-mode", "--alignment_gate_mode",
        choices=["native", "common_match_coverage"],
        default=None,
        help=(
            "Alignment acceptance contract: native preserves provider-specific "
            "gates; common_match_coverage uses shared matched-residue and "
            "interface-coverage gates without a TM-score gate."
        ),
    )
    parser.add_argument(
        "--multiprot-minimum-residue-match-count", "--multiprot_minimum_residue_match_count",
        type=int,
        default=None,
        help="Minimum native MultiProt correspondence count.",
    )
    parser.add_argument(
        "--multiprot-minimum-residue-match-percentage", "--multiprot_minimum_residue_match_percentage",
        type=float,
        default=None,
        help="Minimum native MultiProt interface coverage percentage.",
    )
    parser.add_argument(
        "--surface_backend", "--surface-backend",
        dest="surface_backend",
        choices=["naccess", "freesasa"],
        default=os.environ.get("PRISM_SURFACE_BACKEND", "naccess"),
        help="Surface-area backend; NACCESS remains the default.",
    )
    parser.add_argument(
        "--freesasa_python", "--freesasa-python",
        dest="freesasa_python",
        type=str,
        default=os.environ.get("PRISM_FREESASA_PYTHON"),
        help="Python interpreter containing FreeSASA when --surface_backend=freesasa.",
    )
    parser.add_argument("--gtalign_path", "--gtalign-path", dest="gtalign_path", default="gtalign")
    parser.add_argument(
        "--usalign-path", dest="usalign_path",
        default=os.environ.get("PRISM_USALIGN", "/home/rshadi25/.conda/envs/gtalign_env/bin/USalign"),
        help="USalign executable used when --aligner=usalign.",
    )
    parser.add_argument(
        "--usalign-fast", dest="usalign_fast", action="store_true",
        default=os.environ.get("PRISM_USALIGN_FAST", "").strip().lower() in {"1", "true", "yes", "on"},
        help="Pass USalign -fast for alignment-only pilot/production arms.",
    )
    parser.add_argument("--gtalign_dev_min_length", "--gtalign-dev-min-length", dest="gtalign_dev_min_length", type=int, default=3)
    parser.add_argument("--gtalign_pre_score", "--gtalign-pre-score", dest="gtalign_pre_score", type=float, default=0.0)
    parser.add_argument("--gtalign_speed", "--gtalign-speed", dest="gtalign_speed", type=int, default=0)
    parser.add_argument("--gtalign_refinement", "--gtalign-refinement", dest="gtalign_refinement", type=int, default=3)
    parser.add_argument("--multiprot-workers", type=int, default=int(os.environ.get("PRISM_MULTIPROT_WORKERS", "8")))
    parser.add_argument("--multiprot-path", default=os.environ.get("PRISM_MULTIPROT", "external_tools/multiprot.Linux"))
    parser.add_argument(
        "--multiprot-mode",
        choices=["current", "legacy_compatible"],
        default=os.environ.get("PRISM_MULTIPROT_MODE", "current"),
        help="MultiProt invocation/parser mode; legacy_compatible is opt-in.",
    )
    parser.add_argument(
        "--multiprot-params",
        default=os.environ.get("PRISM_MULTIPROT_PARAMS"),
        help="Optional params.txt copied into each isolated MultiProt run.",
    )
    parser.add_argument(
        "--multiprot-solutions",
        type=int,
        default=int(os.environ.get("PRISM_MULTIPROT_SOLUTIONS", "3")),
        help="Number of legacy MultiProt solutions retained per alignment.",
    )
    parser.add_argument("--pyrosetta-output-dir", default="processed/pyrosetta_refinement")
    parser.add_argument("--pyrosetta-init-options", default="-mute all -constant_seed -jran 12345")
    parser.add_argument("--fiberdock-dir", default=os.environ.get("PRISM_FIBERDOCK_DIR", "external_tools/fiberdock"))
    parser.add_argument(
        "--compare",
        nargs="?",
        const=True,
        type=parse_bool,
        default=os.environ.get("PRISM_COMPARE", "").lower() in ("true", "1", "yes"),
        help="Run DockQ evaluation after refinement (opt-in).",
    )
    parser.add_argument(
        "--dockq-no-align",
        nargs="?",
        const=True,
        type=parse_bool,
        default=False,
        help="Skip superposition before DockQ calculation.",
    )
    parser.add_argument(
        "--compare-jobs",
        type=int,
        default=int(os.environ.get("PRISM_COMPARE_JOBS", "1")),
        help="Parallel workers for DockQ comparison.",
    )
    parser.add_argument(
        "--rank",
        type=parse_bool,
        default=os.environ.get("PRISM_RANK", "").lower() in ("true", "1", "yes"),
        help="Enable candidate ranking before refinement (opt-in).",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=int(os.environ.get("PRISM_TOP_K", "5")),
        help="Maximum candidates to keep per receptor-ligand pair after ranking.",
    )
    parser.add_argument(
        "--rank-min-score",
        type=float,
        default=float(os.environ.get("PRISM_RANK_MIN_SCORE", "0.0")),
        help="Minimum baseline score threshold for ranking (0-1).",
    )
    parser.add_argument(
        "--rank-method",
        choices=["baseline", "prodigy"],
        default=os.environ.get("PRISM_RANK_METHOD", "baseline"),
        help="Candidate ranking scorer. PRODIGY is opt-in and requires a separate installed executable.",
    )
    parser.add_argument(
        "--prodigy-executable",
        type=str,
        default=os.environ.get("PRISM_PRODIGY_EXECUTABLE", "prodigy"),
        help="PRODIGY command or executable path used with --rank-method prodigy.",
    )
    parser.add_argument(
        "--prodigy-output-dir",
        type=str,
        default=os.environ.get("PRISM_PRODIGY_OUTPUT_DIR", "processed/ranking/prodigy"),
        help="Directory for combined PRODIGY inputs, stdout/stderr, and score metadata.",
    )
    parser.add_argument(
        "--prodigy-distance-cutoff",
        type=float,
        default=float(os.environ.get("PRISM_PRODIGY_DISTANCE_CUTOFF", "5.5")),
        help="PRODIGY intermolecular contact distance cutoff in Angstroms.",
    )
    parser.add_argument(
        "--prodigy-acc-threshold",
        type=float,
        default=float(os.environ.get("PRISM_PRODIGY_ACC_THRESHOLD", "0.05")),
        help="PRODIGY accessibility threshold for interface analysis.",
    )
    parser.add_argument(
        "--prodigy-temperature",
        type=float,
        default=float(os.environ.get("PRISM_PRODIGY_TEMPERATURE", "25.0")),
        help="PRODIGY temperature in Celsius for the reported affinity model.",
    )
    parser.add_argument(
        "--prodigy-timeout",
        type=float,
        default=float(os.environ.get("PRISM_PRODIGY_TIMEOUT", "120")),
        help="Maximum seconds allowed for one PRODIGY candidate score.",
    )
    parser.add_argument(
        "--candidate-audit-path",
        type=str,
        default=os.environ.get("PRISM_CANDIDATE_AUDIT_PATH"),
        help="JSONL candidate audit path; ranked runs otherwise receive a fresh run-scoped path.",
    )
    return parser


if __name__ == "__main__":
    main(build_parser().parse_args())
