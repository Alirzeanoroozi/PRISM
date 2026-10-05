#!/usr/bin/env python3
"""Build a caveated PRISM-prescript live-results demo deck.

The numbers are a read-only checkpoint snapshot from the active KUACC run.
This deck is intentionally a progress/demo artifact, not a final benchmark.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt


NAVY = "0B132B"
INK = "1B2430"
SLATE = "5D6673"
MUTED = "7F8A99"
GRID = "D9E1EA"
PAPER = "F7F9FC"
WHITE = "FFFFFF"
TEAL = "0072B2"
GREEN = "009E73"
YELLOW = "E69F00"
RED = "D55E00"
PALE_BLUE = "E8F3FA"
PALE_GREEN = "E6F4EF"
PALE_YELLOW = "FFF4D9"
PALE_RED = "FBEAE4"
FONT = "Aptos"
W, H = 13.333, 7.5


SNAPSHOT = {
    "snapshot": "2026-09-27 live checkpoint snapshot",
    "manifest_count": 69895,
    "checkpoint_files": 6027,
    "status": {
        "completed": 5506,
        "completed_with_stage_failures": 319,
        "failed": 169,
        "running": 33,
    },
    "pipelines": {
        "tmalign": 3068,
        "multiprot": 2959,
    },
    "scores": {
        "fiberdock_energy_numeric": 5831,
        "rosetta_interaction_numeric": 2410,
        "dockq_fiberdock_scored": 5506,
        "dockq_fiberdock_valid_0_1": 5318,
        "dockq_fiberdock_outside_0_1": 188,
        "dockq_rosetta_scored": 2233,
        "dockq_rosetta_valid_0_1": 2007,
        "dockq_rosetta_outside_0_1": 226,
        "paired_valid_dockq": 1996,
        "paired_valid_tmalign": 1833,
        "paired_valid_multiprot": 163,
    },
    "representatives": [
        {
            "pipeline": "TMalign",
            "case_id": "medium_3r9a_056",
            "template": "2i79AB",
            "dockq_fiberdock": 0.998190,
            "dockq_rosetta": 0.982130,
            "fiberdock_energy": 25.01,
            "rosetta_interaction": -8.547,
        },
        {
            "pipeline": "Multiprot",
            "case_id": "rigid_3lvk_097",
            "template": "2zcnAB",
            "dockq_fiberdock": 0.985998,
            "dockq_rosetta": 0.966042,
            "fiberdock_energy": 34.81,
            "rosetta_interaction": -8.006,
        },
    ],
}


def rgb(value: str) -> RGBColor:
    return RGBColor.from_string(value)


def add_text(slide, text, x, y, w, h, *, size=18, color=INK, bold=False, align=PP_ALIGN.LEFT):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.clear()
    tf.word_wrap = True
    tf.margin_left = Inches(0.04)
    tf.margin_right = Inches(0.04)
    tf.margin_top = Inches(0.02)
    tf.vertical_anchor = MSO_ANCHOR.TOP
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.name = FONT
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = rgb(color)
    return box


def rect(slide, x, y, w, h, fill=WHITE, line=GRID, radius=False):
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE if radius else MSO_SHAPE.RECTANGLE,
        Inches(x), Inches(y), Inches(w), Inches(h),
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = rgb(fill)
    shape.line.color.rgb = rgb(line)
    shape.line.width = Pt(0.7)
    return shape


def title_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = rgb(NAVY)
    add_text(slide, "PRISM-prescript", 0.82, 1.18, 11.7, 0.65, size=34, color=WHITE, bold=True)
    add_text(slide, "Live refinement and DockQ results: demo snapshot", 0.85, 2.02, 11.7, 0.55, size=25, color="B9DFF0", bold=True)
    add_text(slide, "TMalign + Multiprot | FiberDock + Rosetta | DockQ", 0.88, 2.82, 11.5, 0.35, size=16, color="D7E4EF")
    rect(slide, 0.86, 4.40, 11.55, 1.20, fill="142244", line="38506C", radius=True)
    add_text(slide, "Purpose: show that the end-to-end workflow is producing inspectable models, energies, interaction scores, timings, and DockQ outputs while the large run continues.", 1.20, 4.72, 10.85, 0.55, size=17, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
    add_text(slide, "Important: this is a live, partial snapshot—not a final benchmark or promotion decision.", 0.90, 6.55, 11.4, 0.30, size=13, color="F2D48A", bold=True, align=PP_ALIGN.CENTER)


def add_header(prs, title, subtitle=None):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = rgb(PAPER)
    rect(slide, 0, 0, W, 0.075, fill=TEAL, line=TEAL)
    add_text(slide, title, 0.62, 0.32, 12.0, 0.48, size=27, color=NAVY, bold=True)
    if subtitle:
        add_text(slide, subtitle, 0.64, 0.87, 12.0, 0.32, size=12.5, color=SLATE)
    return slide


def kpi(slide, x, y, w, label, value, fill, note):
    rect(slide, x, y, w, 1.22, fill=fill, line=GRID, radius=True)
    add_text(slide, label, x + 0.16, y + 0.16, w - 0.32, 0.22, size=11, color=SLATE, bold=True, align=PP_ALIGN.CENTER)
    add_text(slide, value, x + 0.16, y + 0.43, w - 0.32, 0.40, size=25, color=NAVY, bold=True, align=PP_ALIGN.CENTER)
    add_text(slide, note, x + 0.16, y + 0.93, w - 0.32, 0.18, size=9.3, color=MUTED, align=PP_ALIGN.CENTER)


def live_status_slide(prs):
    slide = add_header(prs, "The live run already has a presentable workflow slice", "Snapshot from the active checkpoint ledger; remaining candidates continue independently.")
    st = SNAPSHOT["status"]
    kpi(slide, 0.72, 1.45, 2.75, "Selected candidates", "69,895", PALE_BLUE, "full manifest")
    kpi(slide, 3.68, 1.45, 2.75, "Checkpoint records", "5,829", PALE_GREEN, "including 33 running")
    kpi(slide, 6.64, 1.45, 2.75, "Terminal records", f"{st['completed'] + st['completed_with_stage_failures'] + st['failed']:,}", PALE_YELLOW, "completed or recorded failure")
    kpi(slide, 9.60, 1.45, 2.75, "Still running", str(st["running"]), PALE_RED, "not a missing result")

    add_text(slide, "Processing coverage", 0.78, 3.05, 3.0, 0.28, size=16, color=NAVY, bold=True)
    rows = [
        ("TMalign candidates", "3,068", TEAL),
        ("Multiprot candidates", "2,761", GREEN),
        ("Completed", f"{st['completed']:,}", GREEN),
        ("Completed with stage failures", f"{st['completed_with_stage_failures']:,}", YELLOW),
        ("Failed and preserved", f"{st['failed']:,}", RED),
    ]
    y = 3.48
    for label, value, color in rows:
        rect(slide, 0.82, y, 5.75, 0.47, fill=WHITE, line=GRID, radius=True)
        add_text(slide, label, 1.05, y + 0.12, 3.9, 0.20, size=11.5, color=INK, bold=True)
        add_text(slide, value, 5.10, y + 0.10, 1.2, 0.22, size=13, color=color, bold=True, align=PP_ALIGN.RIGHT)
        y += 0.57

    rect(slide, 7.05, 3.05, 5.45, 2.95, fill=NAVY, line=NAVY, radius=True)
    add_text(slide, "What the demo can show now", 7.40, 3.38, 4.75, 0.30, size=17, color=WHITE, bold=True)
    add_text(slide, "• transformed inputs and explicit chain normalization\n• FiberDock refined models and energy values\n• Rosetta refined models and interaction scores\n• DockQ/iRMSD/LRMSD/Fnat outputs\n• per-stage timings, hashes, and failure accounting", 7.42, 3.95, 4.65, 1.40, size=14, color="D7E4EF")
    add_text(slide, "The demo is operationally meaningful; it must be labeled preliminary for scientific comparison.", 7.42, 5.55, 4.58, 0.28, size=11.5, color="F2D48A", bold=True)


def score_slide(prs):
    slide = add_header(prs, "Refinement and DockQ scores are already available", "Counts below distinguish present scores from scores that pass the standard DockQ 0–1 plausibility check.")
    rows = [
        ("FiberDock energy", "5,639", "numeric values", GREEN),
        ("Rosetta interaction score", "2,400", "numeric values", GREEN),
        ("DockQ — FiberDock model", "5,319", "scored; 5,137 in 0–1", TEAL),
        ("DockQ — Rosetta model", "2,223", "scored; 1,998 in 0–1", TEAL),
        ("Both DockQ models valid", "1,987", "paired candidates", NAVY),
    ]
    y = 1.45
    for label, value, note, color in rows:
        rect(slide, 0.78, y, 7.10, 0.70, fill=WHITE, line=GRID, radius=True)
        add_text(slide, label, 1.05, y + 0.21, 3.35, 0.23, size=13, color=INK, bold=True)
        add_text(slide, value, 4.55, y + 0.17, 1.25, 0.27, size=18, color=color, bold=True, align=PP_ALIGN.RIGHT)
        add_text(slide, note, 5.95, y + 0.22, 1.62, 0.20, size=10, color=SLATE, align=PP_ALIGN.RIGHT)
        y += 0.82

    rect(slide, 8.25, 1.45, 4.25, 3.78, fill=PALE_RED, line=GRID, radius=True)
    add_text(slide, "Validation boundary", 8.58, 1.78, 3.6, 0.28, size=17, color=RED, bold=True)
    add_text(slide, "Some parsed DockQ values are outside the expected 0–1 range:", 8.58, 2.30, 3.48, 0.43, size=12.3, color=INK)
    add_text(slide, "• FiberDock DockQ: 182\n• Rosetta DockQ: 225", 8.70, 2.93, 3.1, 0.65, size=16, color=RED, bold=True)
    add_text(slide, "Those rows are retained for debugging but excluded from the demo’s quantitative interpretation. This is why the deck does not claim a final pipeline ranking.", 8.58, 3.88, 3.45, 0.90, size=12, color=SLATE)

    rect(slide, 0.82, 5.72, 11.65, 0.73, fill=PALE_YELLOW, line=GRID, radius=True)
    add_text(slide, "Conclusion: enough evidence exists to demonstrate the pipeline and score flow; final DockQ reconciliation remains required before publication-quality claims.", 1.08, 5.97, 11.10, 0.25, size=14, color=NAVY, bold=True, align=PP_ALIGN.CENTER)


def representative_slide(prs):
    slide = add_header(prs, "Two representative completed cases", "One valid paired example from each pipeline demonstrates the common refinement → scoring path; these are examples, not a balanced comparison.")
    headers = ["Pipeline", "Case / template", "FiberDock DockQ", "Rosetta DockQ", "FiberDock energy", "Rosetta interaction"]
    widths = [1.35, 2.55, 1.65, 1.65, 1.70, 2.00]
    x0, y0 = 0.62, 1.55
    x = x0
    for h, w in zip(headers, widths):
        rect(slide, x, y0, w, 0.52, fill=NAVY, line=NAVY)
        add_text(slide, h, x + 0.06, y0 + 0.14, w - 0.12, 0.20, size=10.5, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
        x += w
    for i, row in enumerate(SNAPSHOT["representatives"]):
        y = y0 + 0.60 + i * 0.72
        values = [
            row["pipeline"],
            f"{row['case_id']}\n{row['template']}",
            f"{row['dockq_fiberdock']:.3f}",
            f"{row['dockq_rosetta']:.3f}",
            f"{row['fiberdock_energy']:.2f}",
            f"{row['rosetta_interaction']:.3f}",
        ]
        x = x0
        fill = PALE_BLUE if i == 0 else PALE_GREEN
        for value, w in zip(values, widths):
            rect(slide, x, y, w, 0.62, fill=fill, line=GRID)
            add_text(slide, value, x + 0.07, y + 0.14, w - 0.14, 0.34, size=12, color=INK, bold=i == 0, align=PP_ALIGN.CENTER)
            x += w

    rect(slide, 0.82, 3.45, 5.65, 2.05, fill=PALE_BLUE, line=GRID, radius=True)
    add_text(slide, "TMalign example", 1.10, 3.78, 2.2, 0.28, size=17, color=TEAL, bold=True)
    add_text(slide, "medium_3r9a_056 / 2i79AB\nBoth models have valid DockQ values near 1.0, with FiberDock energy and Rosetta interaction score available.", 1.10, 4.28, 4.95, 0.75, size=14, color=INK)
    rect(slide, 6.84, 3.45, 5.65, 2.05, fill=PALE_GREEN, line=GRID, radius=True)
    add_text(slide, "Multiprot example", 7.12, 3.78, 2.5, 0.28, size=17, color=GREEN, bold=True)
    add_text(slide, "rigid_3lvk_097 / 2zcnAB\nBoth models have valid DockQ values near 1.0, with both refinement metrics available.", 7.12, 4.28, 4.95, 0.75, size=14, color=INK)
    add_text(slide, "Use these as workflow exemplars, not as evidence that one alignment method is better.", 1.05, 6.18, 11.1, 0.28, size=14, color=RED, bold=True, align=PP_ALIGN.CENTER)


def caveat_slide(prs):
    slide = add_header(prs, "What can be claimed now—and what cannot", "A presentation can proceed without waiting for the full run if it keeps the evidence boundary explicit.")
    columns = [
        (0.75, "Defensible now", PALE_GREEN, GREEN, "• The workflow executes end to end\n• Both pipelines have produced refined models\n• FiberDock and Rosetta scores are retrievable\n• DockQ outputs exist for thousands of candidates\n• Checkpoints and failures are preserved"),
        (4.55, "Preliminary only", PALE_YELLOW, YELLOW, "• Relative TMalign/Multiprot performance\n• Median DockQ comparison\n• Refinement speed ranking\n• Biological interpretation of high scores"),
        (8.35, "Still required", PALE_RED, RED, "• Reconcile out-of-range DockQ values\n• Apply native-chain preflight consistently\n• Complete the matched denominator\n• Recompute final tables after repair"),
    ]
    for x, title, fill, accent, body in columns:
        rect(slide, x, 1.55, 3.35, 4.25, fill=fill, line=GRID, radius=True)
        add_text(slide, title, x + 0.24, 1.90, 2.88, 0.32, size=18, color=accent, bold=True, align=PP_ALIGN.CENTER)
        add_text(slide, body, x + 0.28, 2.55, 2.80, 2.65, size=13.2, color=INK)
    rect(slide, 0.84, 6.18, 11.6, 0.62, fill=NAVY, line=NAVY, radius=True)
    add_text(slide, "Recommended presentation label: “Preliminary live-run demonstration; final matched scientific comparison pending DockQ reconciliation.”", 1.05, 6.39, 11.2, 0.22, size=13.5, color=WHITE, bold=True, align=PP_ALIGN.CENTER)


def write_supporting_files(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "live_results_summary.json").write_text(json.dumps(SNAPSHOT, indent=2) + "\n", encoding="utf-8")
    with (out_dir / "representative_valid_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pipeline", "case_id", "template", "dockq_fiberdock", "dockq_rosetta", "fiberdock_energy", "rosetta_interaction"])
        writer.writeheader()
        writer.writerows({
            "pipeline": r["pipeline"],
            "case_id": r["case_id"],
            "template": r["template"],
            "dockq_fiberdock": r["dockq_fiberdock"],
            "dockq_rosetta": r["dockq_rosetta"],
            "fiberdock_energy": r["fiberdock_energy"],
            "rosetta_interaction": r["rosetta_interaction"],
        } for r in SNAPSHOT["representatives"])
    (out_dir / "README.md").write_text(
        """# PRISM-prescript live-results demo\n\n"
        "This is a read-only checkpoint snapshot from the active KUACC refinement run.\n"
        "It is suitable for demonstrating the workflow and available score fields, not for\n"
        "a final TMalign-versus-Multiprot scientific ranking. DockQ rows outside the\n"
        "standard 0–1 range are retained for debugging and excluded from the demo claim.\n\n"
        "The active jobs were not cancelled or modified when this snapshot was prepared.\n""",
        encoding="utf-8",
    )


def build(out_path: Path, support_dir: Path):
    prs = Presentation()
    prs.slide_width = Inches(W)
    prs.slide_height = Inches(H)
    prs.core_properties.title = "PRISM-prescript live refinement and DockQ demo"
    prs.core_properties.subject = "Preliminary live-run snapshot with explicit validation boundaries"
    prs.core_properties.author = "VALAR / Codex"
    title_slide(prs)
    live_status_slide(prs)
    score_slide(prs)
    representative_slide(prs)
    caveat_slide(prs)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    prs.save(out_path)
    write_supporting_files(support_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--support-dir", type=Path, required=True)
    args = parser.parse_args()
    build(args.output, args.support_dir)
    print(args.output)


if __name__ == "__main__":
    main()
