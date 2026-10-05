"""
Probe whether a structural aligner can replace TMalign in PRISM.

This script checks for the outputs PRISM currently depends on:
1) TM-score
2) Aligned length
3) A rotation/translation matrix (for coordinate transforms)
4) Residue-level correspondence (for match_dict reconstruction)

It supports a built-in TMalign mode and a generic mode for GTalign-like tools.
"""
import argparse
import re
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path


def write_test_pdb(path: Path, shift=(0.0, 0.0, 0.0)) -> None:
    coords = [
        (1.000, 2.000, 3.000),
        (2.500, 2.100, 3.300),
        (4.000, 2.200, 3.600),
        (5.500, 2.300, 3.900),
    ]
    with open(path, "w") as f:
        for i, (x, y, z) in enumerate(coords, start=1):
            sx, sy, sz = shift
            f.write(
                f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                f"{x+sx:8.3f}{y+sy:8.3f}{z+sz:8.3f}"
                f"  1.00 20.00           C\n"
            )
        f.write("TER\nEND\n")


def run_cmd(cmd, cwd=None, timeout=60):
    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=timeout,
        )
        return {
            "ok": result.returncode == 0,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output": result.stdout + result.stderr,
        }
    except subprocess.TimeoutExpired:
        return {"ok": False, "returncode": None, "stdout": "", "stderr": "Timeout", "output": ""}
    except Exception as exc:
        return {"ok": False, "returncode": None, "stdout": "", "stderr": str(exc), "output": ""}


def parse_tmalign_signals(output: str):
    signals = {
        "has_tm_score": False,
        "has_aligned_length": False,
        "has_rmsd": False,
        "has_residue_mapping_text": False,
    }
    for line in output.splitlines():
        if line.startswith("Aligned length"):
            signals["has_aligned_length"] = True
            if "RMSD=" in line:
                signals["has_rmsd"] = True
        if "TM-score=" in line and "(if normalized" in line:
            signals["has_tm_score"] = True
        if line.startswith('(":"'):
            signals["has_residue_mapping_text"] = True
    return signals


def parse_generic_signals(output: str):
    lowered = output.lower()
    signals = {
        "has_tm_score": bool(re.search(r"\btm[- ]?score\b", lowered)),
        "has_aligned_length": bool(re.search(r"aligned\s+length|ali(?:gned)?\s*len", lowered)),
        "has_rmsd": "rmsd" in lowered,
        "has_residue_mapping_text": False,
    }

    # Heuristics for textual alignment blocks / correspondences
    if '(":"' in output:
        signals["has_residue_mapping_text"] = True
    elif re.search(r"residue\s+pair|correspondence|alignment\s+string", lowered):
        signals["has_residue_mapping_text"] = True

    return signals


def has_tmalign_matrix_format(matrix_path: Path):
    if not matrix_path.exists():
        return False
    rows = set()
    with open(matrix_path, "r") as f:
        for line in f:
            toks = line.strip().split()
            if len(toks) < 5:
                continue
            try:
                row = int(toks[0])
                _ = float(toks[1])
                _ = float(toks[2])
                _ = float(toks[3])
                _ = float(toks[4])
            except ValueError:
                continue
            if row in (0, 1, 2):
                rows.add(row)
    return rows == {0, 1, 2}


def build_builtin_tmalign_cmd(binary: Path, pdb1: Path, pdb2: Path, matrix_out: Path):
    return [str(binary), str(pdb1), str(pdb2), "-m", str(matrix_out)]


def build_template_cmd(template: str, pdb1: Path, pdb2: Path, matrix_out: Path):
    rendered = template.format(pdb1=str(pdb1), pdb2=str(pdb2), matrix=str(matrix_out))
    return shlex.split(rendered)


def main():
    parser = argparse.ArgumentParser(
        description="Probe structural aligner compatibility with PRISM's TMalign-dependent pipeline"
    )
    parser.add_argument("--tool", choices=["tmalign", "gtalign", "custom"], required=True)
    parser.add_argument("--binary", required=True, help="Path to aligner executable")
    parser.add_argument(
        "--cmd-template",
        default=None,
        help=(
            "Custom command template for gtalign/custom mode. "
            "Placeholders: {pdb1}, {pdb2}, {matrix}. Example: "
            "'/path/GTalign {pdb1} {pdb2}'"
        ),
    )
    parser.add_argument(
        "--expect-matrix-file",
        action="store_true",
        help="Require/check an explicit matrix file at {matrix} path",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="Timeout in seconds for aligner execution",
    )
    args = parser.parse_args()

    binary = Path(args.binary).resolve()
    if not binary.exists():
        print(f"ERROR: binary not found: {binary}")
        sys.exit(2)

    with tempfile.TemporaryDirectory(prefix="prism_aligner_probe_") as td:
        work = Path(td)
        pdb1 = work / "probe1.pdb"
        pdb2 = work / "probe2.pdb"
        matrix_out = work / "matrix.out"

        # A shifted copy preserves residue order and should align trivially.
        write_test_pdb(pdb1, shift=(0.0, 0.0, 0.0))
        write_test_pdb(pdb2, shift=(10.0, -2.0, 5.0))

        if args.tool == "tmalign":
            cmd = build_builtin_tmalign_cmd(binary, pdb1, pdb2, matrix_out)
        else:
            if not args.cmd_template:
                print("ERROR: --cmd-template is required for gtalign/custom mode")
                sys.exit(2)
            cmd = build_template_cmd(args.cmd_template, pdb1, pdb2, matrix_out)

        result = run_cmd(cmd, cwd=work, timeout=args.timeout)

        print("== Command ==")
        print(" ".join(shlex.quote(x) for x in cmd))
        print(f"Return code: {result['returncode']}")

        if not result["ok"]:
            print("== STDERR ==")
            print(result["stderr"][:4000])
            print("\nConclusion: NOT POSSIBLE to assess yet (tool did not run successfully).")
            sys.exit(1)

        output = result["output"]
        if args.tool == "tmalign":
            signals = parse_tmalign_signals(output)
            matrix_ok = has_tmalign_matrix_format(matrix_out)
        else:
            signals = parse_generic_signals(output)
            matrix_ok = has_tmalign_matrix_format(matrix_out) if args.expect_matrix_file else False

        print("== Capability Signals ==")
        for key, val in signals.items():
            print(f"{key}: {val}")
        print(f"has_transform_matrix_file: {matrix_ok}")

        # Full PRISM replacement currently needs all of these.
        drop_in_possible = (
            signals["has_tm_score"]
            and signals["has_aligned_length"]
            and signals["has_residue_mapping_text"]
            and matrix_ok
        )
        partial_possible = signals["has_tm_score"] and signals["has_aligned_length"]

        print("== Assessment ==")
        if drop_in_possible:
            print("FULL_REPLACEMENT_POSSIBLE: yes (based on this probe)")
        elif partial_possible:
            print("FULL_REPLACEMENT_POSSIBLE: no")
            print("PARTIAL_USE_POSSIBLE: yes (scoring/prefilter use likely)")
        else:
            print("FULL_REPLACEMENT_POSSIBLE: no")
            print("PARTIAL_USE_POSSIBLE: unclear (missing core scoring signals)")

        # Print a short preview to help inspect unknown tools.
        preview = "\n".join(output.splitlines()[:40])
        print("== Output Preview (first 40 lines) ==")
        print(preview)


if __name__ == "__main__":
    main()
