"""
Benchmark alignment execution patterns relevant to PRISM:

1) Current PRISM style: one TMalign subprocess per query-reference pair
2) GTalign batch style: one GTalign run over many queries and references

This does not modify the pipeline; it measures the practical runtime impact
of replacing the alignment stage invocation strategy.
"""
import argparse
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def write_chain_pdb(path: Path, n_res: int, shift_x=0.0, shift_y=0.0, shift_z=0.0):
    with open(path, "w") as f:
        for i in range(1, n_res + 1):
            x = 1.0 + 1.5 * (i - 1) + shift_x
            y = 2.0 + 0.1 * (i - 1) + shift_y
            z = 3.0 + 0.3 * (i - 1) + shift_z
            f.write(
                f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
            )
        f.write("TER\nEND\n")


def run(cmd, cwd=None, timeout=None):
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=timeout,
    )


def build_dataset(root: Path, n_queries: int, n_refs: int, n_res: int):
    qdir = root / "queries"
    rdir = root / "refs"
    qdir.mkdir(parents=True, exist_ok=True)
    rdir.mkdir(parents=True, exist_ok=True)

    queries = []
    refs = []
    for i in range(n_queries):
        p = qdir / f"q_{i:04d}.pdb"
        write_chain_pdb(p, n_res, shift_x=0.0, shift_y=float(i) * 0.01, shift_z=0.0)
        queries.append(p)
    for j in range(n_refs):
        p = rdir / f"r_{j:04d}.pdb"
        # Translate references so alignments are non-trivial but still exact.
        write_chain_pdb(
            p,
            n_res,
            shift_x=10.0 + float(j) * 0.02,
            shift_y=-2.0,
            shift_z=5.0,
        )
        refs.append(p)
    return queries, refs, qdir, rdir


def bench_tmalign_pairs(tmalign_bin: Path, queries, refs, work: Path, limit_pairs=None):
    work.mkdir(parents=True, exist_ok=True)
    pair_count = 0
    t0 = time.perf_counter()
    for qi, q in enumerate(queries):
        for ri, r in enumerate(refs):
            if limit_pairs is not None and pair_count >= limit_pairs:
                dt = time.perf_counter() - t0
                return {"seconds": dt, "pairs": pair_count, "ok": True, "errors": []}

            matrix_file = work / f"tm_{qi}_{ri}.matrix"
            result = run([str(tmalign_bin), str(q), str(r), "-m", str(matrix_file)], cwd=work, timeout=120)
            pair_count += 1
            if result.returncode != 0:
                dt = time.perf_counter() - t0
                return {
                    "seconds": dt,
                    "pairs": pair_count,
                    "ok": False,
                    "errors": [result.stderr[:1000], result.stdout[:1000]],
                }
    dt = time.perf_counter() - t0
    return {"seconds": dt, "pairs": pair_count, "ok": True, "errors": []}


def bench_gtalign_batch(gtalign_bin: Path, qdir: Path, rdir: Path, outdir: Path, n_refs: int, n_res: int):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(gtalign_bin),
        f"--qrs={qdir}",
        f"--rfs={rdir}",
        "-o",
        str(outdir),
        "-s",
        "0",
        f"--nhits={n_refs}",
        f"--nalns={n_refs}",
        "--dev-min-length=3" if n_res < 20 else f"--dev-min-length={min(20, n_res)}",
    ]
    t0 = time.perf_counter()
    result = run(cmd, timeout=600)
    dt = time.perf_counter() - t0
    out_files = sorted(outdir.glob("*.out"))
    return {
        "seconds": dt,
        "ok": result.returncode == 0,
        "returncode": result.returncode,
        "stdout_preview": result.stdout[:1500],
        "stderr_preview": result.stderr[:1500],
        "output_files": len(out_files),
        "cmd": cmd,
    }


def main():
    parser = argparse.ArgumentParser(description="Compare TMalign per-pair vs GTalign batch runtimes")
    parser.add_argument("--tmalign", required=True, help="Path to TMalign binary")
    parser.add_argument("--gtalign", required=True, help="Path to GTalign binary")
    parser.add_argument("--queries", type=int, default=10, help="Number of query structures")
    parser.add_argument("--refs", type=int, default=10, help="Number of reference structures")
    parser.add_argument("--residues", type=int, default=50, help="Residues per synthetic structure")
    parser.add_argument(
        "--tmalign-limit-pairs",
        type=int,
        default=None,
        help="Optional cap on TMalign pairs to keep runtime short while testing",
    )
    args = parser.parse_args()

    tmalign_bin = Path(args.tmalign).resolve()
    gtalign_bin = Path(args.gtalign).resolve()
    if not tmalign_bin.exists():
        print(f"ERROR: TMalign not found: {tmalign_bin}")
        sys.exit(2)
    if not gtalign_bin.exists():
        print(f"ERROR: GTalign not found: {gtalign_bin}")
        sys.exit(2)

    total_pairs = args.queries * args.refs
    print(f"Dataset: {args.queries} queries x {args.refs} refs = {total_pairs} pairs, {args.residues} residues/structure")

    with tempfile.TemporaryDirectory(prefix="prism_align_bench_") as td:
        root = Path(td)
        queries, refs, qdir, rdir = build_dataset(root, args.queries, args.refs, args.residues)

        print("Running TMalign per-pair benchmark...")
        tm = bench_tmalign_pairs(
            tmalign_bin,
            queries,
            refs,
            root / "tmalign_work",
            limit_pairs=args.tmalign_limit_pairs,
        )
        print(f"TMalign ok={tm['ok']} pairs={tm['pairs']} time={tm['seconds']:.3f}s")
        if tm["ok"] and tm["pairs"]:
            print(f"TMalign throughput: {tm['pairs'] / tm['seconds']:.2f} pairs/s")
        if not tm["ok"]:
            print("TMalign error preview:")
            for e in tm["errors"]:
                if e:
                    print(e)

        print("Running GTalign batch benchmark...")
        gt = bench_gtalign_batch(gtalign_bin, qdir, rdir, root / "gtalign_out", args.refs, args.residues)
        print(f"GTalign ok={gt['ok']} out_files={gt['output_files']} time={gt['seconds']:.3f}s")
        if gt["ok"]:
            print(f"GTalign effective pair rate (all-pairs searched): {total_pairs / gt['seconds']:.2f} pairs/s")
        else:
            print(f"GTalign returncode={gt['returncode']}")
        if gt["stderr_preview"]:
            print("GTalign stderr preview:")
            print(gt["stderr_preview"])
        if gt["stdout_preview"]:
            print("GTalign stdout preview:")
            print(gt["stdout_preview"])

        print("\nInterpretation:")
        print("- TMalign benchmark matches PRISM's current invocation pattern (one process per pair).")
        print("- GTalign benchmark uses batch search, which is where GTalign is designed to be fast.")
        print("- A line-for-line replacement in src/alignment.py without batching may not realize GTalign's speed advantage.")


if __name__ == "__main__":
    main()
