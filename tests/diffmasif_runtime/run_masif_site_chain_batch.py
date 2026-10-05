import argparse
import os
import subprocess
from pathlib import Path
from typing import List


def run_command(args: List[str], cwd: str) -> None:
    result = subprocess.run(
        args,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(args)}\n{result.stdout}")


def slice_items(items: List[str], chunk_size: int, chunk_index: int) -> List[str]:
    start = chunk_index * chunk_size
    end = min(start + chunk_size, len(items))
    return items[start:end]


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masif-root", required=True)
    parser.add_argument("--chain-list", required=True)
    parser.add_argument("--chunk-size", type=int, required=True)
    parser.add_argument("--chunk-index", type=int, required=True)
    parser.add_argument("--log-dir", default="tests/diffmasif_runtime/slurm_logs")
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    masif_root = Path(args.masif_root).resolve()
    chain_list = [
        line.strip()
        for line in Path(args.chain_list).read_text().splitlines()
        if line.strip()
    ]
    batch = slice_items(chain_list, args.chunk_size, args.chunk_index)
    if not batch:
        print("No chains assigned to this chunk.")
        return

    masif_site_dir = masif_root / "data" / "masif_site"
    benchmark_pdb_dir = masif_site_dir / "data_preparation" / "01-benchmark_pdbs"
    pred_dir = masif_site_dir / "output" / "all_feat_3l" / "pred_data"
    log_dir = Path(args.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)

    for pdb_chain in batch:
        pdb_id, chain_id = pdb_chain.split("_")
        local_pdb = Path(f"templates/pdbs/{pdb_id.lower()}.pdb").resolve()
        prep_marker = benchmark_pdb_dir / f"{pdb_chain}.pdb"
        pred_marker = pred_dir / f"pred_{pdb_chain}.npy"
        status_path = log_dir / f"{pdb_chain}.status"
        ensure_parent(status_path)

        if args.skip_existing and prep_marker.exists() and pred_marker.exists():
            status_path.write_text("skipped_existing\n")
            continue

        try:
            run_command(
                ["./data_prepare_one.sh", "--file", str(local_pdb), pdb_chain],
                cwd=str(masif_site_dir),
            )
            run_command(
                ["./predict_site.sh", pdb_chain],
                cwd=str(masif_site_dir),
            )
            status_path.write_text("ok\n")
        except Exception as exc:
            status_path.write_text(f"failed\n{exc}\n")


if __name__ == "__main__":
    main()
