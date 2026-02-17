"""Convert xlsx files to csv. Usage: python convert_xlss_csv.py [file1.xlsx ...]"""

import sys
import os

try:
    import pandas as pd
except ImportError:
    print("pandas is required. Install with: pip install pandas openpyxl", file=sys.stderr)
    sys.exit(1)


def xlsx_to_csv(path: str) -> None:
    """Convert a single xlsx file to one or more csv files (one per sheet)."""
    base, _ = os.path.splitext(path)
    xl = pd.ExcelFile(path, engine="openpyxl")
    for sheet in xl.sheet_names:
        df = pd.read_excel(xl, sheet_name=sheet)
        out_path = f"{base}_{sheet}.csv" if len(xl.sheet_names) > 1 else f"{base}.csv"
        df.to_csv(out_path, index=False)
        print(out_path)


def main() -> None:
    if len(sys.argv) > 1:
        paths = [p for p in sys.argv[1:] if p.endswith(".xlsx")]
        if len(paths) != len(sys.argv) - 1:
            print("Only .xlsx paths are used; others ignored.", file=sys.stderr)
    else:
        benchmark_dir = os.path.dirname(os.path.abspath(__file__))
        paths = [os.path.join(benchmark_dir, f) for f in os.listdir(benchmark_dir) if f.endswith(".xlsx")]

    if not paths:
        print("No .xlsx files found.", file=sys.stderr)
        sys.exit(1)

    for path in paths:
        if not os.path.isfile(path):
            print(f"Skip (not a file): {path}", file=sys.stderr)
            continue
        try:
            xlsx_to_csv(path)
        except Exception as e:
            print(f"Error converting {path}: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
