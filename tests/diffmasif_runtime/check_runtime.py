import argparse
import json
import os
import shutil
from typing import Dict, List


EXPECTED_ENV_VARS = [
    "APBS_BIN",
    "MULTIVALUE_BIN",
    "PDB2PQR_BIN",
    "REDUCE_HET_DICT",
    "PYMESH_PATH",
    "MSMS_BIN",
]

EXPECTED_PATH_HINTS = [
    "source/default_config/chemistry.py",
    "source",
]


def check_exists(path: str) -> bool:
    return os.path.exists(path)


def resolve_executable(path_or_name: str) -> str:
    if os.path.isabs(path_or_name) or os.path.sep in path_or_name:
        return path_or_name if os.path.exists(path_or_name) else ""
    return shutil.which(path_or_name) or ""


def inspect_env_var(name: str) -> Dict[str, object]:
    value = os.environ.get(name, "")
    if not value:
        return {"name": name, "set": False, "resolved": "", "exists": False}
    resolved = resolve_executable(value)
    exists = bool(resolved) or os.path.exists(value)
    return {
        "name": name,
        "set": True,
        "value": value,
        "resolved": resolved,
        "exists": exists,
    }


def inspect_masif_root(masif_root: str) -> Dict[str, object]:
    report = {
        "path": masif_root,
        "exists": os.path.isdir(masif_root),
        "hints": [],
    }
    if report["exists"]:
        for rel_path in EXPECTED_PATH_HINTS:
            full_path = os.path.join(masif_root, rel_path)
            report["hints"].append(
                {
                    "path": rel_path,
                    "exists": os.path.exists(full_path),
                }
            )
    return report


def detect_python_modules() -> Dict[str, bool]:
    modules = {}
    for module_name in ["Bio", "open3d", "pymesh"]:
        try:
            __import__(module_name)
            modules[module_name] = True
        except Exception:
            modules[module_name] = False
    return modules


def summarize_failures(report: Dict[str, object]) -> List[str]:
    failures = []
    if not report["masif_root"]["exists"]:
        failures.append("MaSIF root path does not exist")
    if not report["pdb_exists"]:
        failures.append("Input PDB does not exist")
    if not report["python_modules"].get("Bio", False):
        failures.append("Biopython is not importable in the current python")
    if not report["python_modules"].get("pymesh", False):
        failures.append("PyMesh is not importable in the current python")
    missing_env = [
        item["name"] for item in report["env_vars"]
        if not item["set"] or not item["exists"]
    ]
    if missing_env:
        failures.append("Missing or unresolved env vars: " + ", ".join(missing_env))
    return failures


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masif-root", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    report = {
        "masif_root": inspect_masif_root(args.masif_root),
        "pdb": args.pdb,
        "pdb_exists": os.path.exists(args.pdb),
        "env_vars": [inspect_env_var(name) for name in EXPECTED_ENV_VARS],
        "python_modules": detect_python_modules(),
    }
    report["failures"] = summarize_failures(report)
    report["ready"] = len(report["failures"]) == 0

    if args.output:
        parent = os.path.dirname(args.output)
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(args.output, "w") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)

    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
