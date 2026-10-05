#!/usr/bin/env python3
"""Fail-closed PyRosetta runtime probe for a downstream submission."""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python-wrapper", required=True, type=Path)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    started = datetime.now(timezone.utc).isoformat()
    environment = dict(os.environ)
    environment["PRISM_PYROSETTA_BUNDLE"] = str(args.bundle)
    result = subprocess.run(
        ["bash", str(args.python_wrapper), "-c", "import pyrosetta; print(getattr(pyrosetta, '__version__', 'available'))"],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )
    payload = {
        "status": "available" if result.returncode == 0 else "unavailable",
        "available": result.returncode == 0,
        "started_at": started,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "python_wrapper": str(args.python_wrapper),
        "bundle": str(args.bundle),
        "return_code": result.returncode,
        "version": result.stdout.strip(),
        "error": result.stderr.strip()[-4000:],
        "platform": platform.platform(),
        "python": sys.version,
        "glibc": platform.libc_ver(),
        "fail_closed": True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, sort_keys=True))
    return 0 if result.returncode == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
