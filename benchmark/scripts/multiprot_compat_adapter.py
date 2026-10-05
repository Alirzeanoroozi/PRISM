#!/usr/bin/env python3
"""Create a non-destructive, hash-pinned MultiProt compatibility snapshot.

The pristine ``working_version/multiprot`` tree is never edited.  The
effective snapshot contains the historical algorithm modules plus a narrow
runtime boundary that replaces network PDB downloads, database writers,
HTML/mail writers, and destructive cleanup with local event records.  The
alignment, transformation-filtering, and FiberDock source files are copied
without algorithmic edits and their pristine hashes are recorded.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from pathlib import Path


BOUNDARY_MODULES = {
    "pdbDownload.py",
    "htmlWriter.py",
    "mysqlWriter.py",
    "databaseChecker.py",
    "sendMail.py",
    "checkTemplate.py",
    "templateGenerator.py",
    "startMail.py",
}
PATCHED_MODULES = (
    "mainController.py",
    "structuralAlignment.py",
    "transformationFiltering.py",
    "surfaceExtractor.py",
)
RUNTIME = r'''# Python 2/3-compatible local-only compatibility boundary.
from __future__ import print_function
import json
import os


def _event(kind, value):
    path = os.environ.get("PRISM_COMPAT_EVENT_LOG", "compat_events.jsonl")
    record = {"kind": kind, "value": str(value)}
    with open(path, "a") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def compat_cleanup(value):
    _event("cleanup_suppressed", value)


class HtmlWriter(object):
    def __init__(self, *args):
        _event("html_suppressed", args[0] if args else "")


class MysqlWriter(object):
    def __init__(self, *args):
        _event("database_write_suppressed", args[0] if args else "")


class MailSender(object):
    def __init__(self, *args):
        _event("mail_suppressed", args[0] if args else "")


class DatabaseChecker(object):
    def __init__(self, left, right):
        self.left = left
        self.right = right

    def checker(self):
        _event("database_read_suppressed", len(self.left))
        return self.left, self.right, []


class PDBdownload(object):
    """Read already-staged pair/template lists; never access the network."""
    def __init__(self, work_path):
        self.work_path = os.path.abspath(work_path)

    def PDBdownloader(self):
        pair_path = os.path.join(self.work_path, "lists", "pair_list")
        template_path = os.path.join(self.work_path, "lists", "template_list")
        left = []
        right = []
        templates = []
        if os.path.exists(pair_path):
            with open(pair_path) as handle:
                for line in handle:
                    fields = line.strip().split()
                    if len(fields) == 2 and len(fields[0]) >= 4 and len(fields[1]) >= 4:
                        left.append(self._normalize(fields[0]))
                        right.append(self._normalize(fields[1]))
        if os.path.exists(template_path):
            with open(template_path) as handle:
                templates = [line.strip()[:6] for line in handle if line.strip()]
        _event("network_download_suppressed", len(left) + len(right))
        return left, right, templates

    @staticmethod
    def _normalize(value):
        return value[:4].lower() + "".join(sorted(set(ch for ch in value[4:] if ch.isalnum())))


class TemplateChecker(object):
    """Accept only pre-existing template records; never invoke template generation."""
    def __init__(self, work_path, templates):
        self.work_path = os.path.abspath(work_path)
        self.templates = templates

    def checker(self):
        available = []
        # The historical controller changes into ``run_files`` and passes a
        # job-relative path (``../jobs/<id>``).  The checked-out template
        # manifest remains at the workspace root, not inside every job
        # directory.  Preserve that layout while making the compatibility
        # boundary resolve the same manifest the controller was given.
        default_path = ""
        candidate = self.work_path
        for _ in range(3):
            path = os.path.join(candidate, "template_default")
            if os.path.exists(path):
                default_path = path
                break
            candidate = os.path.dirname(candidate)
        if os.path.exists(default_path):
            with open(default_path) as handle:
                available = [line.strip()[:6] for line in handle if line.strip()]
        selected = [item for item in self.templates if item in available]
        if not selected:
            return [0, []]
        return [1 if len(selected) == len(available) else 2, selected]


class _Cursor(object):
    def execute(self, *args):
        _event("database_query_suppressed", args[0] if args else "")

    def fetchone(self):
        return None

    def fetchall(self):
        return []


class _Connection(object):
    def cursor(self):
        return _Cursor()

    def commit(self):
        return None

    def rollback(self):
        return None

    def close(self):
        return None


def connect(*args, **kwargs):
    _event("database_connect_suppressed", kwargs.get("db", ""))
    return _Connection()
'''


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def replace_exact(text: str, old: str, new: str, *, expected: int = 1) -> str:
    count = text.count(old)
    if count != expected:
        raise ValueError(f"compatibility patch expected {expected} occurrences, found {count}: {old!r}")
    return text.replace(old, new)


def patch_module(name: str, text: str) -> str:
    if name == "mainController.py":
        for line in (
            "from pdbDownload import PDBdownload #pdbdownloader\n",
            "from checkTemplate import TemplateChecker #checks if template exists and creates if needed\n",
            "from htmlWriter import HtmlWriter #create htmlfile\n",
            "from mysqlWriter import MysqlWriter #write to the database\n",
            "from databaseChecker import DatabaseChecker #check if the value exist or not saves time\n",
            "from sendMail import MailSender # check if user provided a mail address and sends mail if there is one\n",
        ):
            text = replace_exact(text, line, "", expected=1)
        text = replace_exact(
            text,
            "import os,sys\n",
            "import os,sys\nfrom compat_runtime import PDBdownload, TemplateChecker, HtmlWriter, MysqlWriter, DatabaseChecker, MailSender, compat_cleanup\n",
            expected=1,
        )
        return replace_exact(text, 'os.system("rm -r %s/*/" % (workPath))', 'compat_cleanup("%s/*/" % (workPath))', expected=1)
    if name == "structuralAlignment.py":
        text = replace_exact(text, "import MySQLdb as mdb\n", "import compat_runtime as mdb\nfrom compat_runtime import compat_cleanup\n", expected=1)
        text = replace_exact(text, 'os.system("rm 2_sets.res")', 'compat_cleanup("2_sets.res")', expected=1)
        text = replace_exact(text, 'os.system("rm log_multiprot.txt")', 'compat_cleanup("log_multiprot.txt")', expected=1)
        return replace_exact(text, 'os.system("rm 2_sol.res")', 'compat_cleanup("2_sol.res")', expected=1)
    if name == "transformationFiltering.py":
        return replace_exact(text, "import MySQLdb as mdb\n", "import compat_runtime as mdb\n", expected=1)
    if name == "surfaceExtractor.py":
        text = replace_exact(text, "import string,os,ConfigParser  \n", "import string,os,ConfigParser  \nfrom compat_runtime import compat_cleanup\n", expected=1)
        text = replace_exact(text, 'os.system("rm sigma.out")', 'compat_cleanup("sigma.out")', expected=1)
        text = replace_exact(text, 'os.system("rm %s" % (asafile))', 'compat_cleanup(asafile)', expected=1)
        return replace_exact(text, 'os.system("rm %s" % (logfile))', 'compat_cleanup(logfile)', expected=1)
    return text


def create_snapshot(source_root: str | Path, output_root: str | Path) -> dict[str, object]:
    source = Path(source_root).resolve()
    output = Path(output_root).resolve()
    if output == source or source in output.parents:
        raise ValueError("compatibility output must not be the source directory or inside it")
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite non-empty compatibility output: {output}")
    run_files = output / "run_files"
    run_files.mkdir(parents=True, exist_ok=True)
    source_records = []
    effective_records = []
    for path in sorted((source / "run_files").glob("*.py")):
        if path.name in BOUNDARY_MODULES:
            source_records.append({"path": str(path), "sha256": sha256_file(path), "role": "excluded_boundary_source"})
            continue
        target = run_files / path.name
        pristine = path.read_text(encoding="utf-8", errors="replace")
        effective = patch_module(path.name, pristine) if path.name in PATCHED_MODULES else pristine
        if path.name in PATCHED_MODULES:
            target.write_text(effective, encoding="utf-8")
        else:
            shutil.copy2(path, target)
        source_records.append({"path": str(path), "sha256": sha256_file(path), "role": "algorithm_or_boundary"})
        effective_records.append({"path": str(target), "sha256": sha256_file(target), "role": "effective"})
    (run_files / "compat_runtime.py").write_text(RUNTIME, encoding="utf-8")
    effective_records.append({"path": str(run_files / "compat_runtime.py"), "sha256": sha256_file(run_files / "compat_runtime.py"), "role": "compatibility_boundary"})
    for name in ("prism.py", "prism.ini", "config.inc", "template_default", "template_engin", "template_list", "progress.inc"):
        path = source / name
        if not path.is_file():
            continue
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        source_records.append({"path": str(path), "sha256": sha256_file(path), "role": "configuration"})
        effective_records.append({"path": str(target), "sha256": sha256_file(target), "role": "effective"})
    multiprot_source = source / "external_tools/multiprot"
    if multiprot_source.is_dir():
        for path in sorted(multiprot_source.rglob("*")):
            if not path.is_file() or path.name == "log_multiprot.txt":
                continue
            relative = path.relative_to(source)
            target = output / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            mode = path.stat().st_mode & 0o777
            if path.name == "multiprot.Linux" or path.suffix == ".Linux":
                mode |= 0o111
            target.chmod(mode)
            role = "executable" if path.name == "multiprot.Linux" else "multiprot_runtime"
            effective_role = "effective_executable" if path.name == "multiprot.Linux" else "effective_multiprot_runtime"
            source_records.append({"path": str(path), "sha256": sha256_file(path), "role": role})
            effective_records.append({"path": str(target), "sha256": sha256_file(target), "role": effective_role})
    manifest = output / "compatibility_manifest.json"
    value = {
        "schema_version": "multiprot-compatibility/v2",
        "pristine_root": str(source),
        "effective_root": str(output),
        "disabled_behaviors": ["database", "mail", "html", "network_download", "destructive_cleanup"],
        "scope": "standalone_multiprot_compatibility_only",
        "missing_full_pipeline_assets": ["NACCESS", "POPS", "FiberDock", "legacy_template_interface_assets"],
        "algorithm_modules_boundary_only": ["structuralAlignment.py", "transformationFiltering.py", "flexibleRefinement.py"],
        "source_records": source_records,
        "effective_records": effective_records,
        "patch_scope": list(PATCHED_MODULES),
        "excluded_generated_files": ["external_tools/multiprot/log_multiprot.txt"],
    }
    manifest.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output / "compatibility_files.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("path", "sha256", "role"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(source_records + effective_records)
    return value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        value = create_snapshot(args.source_root, args.output_root)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"wrote MultiProt compatibility snapshot to {args.output_root}")
    print(f"source_records={len(value['source_records'])} effective_records={len(value['effective_records'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
