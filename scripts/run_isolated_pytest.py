#!/usr/bin/env python3
"""Run pytest and prove production scientific identities did not change."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

from production_db_fingerprint import fingerprint

ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_DB = (ROOT / "drug_opt.db").resolve()


def main() -> int:
    before = fingerprint(PRODUCTION_DB)
    test_root = Path(tempfile.mkdtemp(prefix="drugopt-pytest-wrapper-"))
    test_database = (test_root / "drugopt-test.db").resolve()
    env = os.environ.copy()
    env.pop("DRUGOPT_DATABASE_URL", None)
    env.pop("DRUGOPT_ENV", None)
    env["DRUGOPT_TEST_DATABASE_PATH"] = str(test_database)
    try:
        completed = subprocess.run([str(ROOT / ".venv/bin/pytest"), *sys.argv[1:]], cwd=ROOT, env=env)
    finally:
        shutil.rmtree(test_root, ignore_errors=True)
    after = fingerprint(PRODUCTION_DB)
    report = {
        "production_database": str(PRODUCTION_DB),
        "test_database": str(test_database),
        "physical_paths_distinct": test_database != PRODUCTION_DB,
        "before": before,
        "after": after,
        "unchanged": before["logical_fingerprint"] == after["logical_fingerprint"],
        "pytest_exit_code": completed.returncode,
    }
    report_path = Path(tempfile.gettempdir()) / "drugopt-isolated-pytest-production-fingerprint.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Production fingerprint report: {report_path}")
    if not report["physical_paths_distinct"] or not report["unchanged"]:
        print("FATAL: production database scientific identities changed during pytest", file=sys.stderr)
        return 86
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
