#!/usr/bin/env python3
"""Run a browser/API E2E command against a private production-shaped DB copy."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import json
from pathlib import Path

from production_db_fingerprint import fingerprint

ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_DB = (ROOT / "drug_opt.db").resolve()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8876)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("an E2E command is required after --")
    root = Path(tempfile.mkdtemp(prefix="drugopt-e2e-"))
    database = (root / "drugopt-e2e.db").resolve()
    shutil.copy2(PRODUCTION_DB, database)
    if database == PRODUCTION_DB or database.samefile(PRODUCTION_DB):
        raise RuntimeError("FATAL: TEST DATABASE RESOLVES TO PRODUCTION DATABASE")
    env = os.environ.copy()
    env.update({
        "DRUGOPT_ENV": "e2e",
        "DRUGOPT_DATABASE_URL": f"sqlite:///{database}",
        "DRUGOPT_E2E_DATABASE_PATH": str(database),
        "DRUGOPT_PRODUCTION_DATABASE_PATH": str(PRODUCTION_DB),
        "DRUGOPT_E2E_BASE_URL": f"http://127.0.0.1:{args.port}",
    })
    before = fingerprint(PRODUCTION_DB)
    server_log = root / "uvicorn.log"
    return_code = 1
    with server_log.open("w", encoding="utf-8") as log:
        server = subprocess.Popen(
            [str(ROOT / ".venv/bin/uvicorn"), "backend.main:app", "--host", "127.0.0.1", "--port", str(args.port)],
            cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
        )
        try:
            health = f"http://127.0.0.1:{args.port}/api/health"
            for _ in range(180):
                if server.poll() is not None:
                    raise RuntimeError(f"isolated E2E server exited early; log: {server_log}")
                try:
                    with urllib.request.urlopen(health, timeout=1) as response:
                        if response.status == 200:
                            break
                except Exception:
                    time.sleep(1)
            else:
                raise RuntimeError(f"isolated E2E server did not become healthy; log: {server_log}")
            print(f"E2E environment: e2e\nE2E DB: {database}\nProduction DB: {PRODUCTION_DB}", flush=True)
            return_code = subprocess.run(command, cwd=ROOT, env=env).returncode
        finally:
            server.terminate()
            try:
                server.wait(timeout=20)
            except subprocess.TimeoutExpired:
                server.kill()
    after = fingerprint(PRODUCTION_DB)
    report = {
        "production_database": str(PRODUCTION_DB),
        "e2e_database": str(database),
        "physical_paths_distinct": database != PRODUCTION_DB,
        "before": before,
        "after": after,
        "unchanged": before["logical_fingerprint"] == after["logical_fingerprint"],
        "e2e_exit_code": return_code,
    }
    report_path = Path(tempfile.gettempdir()) / "drugopt-isolated-e2e-production-fingerprint.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Production fingerprint report: {report_path}")
    shutil.rmtree(root, ignore_errors=True)
    if not report["physical_paths_distinct"] or not report["unchanged"]:
        print("FATAL: production database scientific identities changed during E2E", file=sys.stderr)
        return 86
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
