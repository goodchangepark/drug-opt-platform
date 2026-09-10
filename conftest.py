"""Global pytest isolation bootstrap for Stable Core v1."""

from __future__ import annotations

import atexit
import os
import shutil
import sqlite3
import tempfile
from pathlib import Path
from urllib.parse import unquote

REPOSITORY_ROOT = Path(__file__).resolve().parent
PRODUCTION_DB = (REPOSITORY_ROOT / "drug_opt.db").resolve()
_REQUESTED_TEST_DB = os.environ.get("DRUGOPT_TEST_DATABASE_PATH")
if _REQUESTED_TEST_DB:
    TEST_DB = Path(_REQUESTED_TEST_DB).expanduser().resolve()
    _TEST_ROOT = TEST_DB.parent
    _TEST_ROOT.mkdir(parents=True, exist_ok=True)
else:
    _TEST_ROOT = Path(tempfile.mkdtemp(prefix="drugopt-pytest-"))
    TEST_DB = (_TEST_ROOT / "drugopt-test.db").resolve()

if PRODUCTION_DB.exists():
    shutil.copy2(PRODUCTION_DB, TEST_DB)
if TEST_DB == PRODUCTION_DB or (TEST_DB.exists() and PRODUCTION_DB.exists() and TEST_DB.samefile(PRODUCTION_DB)):
    raise RuntimeError("FATAL: TEST DATABASE RESOLVES TO PRODUCTION DATABASE")

os.environ["DRUGOPT_ENV"] = "test"
os.environ["DRUGOPT_DATABASE_URL"] = f"sqlite:///{TEST_DB}"
os.environ["DRUGOPT_TEST_DATABASE_PATH"] = str(TEST_DB)
os.environ["DRUGOPT_PRODUCTION_DATABASE_PATH"] = str(PRODUCTION_DB)

# Bootstrap only the private copy before test modules import global application
# objects. The production file is never bound to this engine.
from backend.database import engine as _test_engine  # noqa: E402
from backend.stable_core import (  # noqa: E402
    ensure_stable_core_schema,
    migrate_legacy_scientific_records,
    migrate_stable_core_v1_002,
    migrate_stable_core_v1_003,
    migrate_stable_core_v1_004,
    migrate_stable_core_v1_005,
    migrate_stable_core_v1_006,
)

ensure_stable_core_schema(_test_engine)

# Block direct sqlite3 file access as well as misconfigured SQLAlchemy engines.
# This closes the historical escape hatch where a test bypassed SessionLocal
# and opened repository-root drug_opt.db by path.
_ORIGINAL_SQLITE_CONNECT = sqlite3.connect


def _guarded_sqlite_connect(database, *args, **kwargs):
    raw = os.fspath(database) if isinstance(database, os.PathLike) else str(database)
    candidate = None
    if raw != ":memory:":
        if raw.startswith("file:"):
            candidate = Path(unquote(raw[5:].split("?", 1)[0])).expanduser().resolve()
        elif not raw.startswith(":"):
            candidate = Path(raw).expanduser().resolve()
    if candidate is not None and (
        candidate == PRODUCTION_DB
        or (candidate.exists() and PRODUCTION_DB.exists() and candidate.samefile(PRODUCTION_DB))
    ):
        raise RuntimeError("FATAL: TEST DATABASE RESOLVES TO PRODUCTION DATABASE")
    return _ORIGINAL_SQLITE_CONNECT(database, *args, **kwargs)


sqlite3.connect = _guarded_sqlite_connect
with _test_engine.begin() as _connection:
    migrate_legacy_scientific_records(_connection)
    migrate_stable_core_v1_002(_connection)
    migrate_stable_core_v1_003(_connection)
    migrate_stable_core_v1_004(_connection)
    migrate_stable_core_v1_005(_connection)
    migrate_stable_core_v1_006(_connection)
ensure_stable_core_schema(_test_engine)


def pytest_report_header(config):
    return [
        "Drug-OPT environment: test",
        f"Drug-OPT test DB: {TEST_DB}",
        f"Drug-OPT production DB (fingerprint only): {PRODUCTION_DB}",
    ]


@atexit.register
def _remove_test_root() -> None:
    shutil.rmtree(_TEST_ROOT, ignore_errors=True)
