"""Environment-scoped database configuration for Drug-OPT.

Stable Core v1 makes the database target explicit and fails closed whenever a
destructive test/E2E process could resolve to the production SQLite file.
The module-level ``engine``/``SessionLocal`` names remain compatibility
adapters; their target is selected once, before application modules import.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from sqlalchemy import create_engine, event
from sqlalchemy.engine import Engine, make_url
from sqlalchemy.orm import DeclarativeBase, sessionmaker

ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_DATABASE_PATH = (ROOT / "drug_opt.db").resolve()
VALID_ENVIRONMENTS = {"production", "test", "e2e", "development"}


class UnsafeDatabaseConfiguration(RuntimeError):
    """Raised before a process can use production storage in a test mode."""


@dataclass(frozen=True)
class DatabaseSettings:
    environment: str
    url: str
    sqlite_path: Path | None
    production_path: Path


def _sqlite_path(url: str) -> Path | None:
    parsed = make_url(url)
    if not parsed.drivername.startswith("sqlite") or not parsed.database or parsed.database == ":memory:":
        return None
    return Path(unquote(parsed.database)).expanduser().resolve()


def _same_existing_file(left: Path, right: Path) -> bool:
    if left == right:
        return True
    try:
        return left.exists() and right.exists() and left.samefile(right)
    except OSError:
        return False


def load_database_settings(environ: dict[str, str] | None = None) -> DatabaseSettings:
    env = environ if environ is not None else os.environ
    environment = str(env.get("DRUGOPT_ENV", "production")).strip().lower()
    if environment not in VALID_ENVIRONMENTS:
        raise UnsafeDatabaseConfiguration(
            f"Unsupported DRUGOPT_ENV={environment!r}; expected one of {sorted(VALID_ENVIRONMENTS)}"
        )
    default_url = f"sqlite:///{PRODUCTION_DATABASE_PATH}"
    url = str(env.get("DRUGOPT_DATABASE_URL", default_url)).strip()
    if not url:
        raise UnsafeDatabaseConfiguration("DRUGOPT_DATABASE_URL must not be empty")
    sqlite_path = _sqlite_path(url)

    pytest_loaded = "pytest" in sys.modules or bool(env.get("PYTEST_CURRENT_TEST"))
    destructive_mode = environment in {"test", "e2e"} or pytest_loaded
    if destructive_mode:
        if environment not in {"test", "e2e"}:
            raise UnsafeDatabaseConfiguration(
                "FATAL: pytest/E2E execution requires DRUGOPT_ENV=test or DRUGOPT_ENV=e2e"
            )
        if sqlite_path is None:
            raise UnsafeDatabaseConfiguration(
                "FATAL: TEST/E2E DATABASE must be an explicit file-backed SQLite database"
            )
        if _same_existing_file(sqlite_path, PRODUCTION_DATABASE_PATH):
            raise UnsafeDatabaseConfiguration(
                "FATAL: TEST DATABASE RESOLVES TO PRODUCTION DATABASE"
            )
    elif environment == "production" and sqlite_path and sqlite_path != PRODUCTION_DATABASE_PATH:
        if not Path(unquote(make_url(url).database or "")).is_absolute():
            raise UnsafeDatabaseConfiguration("Production SQLite URL must use an absolute path")

    return DatabaseSettings(environment=environment, url=url, sqlite_path=sqlite_path, production_path=PRODUCTION_DATABASE_PATH)


def create_database_engine(settings: DatabaseSettings) -> Engine:
    connect_args = {"check_same_thread": False, "timeout": 30.0} if settings.url.startswith("sqlite") else {}
    db_engine = create_engine(settings.url, connect_args=connect_args)
    if settings.url.startswith("sqlite"):
        @event.listens_for(db_engine, "connect")
        def set_sqlite_pragma(dbapi_connection, connection_record):
            cursor = dbapi_connection.cursor()
            cursor.execute("PRAGMA journal_mode=WAL;")
            cursor.execute("PRAGMA busy_timeout=30000;")
            cursor.execute("PRAGMA foreign_keys=ON;")
            cursor.close()
    return db_engine


DATABASE_SETTINGS = load_database_settings()
DATABASE_URL = DATABASE_SETTINGS.url
engine = create_database_engine(DATABASE_SETTINGS)
SessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False)


class Base(DeclarativeBase):
    pass


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
