#!/usr/bin/env python3
"""Controlled, idempotent Stable Core v1 migration.

Rollback is restoration of the verified pre-migration SQLite backup while the
service is stopped. No legacy scientific table is removed by this migration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

from sqlalchemy import create_engine

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.stable_core import (
    MIGRATION_ID,
    SCHEMA_VERSION,
    ensure_stable_core_schema,
    migrate_legacy_scientific_records,
    migrate_stable_core_v1_002,
    migrate_stable_core_v1_003,
    migrate_stable_core_v1_004,
    migrate_stable_core_v1_005,
)

PRODUCTION_DB = (ROOT / "drug_opt.db").resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True)
    parser.add_argument("--verified-backup")
    parser.add_argument("--verified-backup-sha256")
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    database = Path(args.database).expanduser().resolve()
    if database == PRODUCTION_DB:
        if not args.verified_backup or not args.verified_backup_sha256:
            parser.error("production migration requires --verified-backup and --verified-backup-sha256")
        backup = Path(args.verified_backup).expanduser().resolve()
        actual = sha256(backup)
        if actual != args.verified_backup_sha256:
            raise RuntimeError(f"Verified backup hash mismatch: expected {args.verified_backup_sha256}, got {actual}")
    else:
        backup = Path(args.verified_backup).expanduser().resolve() if args.verified_backup else None

    engine = create_engine(f"sqlite:///{database}", connect_args={"check_same_thread": False, "timeout": 30.0})
    ensure_stable_core_schema(engine)
    with engine.begin() as connection:
        migrated = migrate_legacy_scientific_records(connection)
        normalized = migrate_stable_core_v1_002(connection)
        provenance_corrections = migrate_stable_core_v1_003(connection)
        evidence_acceptance = migrate_stable_core_v1_004(connection)
        source_semantics = migrate_stable_core_v1_005(connection)
        connection.exec_driver_sql(f"PRAGMA user_version={SCHEMA_VERSION}")
    # Reinstall the immutable-history barriers immediately after the tightly
    # scoped provenance correction transaction.
    ensure_stable_core_schema(engine)
    report = {
        "migration_id": MIGRATION_ID,
        "database": str(database),
        "verified_backup": str(backup) if backup else None,
        "verified_backup_sha256": args.verified_backup_sha256,
        "migrated": migrated,
        "normalized": normalized,
        "provenance_corrections": provenance_corrections,
        "evidence_acceptance": evidence_acceptance,
        "source_semantics": source_semantics,
        "rollback": {
            "requires_service_stopped": True,
            "procedure": "Replace the migrated database with the verified byte-for-byte backup, then restart the existing service.",
            "warning": "Rollback discards all post-migration transactions and therefore requires an explicit recovery decision.",
        },
    }
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
