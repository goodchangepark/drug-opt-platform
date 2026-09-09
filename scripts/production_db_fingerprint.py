#!/usr/bin/env python3
"""Read-only, identity-focused production DB fingerprint."""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

KEY_TABLES = (
    "projects", "compounds", "compound_versions", "prediction_runs",
    "external_experimental_evidence", "admet_measurements", "prediction_endpoint_snapshots",
    "historical_predictions", "experimental_observations", "current_prediction_snapshots",
    "scientific_mutation_audit", "pk_studies", "pk_observations", "pk_nca_results",
    "pk_parameter_sets", "pk_simulation_runs", "ivive_runs",
)


def _digest_rows(connection: sqlite3.Connection, table: str) -> dict:
    columns = [row[1] for row in connection.execute(f"PRAGMA table_info({table})")]
    order = "id" if "id" in columns else ",".join(columns)
    digest = hashlib.sha256()
    count = 0
    for row in connection.execute(f"SELECT * FROM {table} ORDER BY {order}"):
        digest.update(json.dumps(row, default=str, ensure_ascii=False, separators=(",", ":")).encode())
        digest.update(b"\n")
        count += 1
    return {"count": count, "sha256": digest.hexdigest()}


def fingerprint(path: Path) -> dict:
    resolved = path.expanduser().resolve()
    connection = sqlite3.connect(f"file:{resolved}?mode=ro", uri=True)
    try:
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        result = {
            "database_path": str(resolved),
            "database_inode": resolved.stat().st_ino,
            "tables": {name: _digest_rows(connection, name) for name in KEY_TABLES if name in tables},
            "project_identities": [list(row) for row in connection.execute("SELECT id,name FROM projects ORDER BY id")],
            "prediction_run_ids": [row[0] for row in connection.execute("SELECT id FROM prediction_runs ORDER BY id")],
            "integrity_check": [row[0] for row in connection.execute("PRAGMA integrity_check")],
            "foreign_key_check": [list(row) for row in connection.execute("PRAGMA foreign_key_check")],
        }
    finally:
        connection.close()
    canonical = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    result["logical_fingerprint"] = hashlib.sha256(canonical).hexdigest()
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default=str(Path(__file__).resolve().parents[1] / "drug_opt.db"))
    parser.add_argument("--output")
    args = parser.parse_args()
    result = fingerprint(Path(args.database))
    rendered = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        Path(args.output).write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
