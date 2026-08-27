#!/usr/bin/env python3
"""Audit SQLite integrity and repair only orphan rows tied to confirmed test projects."""

from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "drug_opt.db"
CLEANUP = ROOT / "validation" / "test_projects_cleanup.json"
OUTPUT = ROOT / "validation" / "stage5b4_stabilization_database_audit.json"


def main() -> int:
    cleanup = json.loads(CLEANUP.read_text(encoding="utf-8"))
    deleted_projects = {int(row["project_id"]) for row in cleanup["projects"] if row["deletion_status"] == "DELETED"}
    connection = sqlite3.connect(DB)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys=ON")
    before = [tuple(row) for row in connection.execute("PRAGMA foreign_key_check")]
    repairs = []
    for table, rowid, parent, foreign_key_index in before:
        columns = {row[1] for row in connection.execute(f'PRAGMA table_info("{table}")')}
        if "project_id" not in columns:
            continue
        row = connection.execute(f'SELECT project_id FROM "{table}" WHERE rowid=?', (rowid,)).fetchone()
        legacy_stage5b3_fixture = False
        if row and table == "pk_translational_snapshots" and "frozen_inputs" in columns:
            evidence = connection.execute(
                'SELECT frozen_inputs FROM "pk_translational_snapshots" WHERE rowid=?', (rowid,)
            ).fetchone()[0]
            legacy_stage5b3_fixture = all(marker in evidence for marker in (
                "Mouse IV Bolus PK", "Rat IV Bolus PK", "Dog IV Bolus PK", "Monkey IV Bolus PK",
            ))
        if row and (int(row["project_id"]) in deleted_projects or legacy_stage5b3_fixture):
            connection.execute(f'DELETE FROM "{table}" WHERE rowid=?', (rowid,))
            repairs.append({"table": table, "rowid": rowid, "missing_parent_table": parent,
                            "foreign_key_index": foreign_key_index,
                            "reason": ("orphan belonged to confirmed deleted test project" if not legacy_stage5b3_fixture else
                                       "legacy orphan contains the exact four-species Stage 5B-3 browser fixture study markers"),
                            "project_id": int(row["project_id"])})
    connection.commit()
    integrity = [row[0] for row in connection.execute("PRAGMA integrity_check")]
    after = [tuple(row) for row in connection.execute("PRAGMA foreign_key_check")]
    important = ["projects", "compounds", "compound_versions", "activity_measurements", "admet_predictions",
                 "optimization_runs", "pk_studies", "pk_observations", "pk_nca_results", "pk_ivive_runs",
                 "pk_simulation_runs", "pk_translational_snapshots", "pk_human_prediction_snapshots"]
    existing = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    counts = {table: connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
              for table in important if table in existing}
    connection.close()
    payload = {"timestamp": datetime.now(timezone.utc).isoformat(), "database": DB.name,
               "integrity_check": integrity, "foreign_key_violations_before": before,
               "confirmed_test_orphan_repairs": repairs, "foreign_key_violations_after": after,
               "major_entity_counts": counts,
               "status": "PASS" if integrity == ["ok"] and not after else "FAIL"}
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "integrity": integrity,
                      "fk_before": len(before), "repaired": len(repairs), "fk_after": len(after)}, sort_keys=True))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
