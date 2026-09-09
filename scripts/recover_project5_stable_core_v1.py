#!/usr/bin/env python3
"""Evidence-preserving row-level recovery of protected AMYR Project 5."""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

EXPECTED_RUN_IDS = (64,65,66,67,68,73,74,75,127,128,153,154,155,156,157,158,159,160,161,162)
EXPECTED_COMPOUND_IDS = (11,12,13,14)
EXPECTED_VERSION_IDS = (14,15,16,17)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True)
    parser.add_argument("--source-backup", required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--recovery-backup", required=True)
    parser.add_argument("--recovery-backup-sha256", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    target = Path(args.database).resolve()
    source = Path(args.source_backup).resolve()
    recovery_backup = Path(args.recovery_backup).resolve()
    if sha256(source) != args.source_sha256:
        raise RuntimeError("Project 5 source backup hash mismatch")
    if sha256(recovery_backup) != args.recovery_backup_sha256:
        raise RuntimeError("Pre-recovery safety backup hash mismatch")

    db = sqlite3.connect(target)
    db.execute("PRAGMA foreign_keys=ON")
    db.execute("ATTACH DATABASE ? AS recovery_source", (str(source),))
    try:
        if db.execute("SELECT 1 FROM projects WHERE id=5").fetchone():
            raise RuntimeError("Project 5 already exists; refusing duplicate recovery")
        for table, ids in (("compounds", EXPECTED_COMPOUND_IDS), ("compound_versions", EXPECTED_VERSION_IDS), ("prediction_runs", EXPECTED_RUN_IDS)):
            marks = ",".join("?" for _ in ids)
            existing = db.execute(f"SELECT id FROM {table} WHERE id IN ({marks})", ids).fetchall()
            if existing:
                raise RuntimeError(f"Target IDs already occupied in {table}: {existing}")
        source_runs = tuple(row[0] for row in db.execute(
            "SELECT id FROM recovery_source.prediction_runs WHERE id IN (%s) ORDER BY id" % ",".join("?" for _ in EXPECTED_RUN_IDS),
            EXPECTED_RUN_IDS,
        ))
        if source_runs != EXPECTED_RUN_IDS:
            raise RuntimeError(f"Source backup is missing expected historical runs: {source_runs}")

        db.execute("BEGIN IMMEDIATE")
        db.execute("""
          INSERT INTO projects (
            id,name,target,molecule_type,indication,mechanism_modality,description,created_at,updated_at,is_test_fixture,
            lifecycle_status,protection_policy,archived_at
          ) SELECT id,name,target,molecule_type,indication,mechanism_modality,description,created_at,updated_at,is_test_fixture,
                   'ACTIVE','PROTECTED_REAL_PROJECT',NULL
            FROM recovery_source.projects WHERE id=5
        """)
        db.execute("""
          INSERT INTO compounds (id,project_id,compound_id,cas_number,name,notes,status,current_version,created_at,updated_at)
          SELECT id,project_id,compound_id,cas_number,name,notes,status,current_version,created_at,updated_at
          FROM recovery_source.compounds WHERE project_id=5 ORDER BY id
        """)
        db.execute("""
          INSERT INTO compound_versions (
            id,compound_row_id,version_number,original_smiles,canonical_smiles,isomeric_smiles,inchi,inchikey,change_note,
            properties_json,alerts_json,assessment_json,calculation_json,svg,highlighted_svg,created_at
          ) SELECT id,compound_row_id,version_number,original_smiles,canonical_smiles,isomeric_smiles,inchi,inchikey,change_note,
                   properties_json,alerts_json,assessment_json,calculation_json,svg,highlighted_svg,created_at
            FROM recovery_source.compound_versions WHERE id IN (14,15,16,17) ORDER BY id
        """)
        db.execute("""
          INSERT INTO prediction_runs (id,version_id,stage,model_name,model_version,inputs_hash,outputs_json,provenance_json,confidence,created_at)
          SELECT id,version_id,stage,model_name,model_version,inputs_hash,outputs_json,provenance_json,confidence,created_at
          FROM recovery_source.prediction_runs
          WHERE id IN (64,65,66,67,68,73,74,75,127,128,153,154,155,156,157,158,159,160,161,162)
          ORDER BY id
        """)
        db.execute("""
          INSERT INTO scientific_mutation_audit (
            occurred_at,actor,process,operation,object_type,object_id,project_id_snapshot,transaction_id,
            before_identity_json,after_identity_json,reason,outcome
          ) VALUES (CURRENT_TIMESTAMP,'STABLE_CORE_MIGRATION','recover_project5_stable_core_v1.py','RECOVER',
                    'Project','5',5,'stable-core-v1-project5-recovery','{}',
                    json_object('project_id',5,'compound_count',4,'historical_prediction_count',20),
                    'Controlled row-level recovery from verified pre-deletion backup','COMMITTED')
        """)
        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.execute("DETACH DATABASE recovery_source")
        db.close()

    verify = sqlite3.connect(f"file:{target}?mode=ro", uri=True)
    report = {
        "recovery": "stable-core-v1-project5",
        "source_backup": str(source),
        "source_sha256": args.source_sha256,
        "pre_recovery_backup": str(recovery_backup),
        "pre_recovery_backup_sha256": args.recovery_backup_sha256,
        "project": list(verify.execute("SELECT id,name,lifecycle_status,protection_policy FROM projects WHERE id=5").fetchone()),
        "compound_ids": [row[0] for row in verify.execute("SELECT id FROM compounds WHERE project_id=5 ORDER BY id")],
        "compound_version_ids": [row[0] for row in verify.execute("SELECT id FROM compound_versions WHERE compound_row_id IN (11,12,13,14) ORDER BY id")],
        "historical_prediction_run_ids": [row[0] for row in verify.execute("SELECT id FROM prediction_runs WHERE version_id IN (14,15,16,17) ORDER BY id")],
        "immutable_history_ids": [row[0] for row in verify.execute("SELECT legacy_prediction_run_id FROM historical_predictions WHERE project_id_snapshot=5 ORDER BY legacy_prediction_run_id")],
        "later_run_ids_preserved": [row[0] for row in verify.execute("SELECT id FROM prediction_runs WHERE id IN (167,168) ORDER BY id")],
        "prediction_run_total": verify.execute("SELECT COUNT(*) FROM prediction_runs").fetchone()[0],
        "derived_cache_rows_recovered": 0,
        "integrity_check": verify.execute("PRAGMA integrity_check").fetchone()[0],
        "foreign_key_check": verify.execute("PRAGMA foreign_key_check").fetchall(),
    }
    verify.close()
    if tuple(report["compound_ids"]) != EXPECTED_COMPOUND_IDS or tuple(report["compound_version_ids"]) != EXPECTED_VERSION_IDS:
        raise RuntimeError(f"Project 5 identity verification failed: {report}")
    if tuple(report["historical_prediction_run_ids"]) != EXPECTED_RUN_IDS or tuple(report["immutable_history_ids"]) != EXPECTED_RUN_IDS:
        raise RuntimeError(f"Historical prediction verification failed: {report}")
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
