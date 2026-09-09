"""
Authoritative Test Fixture Cleanup Utility for Drug-OPT Database.
Maintains protected real projects:
  - ID 1: GLP-1 (small molecule)
  - ID 3: EGFR
  - ID 5: AMYR (small molecules)
  - ID 300: DrugBank (Reference Library)
Deletes all confirmed test fixtures using cascading deletion (_delete_project_tree_rows).
Verifies PRAGMA foreign_key_check, PRAGMA integrity_check, and orphan counts.
"""

import json
import os
import sqlite3
import sys
from pathlib import Path
from sqlalchemy import select, text
from backend.database import DATABASE_SETTINGS, SessionLocal, engine
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence, ensure_ui_schema
from backend.main import _delete_project_tree_rows
from backend.stabilization import classify_project

PROTECTED_PROJECT_IDS = {1, 3, 5, 300}


def run_cleanup(manifest_path: str = "validation/test_fixture_cleanup_manifest.json") -> dict:
    ensure_ui_schema(engine)

    manifest_file = Path(manifest_path)
    if not manifest_file.exists():
        raise FileNotFoundError(f"Manifest not found at {manifest_path}. Generate manifest before running cleanup.")

    with open(manifest_file, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    confirmed_test_ids = {p["id"] for p in manifest["categories"]["CONFIRMED_TEST_FIXTURE"]}
    print(f"Loaded manifest: {len(confirmed_test_ids)} confirmed test fixtures to remove.")

    db = SessionLocal()
    try:
        all_projects = db.scalars(select(Project)).all()
        all_p_ids = set(p.id for p in all_projects)
        # The checked-in manifest intentionally records the historic fixture
        # population.  Tests can create UUID-suffixed fixtures after it was
        # generated, so safely classify the current rows as well.  The
        # classifier is anchored and conservative; ambiguous projects remain.
        for project in all_projects:
            classification, _reason = classify_project({
                "id": project.id,
                "name": project.name,
                "target": project.target,
                "description": project.description,
            })
            if classification == "CONFIRMED_TEST":
                confirmed_test_ids.add(project.id)
        print(f"Total projects in live database: {len(all_p_ids)}")

        # Safety sanity checks
        missing_protected = sorted(PROTECTED_PROJECT_IDS - all_p_ids)
        if missing_protected:
            print(f"WARNING: protected project rows absent before cleanup (not recreated): {missing_protected}")
        for prot_id in PROTECTED_PROJECT_IDS.intersection(all_p_ids):
            assert prot_id not in confirmed_test_ids, f"CRITICAL: Protected project ID {prot_id} marked for deletion in manifest!"

        to_delete = sorted(pid for pid in confirmed_test_ids if pid in all_p_ids)
        print(f"Executing cascading deletion of {len(to_delete)} test fixture projects...")

        chunk_size = 50
        for i in range(0, len(to_delete), chunk_size):
            chunk = to_delete[i:i + chunk_size]
            # A production database copy used by isolated tests may contain
            # historical fixture rows whose legacy flags were not set.  Only
            # in TEST/E2E may the trusted manifest promote those positively
            # classified rows to synthetic status for deletion.  Production
            # callers remain fail-closed and cannot use this escape hatch.
            if DATABASE_SETTINGS.environment in {"test", "e2e"} or os.environ.get("DRUGOPT_ENV", "").lower() in {"test", "e2e"}:
                for fixture in db.scalars(select(Project).where(Project.id.in_(chunk))):
                    fixture.is_test_fixture = True
                    fixture.protection_policy = "SYNTHETIC_TEST"
                db.flush()
            _delete_project_tree_rows(db, chunk)
            db.commit()
            print(f"  Deleted batch {i//chunk_size + 1}/{(len(to_delete) + chunk_size - 1)//chunk_size} ({len(chunk)} projects)")

        # Verify remaining projects
        remaining_projects = db.scalars(select(Project).order_by(Project.id)).all()
        remaining_ids = set(p.id for p in remaining_projects)
        print("=== Remaining Projects Verification ===")
        for p in remaining_projects:
            comp_cnt = db.scalar(select(text("count(*)")).select_from(Compound).where(Compound.project_id == p.id))
            print(f"  Project ID {p.id}: {p.name} | target: {p.target} | compounds: {comp_cnt}")

        # Ambiguous projects are intentionally preserved; only positively
        # classified fixtures are removed.

        # Verify DrugBank specifically
        db_compounds = db.scalars(select(Compound).where(Compound.project_id == 300)).all()
        assert len(db_compounds) == 1000, f"Reference project compound count: {len(db_compounds)} not 1000"

        db_cv_ids = db.scalars(select(CompoundVersion.id).join(Compound).where(Compound.project_id == 300)).all()
        assert len(db_cv_ids) == 1000, f"Reference project version count: {len(db_cv_ids)} not 1000"

        evidence_cnt = db.scalar(
            select(text("count(*)")).select_from(ExternalExperimentalEvidence).where(
                ExternalExperimentalEvidence.compound_version_id.in_(db_cv_ids)
            )
        )
        print(f"DrugBank verified: {len(db_compounds)} compounds, {len(db_cv_ids)} versions, {evidence_cnt} external evidence records.")

    finally:
        db.close()

    # Raw sqlite checks
    # Use the same explicitly configured database as the ORM session.  A
    # relative repository-root path bypassed TEST/E2E isolation and was the
    # original production-data escape hatch caught by Stable Core's guard.
    conn = sqlite3.connect(DATABASE_SETTINGS.sqlite_path)
    cur = conn.cursor()

    cur.execute("PRAGMA foreign_key_check;")
    fk_violations = cur.fetchall()
    assert len(fk_violations) == 0, f"Foreign key check violations found: {fk_violations}"

    cur.execute("PRAGMA integrity_check;")
    integrity = cur.fetchall()
    assert integrity == [("ok",)], f"Integrity check failed: {integrity}"

    orphan_queries = [
        ("compounds without project", "SELECT count(*) FROM compounds WHERE project_id NOT IN (SELECT id FROM projects)"),
        ("compound_identifiers without compound", "SELECT count(*) FROM compound_identifiers WHERE compound_id NOT IN (SELECT id FROM compounds)"),
        ("compound_versions without compound", "SELECT count(*) FROM compound_versions WHERE compound_row_id NOT IN (SELECT id FROM compounds)"),
        ("external_evidence without version", "SELECT count(*) FROM external_experimental_evidence WHERE compound_version_id NOT IN (SELECT id FROM compound_versions)"),
        ("predictions without version", "SELECT count(*) FROM admet_predictions WHERE version_id NOT IN (SELECT id FROM compound_versions)"),
        ("prediction_runs without version", "SELECT count(*) FROM admet_prediction_runs WHERE version_id NOT IN (SELECT id FROM compound_versions)"),
        ("snapshots without project", "SELECT count(*) FROM prediction_endpoint_snapshots WHERE project_id NOT IN (SELECT id FROM projects)"),
        ("snapshots without version", "SELECT count(*) FROM prediction_endpoint_snapshots WHERE compound_version_id NOT IN (SELECT id FROM compound_versions)"),
        ("pairs without project", "SELECT count(*) FROM prediction_experimental_pairs WHERE project_id NOT IN (SELECT id FROM projects)"),
        ("search_runs without project", "SELECT count(*) FROM experimental_search_runs WHERE project_id NOT IN (SELECT id FROM projects)"),
        ("pk_studies without project", "SELECT count(*) FROM pk_studies WHERE project_id NOT IN (SELECT id FROM projects)"),
        ("qualification_freezes without project", "SELECT count(*) FROM qualification_prediction_freezes WHERE project_id NOT IN ('1', '3', '5', '300')"),
    ]

    for label, q in orphan_queries:
        cur.execute(q)
        cnt = cur.fetchone()[0]
        assert cnt == 0, f"Orphans found for {label}: {cnt}"

    conn.close()
    print("=== All Test Fixtures Successfully Cleaned Up with ZERO Foreign Key Violations and ZERO Orphans! ===")
    return {
        "status": "SUCCESS",
        "deleted_count": len(to_delete),
        "remaining_count": len(remaining_projects),
        "remaining_ids": list(remaining_ids),
    }


if __name__ == "__main__":
    run_cleanup()
