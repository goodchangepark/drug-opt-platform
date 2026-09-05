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
import sqlite3
import sys
from pathlib import Path
from sqlalchemy import select, text
from backend.database import SessionLocal, engine
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence, ensure_ui_schema
from backend.main import _delete_project_tree_rows

PROTECTED_PROJECT_IDS = {1, 3, 5, 300}


def run_cleanup(manifest_path: str = "validation/test_fixture_cleanup_manifest.json") -> dict:
    ensure_ui_schema(engine)

    manifest_file = Path(manifest_path)
    if not manifest_file.exists():
        raise FileNotFoundError(f"Manifest not found at {manifest_path}. Generate manifest before running cleanup.")

    with open(manifest_file, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    confirmed_test_ids = [p["id"] for p in manifest["categories"]["CONFIRMED_TEST_FIXTURE"]]
    print(f"Loaded manifest: {len(confirmed_test_ids)} confirmed test fixtures to remove.")

    db = SessionLocal()
    try:
        all_projects = db.scalars(select(Project)).all()
        all_p_ids = set(p.id for p in all_projects)
        print(f"Total projects in live database: {len(all_p_ids)}")

        # Safety sanity checks
        for prot_id in PROTECTED_PROJECT_IDS:
            assert prot_id in all_p_ids, f"CRITICAL: Protected project ID {prot_id} missing from database!"
            assert prot_id not in confirmed_test_ids, f"CRITICAL: Protected project ID {prot_id} marked for deletion in manifest!"

        to_delete = [pid for pid in confirmed_test_ids if pid in all_p_ids]
        print(f"Executing cascading deletion of {len(to_delete)} test fixture projects...")

        chunk_size = 50
        for i in range(0, len(to_delete), chunk_size):
            chunk = to_delete[i:i + chunk_size]
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

        assert remaining_ids == PROTECTED_PROJECT_IDS, f"Mismatch in remaining projects: {remaining_ids} != {PROTECTED_PROJECT_IDS}"

        # Verify DrugBank specifically
        db_compounds = db.scalars(select(Compound).where(Compound.project_id == 300)).all()
        assert len(db_compounds) in (150, 200), f"DrugBank compound count: {len(db_compounds)} not in (150, 200)"

        db_cv_ids = db.scalars(select(CompoundVersion.id).join(Compound).where(Compound.project_id == 300)).all()
        assert len(db_cv_ids) in (150, 200), f"DrugBank version count: {len(db_cv_ids)} not in (150, 200)"

        evidence_cnt = db.scalar(
            select(text("count(*)")).select_from(ExternalExperimentalEvidence).where(
                ExternalExperimentalEvidence.compound_version_id.in_(db_cv_ids)
            )
        )
        print(f"DrugBank verified: {len(db_compounds)} compounds, {len(db_cv_ids)} versions, {evidence_cnt} external evidence records.")

    finally:
        db.close()

    # Raw sqlite checks
    conn = sqlite3.connect("drug_opt.db")
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
