import sqlite3
import sys
from sqlalchemy import select, text
from backend.database import SessionLocal, engine
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence, ensure_ui_schema
from backend.main import _delete_project_tree_rows

def main():
    print("=== Step 1: Ensure UI schema migration for is_test_fixture ===")
    ensure_ui_schema(engine)

    db = SessionLocal()
    try:
        # Check current projects
        real_ids = {1, 3, 5, 300}
        all_projects = db.scalars(select(Project)).all()
        print(f"Total projects in DB: {len(all_projects)}")
        
        test_project_ids = [p.id for p in all_projects if p.id not in real_ids]
        print(f"Test projects to clean up: {len(test_project_ids)}")

        # Clean up test projects in batches
        print("=== Step 2: Cascade cleanup of test fixture projects ===")
        chunk_size = 50
        for i in range(0, len(test_project_ids), chunk_size):
            chunk = test_project_ids[i:i + chunk_size]
            _delete_project_tree_rows(db, chunk)
            db.commit()
            print(f"Deleted chunk {i//chunk_size + 1}/{(len(test_project_ids) + chunk_size - 1)//chunk_size} ({len(chunk)} projects)")

        # Verify real projects
        remaining_projects = db.scalars(select(Project).order_by(Project.id)).all()
        print("=== Step 3: Verify Remaining Protected Projects ===")
        for p in remaining_projects:
            comp_cnt = db.scalar(select(text("count(*)")).select_from(Compound).where(Compound.project_id == p.id))
            print(f"Project ID {p.id}: {p.name} | compounds: {comp_cnt} | is_test_fixture: {p.is_test_fixture}")
        
        assert set(p.id for p in remaining_projects) == real_ids, f"Mismatch in remaining projects: {[p.id for p in remaining_projects]}"

        # Verify DrugBank specifically
        drugbank_compounds = db.scalars(select(Compound).where(Compound.project_id == 300)).all()
        assert len(drugbank_compounds) == 80, f"DrugBank compound count: {len(drugbank_compounds)} != 80"
        
        db_cv_ids = db.scalars(select(CompoundVersion.id).join(Compound).where(Compound.project_id == 300)).all()
        assert len(db_cv_ids) == 80, f"DrugBank version count: {len(db_cv_ids)} != 80"
        
        evidence_cnt = db.scalar(select(text("count(*)")).select_from(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id.in_(db_cv_ids)))
        print(f"DrugBank verified: 80 compounds, 80 versions, {evidence_cnt} external evidence records.")

    finally:
        db.close()

    # Step 4: Raw sqlite foreign key and orphan check
    print("=== Step 4: PRAGMA foreign_key_check and Orphan Check ===")
    conn = sqlite3.connect("drug_opt.db")
    cur = conn.cursor()
    cur.execute("PRAGMA foreign_key_check;")
    violations = cur.fetchall()
    print(f"PRAGMA foreign_key_check violations: {len(violations)}")
    assert len(violations) == 0, f"Violations found: {violations}"

    cur.execute("PRAGMA integrity_check;")
    integrity = cur.fetchall()
    print(f"PRAGMA integrity_check: {integrity}")
    assert integrity == [("ok",)], f"Integrity check failed: {integrity}"

    orphan_queries = [
        ("compounds without project", "SELECT count(*) FROM compounds WHERE project_id NOT IN (SELECT id FROM projects)"),
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
        print(f"Orphan count [{label}]: {cnt}")
        assert cnt == 0, f"Orphans found for {label}: {cnt}"

    conn.close()
    print("=== All Test Fixtures Successfully Cleaned Up with ZERO Orphans! ===")

if __name__ == "__main__":
    main()
