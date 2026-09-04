import sys
import sqlite3
from fastapi.testclient import TestClient
from sqlalchemy import select, text
from backend.main import app, _delete_project_tree_rows
from backend.database import SessionLocal, engine
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence, ensure_ui_schema

def run_verification():
    print("=================================================================")
    print(" DRUG-OPT PRODUCTION DATA INTEGRITY & WORKFLOW RECOVERY v3.3.1  ")
    print("=================================================================")

    ensure_ui_schema(engine)
    client = TestClient(app)

    db = SessionLocal()
    test_project_id = None
    try:
        # Create single isolated test fixture project
        print("\n--- Initializing Isolated Test Fixture Project ---")
        proj_res = client.post("/api/projects", json={
            "name": "E2E_Data_Integrity_Verification_Fixture",
            "target": "TEST_TARGET",
            "molecule_type": "Small Molecule",
            "indication": "Temporary Verification Fixture",
            "description": "Transient project for E2E integrity validation v3.3.1",
            "is_test_fixture": True
        })
        assert proj_res.status_code == 201, f"Failed to create test project: {proj_res.text}"
        test_project = proj_res.json()
        test_project_id = test_project["id"]
        print(f"Created test fixture project ID: {test_project_id} (is_test_fixture={test_project.get('is_test_fixture')})")

        # TEST A: CAS-based Compound Addition
        print("\n--- Test A: CAS-based Compound Input & Auto-Resolution ---")
        cas_res = client.post("/api/structure/resolve-cas", json={"cas_number": "15687-27-1"})
        assert cas_res.status_code == 200, f"CAS resolution failed: {cas_res.text}"
        cas_data = cas_res.json()
        assert cas_data.get("found") is True
        assert "CC(C)C" in cas_data.get("smiles", "")
        assert cas_data.get("svg") is not None
        print(f"Resolved CAS 15687-27-1 -> SMILES: {cas_data['smiles']}")
        print(f"2D Depiction SVG verified: {len(cas_data['svg'])} bytes")

        # Create compound using CAS without SMILES to verify backend auto-resolution
        c_a_res = client.post(f"/api/projects/{test_project_id}/compounds", json={
            "name": "Ibuprofen",
            "compound_id": "TEST-IBU-001",
            "cas_number": "15687-27-1",
            "smiles": "",
            "notes": "Added via CAS auto-resolution",
            "calculate": True
        })
        assert c_a_res.status_code == 201, f"Create compound via CAS failed: {c_a_res.text}"
        c_a = c_a_res.json()
        ibu_row_id = c_a["row_id"]
        assert c_a["version"] is not None, "Version was not created by CAS auto-resolution"
        assert c_a["version"]["canonical_smiles"] == cas_data["smiles"]
        assert c_a["version"]["svg"] is not None
        print(f"Test A passed: Compound {c_a['compound_id']} saved with version v{c_a['current_version']} and SVG depiction.")

        # TEST B: SMILES direct input & live validation
        print("\n--- Test B: SMILES Direct Input & Live Validation ---")
        aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
        val_res = client.post("/api/structure/validate", json={"smiles": aspirin_smiles})
        assert val_res.status_code == 200, f"SMILES validation failed: {val_res.text}"
        val_data = val_res.json()
        assert val_data.get("valid") is True
        assert val_data.get("svg") is not None
        print(f"Validated SMILES '{aspirin_smiles}' -> MW: {val_data['properties']['molecular_weight']}, SVG present.")

        c_b_res = client.post(f"/api/projects/{test_project_id}/compounds", json={
            "name": "Aspirin",
            "compound_id": "TEST-ASP-002",
            "cas_number": "50-78-2",
            "smiles": aspirin_smiles,
            "notes": "Added via direct SMILES entry",
            "calculate": True
        })
        assert c_b_res.status_code == 201, f"Create compound via SMILES failed: {c_b_res.text}"
        c_b = c_b_res.json()
        assert c_b["version"] is not None
        assert c_b["version"]["canonical_smiles"] == val_data["identity"]["canonical_smiles"]
        print(f"Test B passed: Compound {c_b['compound_id']} saved with version v{c_b['current_version']}.")

        # TEST C: Ketcher Structure Drawing Simulation
        print("\n--- Test C: Structure Drawing Simulation (Ketcher Export) ---")
        phenol_smiles = "c1ccccc1O"
        val_c = client.post("/api/structure/validate", json={"smiles": phenol_smiles})
        assert val_c.status_code == 200
        c_c_res = client.post(f"/api/projects/{test_project_id}/compounds", json={
            "name": "Phenol",
            "compound_id": "TEST-PHE-003",
            "smiles": phenol_smiles,
            "notes": "Added via Ketcher structure drawing",
            "calculate": True
        })
        assert c_c_res.status_code == 201
        c_c = c_c_res.json()
        assert c_c["version"]["canonical_smiles"] == val_c.json()["identity"]["canonical_smiles"]
        print(f"Test C passed: Structure drawing export saved as {c_c['compound_id']} with 2D depiction.")

        # TEST D: Search Experimental Evidence Persistence E2E
        print("\n--- Test D: Search Experimental Evidence Persistence E2E ---")
        harvest_res = client.post(f"/api/compounds/{ibu_row_id}/experimental-harvest/preview", json={
            "confirm_public_identifier_search": True,
            "cas": "15687-27-1",
            "name": "Ibuprofen",
            "sources": ["PubChem PUG View", "ChEMBL"]
        })
        assert harvest_res.status_code == 200, f"Evidence harvest failed: {harvest_res.text}"
        harvest_data = harvest_res.json()
        records_harvested = harvest_data.get("records", [])
        print(f"Harvested {len(records_harvested)} external experimental records for Ibuprofen.")
        assert len(records_harvested) > 0, "No records returned from public evidence harvest"

        # Hydrate workspace initially
        ibu_vid = c_a["version"]["id"]
        ws_res_1 = client.get(f"/api/compound-versions/{ibu_vid}/workspace")
        assert ws_res_1.status_code == 200
        ws_1 = ws_res_1.json()
        ev_1 = ws_1.get("external_experimental_evidence", [])
        assert len(ev_1) > 0, "External evidence not present in workspace"
        print(f"Workspace initial hydration: {len(ev_1)} external evidence records loaded from DB.")

        # Reopen simulation ("Leave and Re-enter")
        ws_res_2 = client.get(f"/api/compound-versions/{ibu_vid}/workspace")
        assert ws_res_2.status_code == 200
        ws_2 = ws_res_2.json()
        ev_2 = ws_2.get("external_experimental_evidence", [])
        assert len(ev_2) == len(ev_1), f"Persistence failure: {len(ev_2)} != {len(ev_1)}"
        print(f"Test D passed: Reopen verified without re-search: exactly {len(ev_2)} records persisted.")

        # TEST E: Prediction Execution & Snapshot Persistence
        print("\n--- Test E: Prediction Execution & Snapshot Persistence E2E ---")
        pred_res = client.post(f"/api/compounds/{ibu_row_id}/predict-workflow", json={})
        assert pred_res.status_code == 202, f"Prediction workflow failed: {pred_res.text}"
        pred_data = pred_res.json()
        assert pred_data.get("status") in {"COMPLETE", "RUNNING"}
        print(f"Prediction run executed: status={pred_data.get('status')}")

        # Verify workspace snapshot persistence
        ws_pred_1 = client.get(f"/api/compound-versions/{ibu_vid}/workspace").json()
        persisted_runs = ws_pred_1.get("experimental_prediction_runs", []) or ws_pred_1.get("prediction_audit", [])
        assert len(persisted_runs) > 0, "Prediction audit / run missing from workspace"
        print(f"Prediction audit records in workspace: {len(persisted_runs)}")

        # Re-enter without rerun: verify reused_existing_run
        pred_rerun = client.post(f"/api/compounds/{ibu_row_id}/predict-workflow", json={"force_rerun": False})
        assert pred_rerun.status_code == 202
        assert pred_rerun.json().get("reused_existing_run") is True, "Reused existing run flag not set"
        print("Test E passed: Idempotent re-entry safely reused frozen prediction without silent recalculation.")

        # TEST F: DrugBank Reference Project Integrity
        print("\n--- Test F: DrugBank Reference Project (ID 300) Integrity ---")
        db_res = client.get("/api/projects/300")
        assert db_res.status_code == 200, f"Failed to retrieve DrugBank project: {db_res.text}"
        db_proj = db_res.json()
        assert db_proj["name"] == "DrugBank"
        assert "GLOBAL_MODEL_DEVELOPMENT" in db_proj["indication"]
        compounds_80 = db_proj.get("compounds", [])
        assert len(compounds_80) == 80, f"DrugBank compound count: {len(compounds_80)} != 80"

        # Check first and last compound
        c_first = compounds_80[0]
        c_last = compounds_80[-1]
        print(f"DrugBank first compound: {c_first['compound_id']} ({c_first['name']})")
        print(f"DrugBank 80th compound: {c_last['compound_id']} ({c_last['name']})")
        assert c_first["name"] and c_last["name"]

        # Detail workspace for first DrugBank drug
        db_first_detail = client.get(f"/api/compounds/{c_first['row_id']}?include_versions=true").json()
        assert db_first_detail["version"] is not None
        db_first_vid = db_first_detail["version"]["id"]
        db_first_ws = client.get(f"/api/compound-versions/{db_first_vid}/workspace").json()
        db_first_ev = db_first_ws.get("external_experimental_evidence", [])
        print(f"DrugBank first compound workspace: {len(db_first_ev)} external evidence records loaded.")

        # Total DrugBank evidence check in DB
        db_ev_cnt = db.scalar(text("""
            SELECT count(eee.id)
            FROM external_experimental_evidence eee
            JOIN compound_versions cv ON cv.id = eee.compound_version_id
            JOIN compounds c ON c.id = cv.compound_row_id
            WHERE c.project_id = 300
        """))
        assert db_ev_cnt == 577, f"DrugBank evidence count mismatch: {db_ev_cnt} != 577"
        print(f"DrugBank evidence count verified: {db_ev_cnt} records across 80 compounds.")

        # Total DrugBank qualification freezes check
        db_freeze_cnt = db.scalar(text("SELECT count(*) FROM qualification_prediction_freezes WHERE project_id = '300'"))
        assert db_freeze_cnt == 54, f"DrugBank qualification freeze count mismatch: {db_freeze_cnt} != 54"
        db_endpoint_cnt = db.scalar(text("SELECT count(distinct endpoint_id) FROM qualification_prediction_freezes WHERE project_id = '300'"))
        assert db_endpoint_cnt == 18, f"DrugBank unique endpoints mismatch: {db_endpoint_cnt} != 18"
        print(f"DrugBank qualification prediction freezes verified: {db_freeze_cnt} records across {db_endpoint_cnt} endpoints.")
        print("Test F passed: DrugBank 80 compounds, evidence (577), and qualification freezes (18) 100% intact.")

    finally:
        # CLEANUP: Delete the temporary test project
        if test_project_id:
            print("\n--- Automatic Cleanup of Test Fixture Project ---")
            _delete_project_tree_rows(db, [test_project_id])
            db.commit()
            print(f"Successfully deleted test fixture project {test_project_id}.")

        db.close()

    # Post-cleanup Integrity & Orphan Check
    print("\n--- Final Database Foreign-Key & Orphan Audit ---")
    conn = sqlite3.connect("drug_opt.db")
    cur = conn.cursor()

    cur.execute("PRAGMA foreign_key_check;")
    fk_violations = cur.fetchall()
    print(f"PRAGMA foreign_key_check: {len(fk_violations)} violations")
    assert len(fk_violations) == 0, f"FK violations: {fk_violations}"

    cur.execute("PRAGMA integrity_check;")
    integrity = cur.fetchall()
    print(f"PRAGMA integrity_check: {integrity}")
    assert integrity == [("ok",)]

    orphan_checks = [
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

    for label, query in orphan_checks:
        cur.execute(query)
        cnt = cur.fetchone()[0]
        print(f"Orphan count [{label}]: {cnt}")
        assert cnt == 0, f"Orphans found: {label} = {cnt}"

    # Verify protected projects
    cur.execute("SELECT id, name FROM projects ORDER BY id")
    remaining_projs = cur.fetchall()
    print(f"\nRemaining Production Projects ({len(remaining_projs)}):")
    for pid, name in remaining_projs:
        print(f" - ID {pid}: {name}")
    assert set(p[0] for p in remaining_projs) == {1, 3, 5, 300}, f"Mismatch: {remaining_projs}"

    conn.close()
    print("\n=================================================================")
    print(" ALL TESTS A-F PASSED! ZERO ORPHANS, ZERO CORRUPTION, 100% CLEAN ")
    print("=================================================================")

if __name__ == "__main__":
    run_verification()
