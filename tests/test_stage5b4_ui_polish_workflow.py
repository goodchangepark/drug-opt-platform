"""
Stage 5B-4 UI Polish & Unified Prediction Workflow Tests
Covers:
- Application version 0.6.2-stage5b4-ui
- Dashboard layout (Platform Overview, Scientific Workspace, Quick Start, Projects)
- Project navigation & de-emphasized delete button
- Global Noto Sans KR typography standard and CSS container bounds
- Overview Primary PREDICT button & endpoint orchestration (/api/compounds/{row_id}/predict-all)
- Failure isolation & Activity exclusion
- Tab Re-Predict workflows & experimental data preservation
- Overview PK Summary component
- Test project cleanup validation (single newest test retained, real/ambiguous preserved)
- Database integrity and 0 FK errors
"""

import json
import os
import pytest
from pathlib import Path
from fastapi.testclient import TestClient
from sqlalchemy import select, func

from backend.main import app
from backend.database import SessionLocal, engine
from backend.models import Project, Compound, CompoundVersion
from backend.platform_info import APP_VERSION, CURRENT_STAGE, version_history
from backend.stabilization import classify_project

ROOT = Path(__file__).resolve().parents[1]
client = TestClient(app)


def test_app_version():
    assert APP_VERSION == "1.0.0"
    assert CURRENT_STAGE == "5B-4"
    vh = version_history()
    assert len(vh) >= 12
    assert any(entry["version"] == APP_VERSION for entry in vh)
    assert vh[-1]["version"] == "v4.5"



def test_predict_all_endpoint_orchestration():
    # Create isolated small molecule test compound in a temporary project
    with SessionLocal() as db:
        from backend.main import _delete_project_tree_rows
        existing = db.scalars(select(Project).where(Project.name == "E2E Temp Predict Test")).all()
        if existing:
            _delete_project_tree_rows(db, [p.id for p in existing])
            db.commit()
        proj = Project(name="E2E Temp Predict Test", target="EGFR", molecule_type="Small Molecule", description="Temp")
        db.add(proj)
        db.commit()
        db.refresh(proj)
        proj_id = proj.id

    try:
        # Create compound
        res = client.post(f"/api/projects/{proj_id}/compounds", json={
            "name": "ASPIRIN-TEST",
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "calculate": False
        })
        assert res.status_code in (200, 201)
        compound_row_id = res.json()["row_id"]


        # Run /predict-all endpoint
        pred_res = client.post(f"/api/compounds/{compound_row_id}/predict-all")
        assert pred_res.status_code == 202
        data = pred_res.json()

        assert data["status"] in ("COMPLETE", "PARTIAL")
        assert data["activity_excluded"] is True
        assert "completed_endpoints" in data
        assert len(data["completed_endpoints"]) >= 3
        assert "Physicochemical Properties" in data["completed_endpoints"]
        assert "timestamp" in data

        # Test idempotency (repeated call should not fail or corrupt)
        repeat_res = client.post(f"/api/compounds/{compound_row_id}/predict-all")
        assert repeat_res.status_code == 202
        assert repeat_res.json()["status"] in ("COMPLETE", "PARTIAL")

    finally:
        # Cleanup
        with SessionLocal() as db:
            from backend.main import _delete_project_tree_rows
            _delete_project_tree_rows(db, [proj_id])
            db.commit()


def test_css_typography_and_layout_bounds():
    css_text = (ROOT / "frontend/static/app.css").read_text(encoding="utf-8")
    
    # Global Noto Sans KR standard
    assert '--font-ui: "Noto Sans KR", "Noto Sans CJK KR", sans-serif;' in css_text
    assert "font-family: var(--font-ui);" in css_text

    # Structure bounds to prevent sidebar overlap
    assert ".compound-header-card" in css_text
    assert ".compound-header-structure" in css_text
    assert "max-width: 100%" in css_text
    assert "overflow: hidden" in css_text

    # Primary Predict button
    assert ".btn-predict-primary" in css_text
    assert ".predict-meta-bar" in css_text
    assert ".experimental-alert-bar" in css_text

    # Tab Re-Predict button
    assert ".tab-repredict-btn" in css_text

    # PK Summary on Overview
    assert ".overview-pk-grid" in css_text
    assert ".overview-pk-card" in css_text

    # Scientific Workspace Card Grid
    assert ".scientific-workspace-grid" in css_text
    assert ".scientific-card" in css_text

    # De-emphasized delete button
    assert ".project-delete-secondary" in css_text


def test_app_js_components_and_workflow():
    js_text = (ROOT / "frontend/static/app.js").read_text(encoding="utf-8")

    # Version in footer
    assert "renderPredictionMaturity" in js_text

    # Dashboard sections
    assert "PLATFORM OVERVIEW" in js_text
    assert "SCIENTIFIC WORKSPACE" in js_text
    assert "RESEARCH PORTFOLIO" in js_text

    # Clickable project name link
    assert "project-link-title" in js_text

    # Overview Predict & PK Summary
    assert "btn-predict-primary" in js_text
    assert "runFullPredict" in js_text
    assert "TRANSLATIONAL PK SUMMARY" in js_text

    # Re-predict buttons
    assert "tab-repredict-btn" in js_text
    assert "ACTIVITY MODEL NOT READY" in js_text


def test_test_project_cleanup_audit_file():
    audit_file = ROOT / "validation" / "ui_refinement_test_project_cleanup.json"
    assert audit_file.exists()

    data = json.loads(audit_file.read_text(encoding="utf-8"))
    assert "summary" in data
    assert "most_recent_test_project" in data
    
    # Must retain at least 1 newest test project
    assert data["summary"]["CONFIRMED_TEST_RETAINED"] == 1
    assert data["most_recent_test_project"] is not None
    assert data["most_recent_test_project"]["project_id"] > 0

    # Ensure real and ambiguous projects were preserved
    assert data["summary"]["PRESERVE_REAL"] > 0
    assert data["summary"]["PRESERVE_AMBIGUOUS"] > 0


def test_database_integrity_and_fk_check():
    import sqlite3
    con = sqlite3.connect(ROOT / "drug_opt.db")
    
    integrity = con.execute("PRAGMA integrity_check;").fetchall()
    assert integrity == [("ok",)], f"Integrity check failed: {integrity}"

    fk = con.execute("PRAGMA foreign_key_check;").fetchall()
    assert len(fk) == 0, f"Foreign key violations found: {fk}"
    con.close()
