"""
Stage 5B-4 Refinement 3: Dashboard Redesign & Compound Save Workflow Restoration Tests

Covers:
- Application version 0.6.3-stage5b4-ui and version history registry
- API endpoints: /api/health, /api/dashboard, /api/help/registry
- Compound save contract (validation, standardization, CompoundVersion creation, property calculation)
- Compound persistence and project isolation across session reloads
- Robust error handling for empty names, invalid structures, duplicate IDs
- Dashboard capability registry sync and layout structure
- Static assets verification (Noto Sans KR, #f4f7fb background, status badge styles, Ketcher fallback)
- Database schema integrity (PRAGMA integrity_check, PRAGMA foreign_key_check)
"""

import json
import os
import pytest
from pathlib import Path
from fastapi.testclient import TestClient
from sqlalchemy import select, text

from backend.main import app
from backend.database import SessionLocal, engine
from backend.models import Project, Compound, CompoundVersion, PropertyCalculation
from backend.platform_info import APP_VERSION, CURRENT_STAGE, version_history, latest_release_date
from backend.capabilities import build_capability_summary

ROOT = Path(__file__).resolve().parents[1]
client = TestClient(app)


def test_app_version_0_6_3():
    assert APP_VERSION == "1.0.0"
    assert CURRENT_STAGE == "5B-4"
    vh = version_history()
    assert len(vh) >= 13
    v063 = next(entry for entry in vh if entry["version"] == "0.6.3-stage5b4-ui")
    assert v063["stage"] == "Stage 5B-4 Refinement 3"
    assert "Dashboard Redesign & Compound Save" in v063["milestone"]
    assert vh[-1]["version"] == "v4.2"


def test_api_health():
    response = client.get("/api/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["version"] == APP_VERSION
    assert data["step"] == "5B-4"


def test_api_help_registry():
    response = client.get("/api/help/registry")
    assert response.status_code == 200
    data = response.json()
    assert data["application"]["version"] == APP_VERSION
    assert "capability_summary" in data
    assert "pk_method_registry" in data


def _clean_project_by_name(name: str):
    with SessionLocal() as db:
        projs = db.scalars(select(Project).where(Project.name == name)).all()
        for p in projs:
            client.request("DELETE", f"/api/projects/{p.id}", json={"confirmation_name": p.name})


def test_compound_save_workflow_with_immediate_properties():
    _clean_project_by_name("__TEST_SAVE_WORKFLOW_PROJ__")
    # 1. Create a clean test project
    create_proj = client.post("/api/projects", json={
        "name": "__TEST_SAVE_WORKFLOW_PROJ__",
        "target": "EGFR",
        "molecule_type": "Small Molecule",
        "description": "Compound save validation"
    })
    assert create_proj.status_code == 201
    proj_id = create_proj.json()["id"]

    try:
        # 2. Save a valid small molecule compound (Gefitinib SMILES)
        gefitinib_smiles = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
        save_resp = client.post(f"/api/projects/{proj_id}/compounds", json={
            "name": "Gefitinib-Lead",
            "compound_id": "GEF-001",
            "smiles": gefitinib_smiles,
            "notes": "Primary lead candidate",
            "calculate": True
        })
        assert save_resp.status_code == 201
        saved = save_resp.json()
        assert saved["name"] == "Gefitinib-Lead"
        assert saved["compound_id"] == "GEF-001"
        assert saved["version"]["canonical_smiles"] is not None
        assert (saved["version"]["properties"].get("formula") or saved["version"]["properties"].get("molecular_formula")) == "C22H24ClFN4O3"
        assert round(saved["version"]["properties"]["molecular_weight"], 1) == 446.9
        assert saved["version"]["properties"] is not None
        assert saved["status"] in ("CALCULATED", "STRUCTURE_READY")

        compound_row_id = saved["id"] if "id" in saved else saved.get("row_id")

        # 3. Verify database persistence in fresh session
        with SessionLocal() as db:
            c = db.get(Compound, compound_row_id)
            assert c is not None
            assert c.project_id == proj_id
            assert c.name == "Gefitinib-Lead"
            assert c.current_version == 1

            versions = db.scalars(select(CompoundVersion).where(CompoundVersion.compound_row_id == compound_row_id)).all()
            assert len(versions) == 1
            v = versions[0]
            assert v.canonical_smiles == saved["version"]["canonical_smiles"]
            assert v.calculation_json["provenance"]["type"] == "Calculated"
            assert "engine" in v.calculation_json["provenance"]
            assert (v.properties_json.get("formula") or v.properties_json.get("molecular_formula")) == "C22H24ClFN4O3"
            assert v.properties_json["molecular_weight"] > 400

            props = db.scalars(select(PropertyCalculation).where(PropertyCalculation.version_id == v.id)).all()
            assert len(props) > 0

        # 4. Verify fetch compound detail via API
        detail_resp = client.get(f"/api/compounds/{compound_row_id}?include_versions=true")
        assert detail_resp.status_code == 200
        detail = detail_resp.json()
        assert detail["name"] == "Gefitinib-Lead"
        assert detail["version"]["canonical_smiles"] is not None

    finally:
        # Cleanup test project
        client.request("DELETE", f"/api/projects/{proj_id}", json={"confirmation_name": "__TEST_SAVE_WORKFLOW_PROJ__"})


def test_compound_save_error_handling():
    _clean_project_by_name("__TEST_SAVE_ERRORS_PROJ__")
    create_proj = client.post("/api/projects", json={
        "name": "__TEST_SAVE_ERRORS_PROJ__",
        "target": "BRAF",
        "molecule_type": "Small Molecule",
        "description": "Error handling validation"
    })
    assert create_proj.status_code == 201
    proj_id = create_proj.json()["id"]

    try:
        # Case 1: Empty compound name
        err1 = client.post(f"/api/projects/{proj_id}/compounds", json={
            "name": "   ",
            "compound_id": "",
            "smiles": "c1ccccc1",
            "calculate": True
        })
        assert err1.status_code == 400
        assert "Compound name is required" in err1.json()["detail"]

        # Case 2: Invalid chemical SMILES
        err2 = client.post(f"/api/projects/{proj_id}/compounds", json={
            "name": "Invalid-Mol",
            "compound_id": "INV-01",
            "smiles": "INVALID_NOT_A_SMILES_STRING",
            "calculate": True
        })
        assert err2.status_code == 400
        assert "Structure could not be standardized" in err2.json()["detail"]

        # Case 3: Duplicate compound ID in same project
        ok1 = client.post(f"/api/projects/{proj_id}/compounds", json={
            "name": "Valid-Mol-1",
            "compound_id": "MOL-DUP-01",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "calculate": True
        })
        assert ok1.status_code == 201

        err3 = client.post(f"/api/projects/{proj_id}/compounds", json={
            "name": "Valid-Mol-2",
            "compound_id": "MOL-DUP-01",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "calculate": True
        })
        assert err3.status_code == 409
        assert "Compound ID already exists in project" in err3.json()["detail"]

    finally:
        client.request("DELETE", f"/api/projects/{proj_id}", json={"confirmation_name": "__TEST_SAVE_ERRORS_PROJ__"})


def test_dashboard_redesign_capability_structure():
    response = client.get("/api/dashboard")
    assert response.status_code == 200
    data = response.json()

    assert "totals" in data
    assert "projects" in data["totals"]
    assert "compounds" in data["totals"]
    assert "capability_summary" in data

    caps = data["capability_summary"]
    assert caps["stage"] == "5B-4"
    assert len(caps["groups"]) == 7

    group_titles = [g["title"] for g in caps["groups"]]
    expected_titles = [
        "Structure & Chemistry",
        "Activity & SAR",
        "ADME",
        "CYP & Transporters",
        "Safety / Toxicology",
        "Optimization",
        "PK / DMPK"
    ]
    assert group_titles == expected_titles

    # Verify model registry live count
    models = data["model_registry"]
    operational_models = [m for m in models if m.get("status") == "READY" or m.get("availability") == "READY"]
    assert len(operational_models) >= 15


def test_static_files_standards():
    app_js = (ROOT / "frontend/static/app.js").read_text(encoding="utf-8")
    app_css = (ROOT / "frontend/static/app.css").read_text(encoding="utf-8")

    # Typography & Palette
    assert "Noto Sans KR" in app_css
    assert "--bg: #f4f7fb;" in app_css
    assert "grid-template-columns: 290px minmax(0, 1fr);" in app_css

    # Dashboard CSS Classes
    assert ".dashboard-hero" in app_css
    assert ".platform-capabilities-line" in app_css
    assert ".dashboard-stats-grid" in app_css
    assert ".scientific-workspace-section" in app_css
    assert ".scientific-workspace-grid" in app_css
    assert ".scientific-card" in app_css
    assert ".scientific-card-row" in app_css

    # Save logic & Ketcher fallback in app.js
    assert "getSmiles" in app_js
    assert "savingCompound" in app_js
    assert "calculate:true" in app_js.replace(" ", "")

    # Version in app.js
    assert "renderPredictionMaturity" in app_js


def test_database_integrity():
    with engine.connect() as conn:
        integrity = conn.execute(text("PRAGMA integrity_check;")).fetchall()
        assert len(integrity) == 1
        assert integrity[0][0] == "ok"

        fk_check = conn.execute(text("PRAGMA foreign_key_check;")).fetchall()
        assert len(fk_check) == 0
