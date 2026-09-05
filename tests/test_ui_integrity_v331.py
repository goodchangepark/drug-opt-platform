import sqlite3
import pytest
from fastapi.testclient import TestClient
from backend.main import app
from backend.prediction_maturity import get_maturity_statistics, get_endpoint_maturity_registry

client = TestClient(app)

def test_endpoint_maturity_taxonomy_50():
    """Verify authoritative 50-endpoint maturity registry and statistics."""
    stats = get_maturity_statistics()
    assert stats["total_endpoints"] == 50
    assert stats["level_breakdown"]["level_1_base"] == 23
    assert stats["level_breakdown"]["level_2_validated_base"] == 15
    assert stats["level_breakdown"]["level_3_validated_multi_model"] == 3
    assert stats["level_breakdown"]["level_4_production_validated"] == 9
    assert stats["level_breakdown"]["level_5_mature"] == 0

    resp = client.get("/api/prediction-engine/endpoint-maturity")
    assert resp.status_code == 200
    data = resp.json()
    assert data["total_endpoints"] == 50
    assert len(data["endpoints"]) == 50

def test_prediction_engine_current_baseline():
    """Verify v3.3.1 baseline protection and policy hash."""
    resp = client.get("/api/prediction-engine/current")
    assert resp.status_code == 200
    data = resp.json()
    assert data["current_production_engine"]["engine_id"] == "drugopt-prediction-engine-v3@3.3.1"
    assert data["current_production_engine"]["engine_version"] == "3.3.1"
    assert data["current_production_engine"]["status"] == "PRODUCTION_DEFAULT"
    assert data["current_production_engine"]["policy_hash"] == "4647810a58bdbdbc700e4f5c26c5a187032e5cebc80bee6b0d64738f640954a9"
    assert data["endpoint_maturity"]["total_endpoints"] == 50
    assert data["endpoint_maturity"]["level_breakdown"]["level_4_production_validated"] == 9

def test_drugbank_150_cas_hydration():
    """Verify 150/150 DrugBank compounds have CAS numbers and canonical identifiers."""
    conn = sqlite3.connect("drug_opt.db")
    c = conn.cursor()
    c.execute("SELECT id, name, cas_number FROM compounds WHERE project_id = 300")
    db_rows = c.fetchall()
    assert len(db_rows) == 150
    missing_cas = [r for r in db_rows if not r[2] or r[2].strip() == ""]
    assert len(missing_cas) == 0, f"Found {len(missing_cas)} compounds missing CAS: {missing_cas[:5]}"

    # Check compound_identifiers for DrugBank compounds
    c.execute("""
        SELECT ci.identifier_type, count(*)
        FROM compound_identifiers ci
        JOIN compounds c ON ci.compound_id = c.id
        WHERE c.project_id = 300
        GROUP BY ci.identifier_type
    """)
    id_counts = dict(c.fetchall())
    assert id_counts.get("CAS") == 150
    assert id_counts.get("DRUGBANK_ID") == 150
    assert id_counts.get("CHEMBL_ID") == 150
    assert id_counts.get("PUBCHEM_CID") == 150
    assert id_counts.get("UNII") == 150
    conn.close()

def test_historical_prediction_runs_protected():
    """Verify historical prediction runs 1-128 were not modified or deleted."""
    conn = sqlite3.connect("drug_opt.db")
    c = conn.cursor()
    c.execute("SELECT count(*) FROM prediction_runs WHERE id <= 128")
    hist_count = c.fetchone()[0]
    assert hist_count == 102, f"Expected 102 historical runs up to ID 128, found {hist_count}"
    conn.close()

def test_all_compounds_in_projects_1_3_5_have_v331_runs():
    """Verify all 15 compounds in projects 1, 3, and 5 have valid v3.3.1 prediction runs."""
    conn = sqlite3.connect("drug_opt.db")
    c = conn.cursor()
    c.execute("""
        SELECT c.project_id, c.id, c.compound_id, cv.id as version_id
        FROM compounds c
        JOIN compound_versions cv ON cv.compound_row_id = c.id AND cv.version_number = c.current_version
        WHERE c.project_id IN (1, 3, 5)
        ORDER BY c.project_id, c.id
    """)
    active_compounds = c.fetchall()
    assert len(active_compounds) == 15

    for pid, cid, clabel, vid in active_compounds:
        c.execute("""
            SELECT id, model_version, stage
            FROM prediction_runs
            WHERE version_id = ? AND model_version = '3.3.1'
        """, (vid,))
        run = c.fetchone()
        assert run is not None, f"Compound {cid} ({clabel}) in project {pid} missing v3.3.1 prediction run"
        assert run[2] is not None
    conn.close()

def test_compound_workspace_endpoint_comparison_maturity():
    """Verify that compound workspace returns endpoint comparison with proper maturity."""
    resp = client.get("/api/compound-versions/1/workspace")
    assert resp.status_code == 200
    data = resp.json()
    assert "endpoint_comparison" in data
    ep_comp = data["endpoint_comparison"]
    assert "scientific_rows" in ep_comp
    rows = ep_comp["scientific_rows"]
    assert len(rows) > 0

    sol_row = next((r for r in rows if r["canonical_endpoint"] == "SOLUBILITY_GENERIC"), None)
    assert sol_row is not None
    assert sol_row.get("prediction") is not None
    assert sol_row["prediction"]["maturity"]["level"] == 4
    assert sol_row["prediction"]["maturity"]["label"] == "Production Validated"

    caco2_row = next((r for r in rows if r["canonical_endpoint"] == "CACO2_PAPP_AB"), None)
    assert caco2_row is not None
    assert caco2_row.get("prediction") is not None
    assert caco2_row["prediction"]["maturity"]["level"] == 4
    assert caco2_row["prediction"]["maturity"]["label"] == "Production Validated"
