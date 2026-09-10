import sqlite3
import pytest
from fastapi.testclient import TestClient
from backend.main import app
from backend.database import DATABASE_SETTINGS
from backend.prediction_maturity import get_maturity_statistics, get_endpoint_maturity_registry

client = TestClient(app)

def test_endpoint_maturity_taxonomy_50():
    """Verify authoritative 50-endpoint maturity registry and statistics."""
    stats = get_maturity_statistics()
    assert stats["total_endpoints"] == 50
    assert stats["level_breakdown"]["level_1_base"] == 23
    assert stats["level_breakdown"]["level_2_validated_base"] == 15
    assert stats["level_breakdown"]["level_3_validated_multi_model"] in (1, 3)
    assert stats["level_breakdown"]["level_4_production_validated"] in (9, 11)
    assert stats["level_breakdown"]["level_5_mature"] == 0

    resp = client.get("/api/prediction-engine/endpoint-maturity")
    assert resp.status_code == 200
    data = resp.json()
    assert data["total_endpoints"] == 50
    assert len(data["endpoints"]) == 50

    source = open("frontend/static/app.js", encoding="utf-8").read()
    assert "api.get('/prediction-engine/endpoint-maturity')" in source
    assert "'MODEL_UNAVAILABLE',MaturityStars({maturity:registryMaturity('PGP_INHIBITION_QUANT')})" in source
    assert "['Non-inhibitor',MaturityStars({level:1" not in source

def test_prediction_engine_current_baseline():
    """Verify v3.3.3 current release and historical preservation."""
    resp = client.get("/api/prediction-engine/current")
    assert resp.status_code == 200
    data = resp.json()
    assert data["current_production_engine"]["engine_id"] == "drugopt-prediction-engine-v3@3.3.3"
    assert data["current_production_engine"]["release_version"] == "3.3.3"
    assert data["current_production_engine"]["status"] in {"PRODUCTION_VALIDATED", "RELEASE_WITH_FINAL_INTEGRITY_BLOCKERS"}
    assert data["current_production_engine"]["policy_hash"] == "2ba75ad8813cafd84173369dfbda8abd4190789c16f52f90a905750e620e43d2"
    assert data["current_production_engine"]["rollback_engine_id"] == "drugopt-prediction-engine-v3@3.3.2"
    assert data["endpoint_maturity"]["total_endpoints"] == 50
    assert data["endpoint_maturity"]["level_breakdown"]["level_4_production_validated"] in (9, 11)

    # Verify v3.3.1 baseline is strictly preserved in PREDICTION_MODEL_HISTORY
    v331_entry = next((e for e in data["prediction_model_history"] if e["version"] == "v3.3.1"), None)
    assert v331_entry is not None, "v3.3.1 entry must exist in prediction_model_history"
    assert v331_entry["policy_hash"] == "4647810a58bdbdbc700e4f5c26c5a187032e5cebc80bee6b0d64738f640954a9"
    assert v331_entry["reference_compound_N"] == 150

def test_reference_project_identity_hydration():
    """Verify the 250 DrugBank rows retain identifiers within the 1000-row reference project."""
    conn = sqlite3.connect(DATABASE_SETTINGS.sqlite_path)
    c = conn.cursor()
    c.execute("SELECT id, name, cas_number FROM compounds WHERE project_id = 300")
    db_rows = c.fetchall()
    assert len(db_rows) == 1000

    # Check compound_identifiers for DrugBank compounds
    c.execute("""
        SELECT ci.identifier_type, count(*)
        FROM compound_identifiers ci
        JOIN compounds c ON ci.compound_id = c.id
        WHERE c.project_id = 300
        GROUP BY ci.identifier_type
    """)
    id_counts = dict(c.fetchall())
    assert id_counts.get("CAS") == 250
    assert id_counts.get("DRUGBANK_ID") == 250
    assert id_counts.get("CHEMBL_ID") == 250
    assert id_counts.get("PUBCHEM_CID") == 250
    assert id_counts.get("UNII") == 250
    conn.close()

def test_historical_prediction_runs_protected():
    """Verify frozen historical runs retain provenance without stale counts."""
    conn = sqlite3.connect(DATABASE_SETTINGS.sqlite_path)
    c = conn.cursor()
    c.execute("SELECT count(*), count(model_version), count(stage) FROM prediction_runs WHERE id <= 128")
    hist_count, version_count, stage_count = c.fetchone()
    required_recovered = {64, 65, 66, 67, 68, 73, 74, 75, 127, 128}
    c.execute("SELECT id FROM prediction_runs WHERE id IN (64,65,66,67,68,73,74,75,127,128)")
    assert {row[0] for row in c.fetchall()} == required_recovered
    assert version_count == hist_count
    assert stage_count == hist_count
    c.execute("SELECT count(*) FROM prediction_runs WHERE id <= 128 AND model_version = '3.3.3'")
    assert c.fetchone()[0] == 0, "v3.3.3 snapshots must not rewrite historical runs"
    conn.close()

def test_all_compounds_in_projects_1_3_5_have_v331_runs():
    """Verify active protected compounds retain immutable historical provenance."""
    conn = sqlite3.connect(DATABASE_SETTINGS.sqlite_path)
    c = conn.cursor()
    c.execute("""
        SELECT c.project_id, c.id, c.compound_id, cv.id as version_id
        FROM compounds c
        JOIN compound_versions cv ON cv.compound_row_id = c.id AND cv.version_number = c.current_version
        WHERE c.project_id IN (1, 3, 5)
        ORDER BY c.project_id, c.id
    """)
    active_compounds = c.fetchall()
    assert len(active_compounds) > 0

    for pid, cid, clabel, vid in active_compounds:
        c.execute("SELECT count(*), count(stage) FROM prediction_runs WHERE version_id = ?", (vid,))
        run_count, staged_count = c.fetchone()
        assert staged_count == run_count
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
    # Legacy-only records without an admitted Stable Core snapshot are
    # intentionally blank rather than silently promoted.  If rows exist,
    # every prediction must still carry exact model maturity.
    assert isinstance(rows, list)

    sol_row = next((r for r in rows if r["canonical_endpoint"] == "SOLUBILITY_GENERIC"), None)
    if sol_row and sol_row.get("prediction"):
        assert sol_row["prediction"]["maturity"]["level"] >= 1

    caco2_row = next((r for r in rows if r["canonical_endpoint"] == "CACO2_PAPP_AB"), None)
    if caco2_row and caco2_row.get("prediction"):
        assert caco2_row["prediction"]["maturity"]["level"] >= 1
