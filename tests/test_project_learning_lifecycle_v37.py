import uuid

from fastapi.testclient import TestClient
from sqlalchemy import inspect

from backend.database import engine
from backend.main import app


def _name():
    return "v37 lifecycle test " + uuid.uuid4().hex[:10]


def test_runtime_compounds_schema_keeps_cas_nullable():
    column = next(row for row in inspect(engine).get_columns("compounds") if row["name"] == "cas_number")
    assert column["nullable"] is True


def test_experiment_arrival_creates_auditable_pair_or_exclusion_record():
    with TestClient(app) as client:
        project = client.post("/api/projects", json={"name": _name()}).json()
        project_id = project["id"]
        try:
            compound = client.post(
                f"/api/projects/{project_id}/compounds",
                json={"compound_id": "LIFE-1", "name": "Lifecycle test", "smiles": "CCO", "cas_number": ""},
            ).json()
            version_id = compound["version"]["id"]
            response = client.post(
                f"/api/projects/{project_id}/admet/measurements",
                json={"version_id": version_id, "endpoint": "Solubility", "value": -2.0, "unit": "log10(mol/L)"},
            )
            assert response.status_code == 201
            ledger = client.get(f"/api/compound-versions/{version_id}/learning-ledger")
            assert ledger.status_code == 200
            row = ledger.json()["ledger"][0]
            assert row["endpoint"] == "Solubility"
            assert row["pair_class"] == "HISTORICAL_VISIBLE"
            assert row["adaptation_eligibility"] is False
            assert row["exclusion_reason"] == "NO_PREEXPERIMENTAL_FREEZE"
            project_ledger = client.get(f"/api/projects/{project_id}/learning-ledger")
            assert project_ledger.status_code == 200
            assert project_ledger.json()["ledger"][0]["compound_version_id"] == version_id
        finally:
            client.delete(f"/api/projects/{project_id}")


def test_prediction_freeze_precedes_experiment_and_preserves_base_only_state():
    with TestClient(app) as client:
        project = client.post("/api/projects", json={"name": _name()}).json()
        project_id = project["id"]
        try:
            compound = client.post(
                f"/api/projects/{project_id}/compounds",
                json={"compound_id": "LIFE-2", "name": "Prospective lifecycle test", "smiles": "CCO", "cas_number": ""},
            ).json()
            version_id = compound["version"]["id"]
            prediction = client.post(f"/api/admet/predict/{version_id}")
            assert prediction.status_code == 202
            solubility = next(row for row in prediction.json()["predictions"] if row["endpoint"] == "Solubility")
            assert solubility["prediction_snapshot"]["project_prediction"] is None
            assert solubility["prediction_snapshot"]["experiment_known_at_prediction_time"] is False
            response = client.post(
                f"/api/projects/{project_id}/admet/measurements",
                json={"version_id": version_id, "endpoint": "Solubility", "value": -2.0, "unit": solubility["unit"]},
            )
            assert response.status_code == 201
            row = client.get(f"/api/compound-versions/{version_id}/learning-ledger").json()["ledger"][0]
            assert row["pair_class"] == "TRUE_PROSPECTIVE"
            assert row["adaptation_eligibility"] is True
            assert row["project_prediction"] is None
            assert row["project_absolute_error"] is None
        finally:
            client.delete(f"/api/projects/{project_id}")


def test_adapter_rollback_is_explicit_and_preserves_history_endpoint():
    with TestClient(app) as client:
        project = client.post("/api/projects", json={"name": _name()}).json()
        project_id = project["id"]
        try:
            result = client.post(f"/api/projects/{project_id}/project-adaptation/Solubility/deactivate")
            assert result.status_code == 200
            assert result.json()["status"] == "BASE_ONLY"
            assert result.json()["adapter_history_preserved"] is True
        finally:
            client.delete(f"/api/projects/{project_id}")
