"""Regression contracts for Predict -> Stable Core -> single/compare display."""

from __future__ import annotations

from fastapi.testclient import TestClient
from sqlalchemy import select

from backend.admet import ADMETModelRegistry, ADMETPrediction
from backend.database import SessionLocal
from backend.main import app
from backend.models import Compound, CompoundVersion
from backend.stable_core import CurrentPredictionSnapshot, admit_current_prediction


def _two_cached_rlm_versions(db):
    rows = list(db.execute(
        select(CompoundVersion, Compound)
        .join(Compound, Compound.id == CompoundVersion.compound_row_id)
        .join(ADMETPrediction, ADMETPrediction.version_id == CompoundVersion.id)
        .join(ADMETModelRegistry, ADMETModelRegistry.id == ADMETPrediction.model_id)
        .where(
            ADMETModelRegistry.endpoint_name == "RLM intrinsic clearance",
            ADMETModelRegistry.model_name == "OpenADMET CheMeleon RLM intrinsic clearance",
            ADMETPrediction.execution_status == "SUCCESS",
            Compound.status != "ARCHIVED",
        )
        .order_by(Compound.project_id, Compound.id)
    ))
    by_project = {}
    for version, compound in rows:
        by_project.setdefault(compound.project_id, {})[compound.id] = (version, compound)
    return next(list(values.values())[:2] for values in by_project.values() if len(values) >= 2)


def test_predict_publishes_admitted_snapshot_and_compare_uses_identical_authority():
    with SessionLocal() as db:
        (first_version, first_compound), (second_version, second_compound) = _two_cached_rlm_versions(db)
        project_id = first_compound.project_id

    client = TestClient(app)
    first_result = client.post(f"/api/admet/predict/{first_version.id}")
    second_result = client.post(f"/api/admet/predict/{second_version.id}")
    assert first_result.status_code == second_result.status_code == 202
    first_publication = first_result.json()["current_publication"]
    published = next(row for row in first_publication if row["endpoint"] == "RLM_CLINT")
    assert published["status"] == "CALCULATED_AND_PUBLISHED"

    single = client.get(f"/api/compound-versions/{first_version.id}/scientific-tabs/admet").json()
    single_row = next(row for row in single["rows"] if row["canonical_endpoint"] == "RLM_CLINT")
    assert single_row["prediction"]["snapshot_id"] == published["snapshot_id"]
    assert single_row["prediction"]["engine_version"] == "drugopt-prediction-engine-v3@3.3.3"

    reloaded = client.get(f"/api/compound-versions/{first_version.id}/scientific-tabs/admet").json()
    reloaded_row = next(row for row in reloaded["rows"] if row["canonical_endpoint"] == "RLM_CLINT")
    assert reloaded_row["prediction"] == single_row["prediction"]

    comparison = client.get(
        f"/api/projects/{project_id}/compare?ids={first_compound.id},{second_compound.id}"
    ).json()
    compared = next(row for row in comparison["compounds"] if row["row_id"] == first_compound.id)
    assert compared["prediction_snapshot_ids"]["RLM"] == published["snapshot_id"]
    assert compared["prediction_metadata"]["RLM"]["value"] == single_row["prediction"]["value"]
    assert compared["prediction_metadata"]["RLM"]["unit"] == single_row["prediction"]["unit"]
    assert compared["prediction_metadata"]["RLM"]["species"] == single_row["species"]
    assert compared["prediction_metadata"]["RLM"]["model_id"] == single_row["prediction"]["model_id"]
    assert compared["prediction_metadata"]["RLM"]["model_version"] == single_row["prediction"]["model_version"]
    assert compared["prediction_metadata"]["RLM"]["engine_version"] == single_row["prediction"]["engine_version"]
    assert compared["prediction_metadata"]["RLM"]["mode"] == single_row["prediction"]["mode"]
    assert compared["prediction_metadata"]["RLM"]["context"] == single_row["context"]
    assert compared["prediction_metadata"]["RLM"]["maturity"] == single_row["maturity"]

    with SessionLocal() as db:
        snapshot = db.get(CurrentPredictionSnapshot, published["snapshot_id"])
        assert admit_current_prediction(db, snapshot).eligible is True


def test_calculated_unrouted_result_is_explicitly_not_eligible():
    with SessionLocal() as db:
        version, _compound = _two_cached_rlm_versions(db)[0]
    payload = TestClient(app).post(f"/api/admet/predict/{version.id}").json()
    solubility = next(row for row in payload["current_publication"] if row["endpoint"] == "SOLUBILITY_GENERIC")
    assert solubility == {
        "endpoint": "SOLUBILITY_GENERIC",
        "status": "CALCULATED_BUT_NOT_ELIGIBLE",
        "reason": "ROUTED_IMPLEMENTATION_MISMATCH",
    }
